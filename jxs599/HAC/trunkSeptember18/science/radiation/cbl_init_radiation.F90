!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Computes radiation absorbed by canopy and soil surface
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b
!
!
! ==============================================================================

MODULE cbl_init_radiation_module

  USE cable_data_module, ONLY : irad_type, point2constants

  IMPLICIT NONE

  PUBLIC init_radiation
  PRIVATE

  TYPE ( irad_type ) :: C


CONTAINS

  SUBROUTINE init_radiation( met, rad, veg, canopy, &
mp,                    &  
nrb,                   &
Clai_thresh,           &
Ccoszen_tols,          &
jls_standalone,        &
jls_radiation ,        &
veg_mask,              &
sunlit_mask,           &
sunlit_veg_mask,       &
reducedLAIdue2snow,    &
coszen,                &
ExtCoeff_beam,         &
ExtCoeff_dif,          &
EffExtCoeff_beam,      &
EffExtCoeff_dif,       &
VegXfang,              &
VegTaul,               &
VegRefl,               &
c1,                    &
rhoch,                 &
metDoY,                &
SW_down,               &
RadFbeam,              &
xk,                    &
CGauss_w,              &
Cpi,                   &
Cpi180,                &
subr_name              &
                         )
 
    ! Alternate version of init_radiation that uses only the
    ! zenith angle instead of the fluxes. This means it can be called
    ! before the cable albedo calculation.
    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
         veg_parameter_type!!, nrb, mp
    USE cable_common_module

USE cbl_spitter_module, ONLY : spitter
USE cbl_rhoch_module, ONLY : calc_rhoch

    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (met_type),       INTENT(INOUT) :: met
    TYPE (canopy_type),    INTENT(IN)    :: canopy
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg
    INTEGER :: ictr

!re-decl input args
integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
real :: Clai_thresh             !threshold LAI below which considered UN-vegetated
real :: Ccoszen_tols            !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)
real :: Cpi                     !PI - from cable_math_constants originally
real :: Cpi180                  !PI in radians - from cable_math_constants originally
LOGICAL :: jls_standalone       !runtime switch defined in cable_*main routines signifying this is jules_standalone
LOGICAL :: jls_radiation        !runtime switch defined in cable_*main routines signifying this is the radiation pathway 
!masks
logical :: veg_mask(mp)         !vegetated mask [formed by comparrisson of LAI CLAI_thresh ]
logical :: sunlit_mask(mp)      !sunlit mask [formed by comparrisson of coszen to coszen_tols i.e. is the sun up]
logical :: sunlit_veg_mask(mp)  !combined mask - BOTH sunlit and vegetated

REAL :: reducedLAIdue2snow(mp)         !Effective LAI given (potential sno coverage)
REAL :: coszen(mp)              ! cosine zenith angle of sun

REAL :: ExtCoeff_beam(mp)       !"raw" Extinction co-efficient for Direct Beam component of SW radiation
REAL :: ExtCoeff_dif(mp)        !"raw"Extinction co-efficient for Diffuse component of SW radiation
REAL :: EffExtCoeff_beam(mp,nrb)!Effective Extinction co-efficient for Direct Beam component of SW radiation
REAL :: EffExtCoeff_dif(mp,nrb) !Effective Extinction co-efficient for Diffuse component of SW radiation

integer :: metDoY(mp)           !Day of the Year [formerly met%doy]
REAL :: SW_down(mp,nrb)         !Downward SW radiation [formerly met%fsd]
REAL :: RadFbeam(mp,nrb)        !Beam Fraction of Downward SW radiation [formerly rad%fbeam]

!vaegetation parameters input via namelist
REAL :: VegXfang(mp)
REAL :: VegTaul(mp,nrb)
REAL :: VegRefl(mp,nrb)

!co-efficients used throughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)              ! extinct. coef.for beam rad. and black leaves

character(len=*) :: subr_name !where am i called from

!local_vars - common scaling co-efficients used throughout init_radiation
REAL :: xvlai2(mp,nrb) ! 2D vlai
REAL :: xphi1(mp)      ! leaf angle parmameter 1
REAL :: xphi2(mp)      ! leaf angle parmameter 2

!build hack
    REAL, DIMENSION(nrb) ::                                                     &
         cos3       ! cos(15 45 75 degrees)

    cos3 = COS(cPI180 * (/ 15.0, 45.0, 75.0 /))

    ! See Sellers 1985, eq.13 (leaf angle parameters):
    WHERE (reducedLAIdue2snow> CLAI_THRESH)
       xphi1 = 0.5 - vegxfang * (0.633 + 0.33 * vegxfang)
       xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
    END WHERE

    ! 2 dimensional LAI
    xvlai2 = SPREAD(reducedLAIdue2snow, 2, 3)

    ! Extinction coefficient for beam radiation and black leaves;
    ! eq. B6, Wang and Leuning, 1998
    WHERE (xvlai2 > CLAI_THRESH) ! vegetated
       xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
    ELSEWHERE ! i.e. bare soil
       xk = 0.0
    END WHERE


    WHERE (reducedLAIdue2snow> CLAI_THRESH ) ! vegetated

       ! Extinction coefficient for diffuse radiation for black leaves:
       ExtCoeff_dif = -LOG( SUM(                                                   &
            SPREAD( CGAUSS_W, 1, mp ) * EXP( -xk * xvlai2 ), 2) )       &
            / reducedLAIdue2snow

    ELSEWHERE ! i.e. bare soil
       ExtCoeff_dif = 0.7
    END WHERE


call calc_rhoch( c1,rhoch, mp, nrb, vegtaul, vegrefl )

rad%extkd = ExtCoeff_dif


    ! Canopy REFLection of diffuse radiation for black leaves:
    DO ictr=1,nrb

       rad%rhocdf(:,ictr) = rhoch(:,ictr) *  2. *                                &
            ( cGAUSS_W(1) * xk(:,1) / ( xk(:,1) + rad%extkd(:) )&
            + cGAUSS_W(2) * xk(:,2) / ( xk(:,2) + rad%extkd(:) )&
            + cGAUSS_W(3) * xk(:,3) / ( xk(:,3) + rad%extkd(:) ) )

    ENDDO

!H!       ! Define beam fraction, fbeam:
!H!       radfbeam(:,1) = spitter(mp,Cpi, real(metDoY), coszen, SW_down(:,1))
!H!       radfbeam(:,2) = spitter(mp,Cpi, real(metDoY), coszen, SW_down(:,2))
!H!
!H!       ! coszen is set during met data read in.
!H!
!H!       WHERE (coszen <1.0e-2)
!H!          RadFbeam(:,1) = 0.0
!H!          RadFbeam(:,2) = 0.0
!H!       END WHERE
!H!rad%fbeam = radfbeam

       ! Define beam fraction, fbeam:
       rad%fbeam(:,1) = spitter(mp,cpi, met%doy, met%coszen, met%fsd(:,1))
       rad%fbeam(:,2) = spitter(mp,cpi, met%doy, met%coszen, met%fsd(:,2))

       ! coszen is set during met data read in.

       WHERE (met%coszen <1.0e-2)
          rad%fbeam(:,1) = 0.0
          rad%fbeam(:,2) = 0.0
       END WHERE


    ! In gridcells where vegetation exists....

    !!vh !! include RAD_THRESH in condition
    WHERE (canopy%vlaiw > cLAI_THRESH .AND. met%coszen > 1.e-6 )

       ! SW beam extinction coefficient ("black" leaves, extinction neglects
       ! leaf SW transmittance and REFLectance):
       rad%extkb = xphi1 / met%coszen + xphi2

    ELSEWHERE ! i.e. bare soil
       rad%extkb = 0.5
    END WHERE

    WHERE ( ABS(rad%extkb - rad%extkd)  < 0.001 )
       rad%extkb = rad%extkd + 0.001
    END WHERE

    WHERE( met%coszen < 1.e-6 )
       ! higher value precludes sunlit leaves at night. affects
       ! nighttime evaporation - Ticket #90
       rad%extkb=1.0e5
    END WHERE

!H!    ! In gridcells where vegetation exists....
!H!
!H!    !!vh !! include RAD_THRESH in condition
!H!    WHERE (reducedLAIdue2snow> CLAI_THRESH .AND. coszen > 1.e-6 )
!H!
!H!       ! SW beam extinction coefficient ("black" leaves, extinction neglects
!H!       ! leaf SW transmittance and REFLectance):
!H!       ExtCoeff_beam= xphi1 / coszen + xphi2
!H!
!H!    ELSEWHERE ! i.e. bare soil
!H!       ExtCoeff_beam= 0.5
!H!    END WHERE
!H!
!H!    WHERE ( ABS(ExtCoeff_beam - ExtCoeff_dif)  < 0.001 )
!H!       ExtCoeff_beam= ExtCoeff_dif+ 0.001
!H!    END WHERE
!H!
!H!    WHERE( coszen < 1.e-6 )
!H!       ! higher value precludes sunlit leaves at night. affects
!H!       ! nighttime evaporation - Ticket #90
!H!       ExtCoeff_beam=1.0e5
!H!    END WHERE


  END SUBROUTINE init_radiation

END MODULE cbl_init_radiation_module
