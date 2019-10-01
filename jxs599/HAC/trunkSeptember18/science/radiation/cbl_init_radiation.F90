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

!local vars
REAL :: cos3(nrb)      ! cos(15 45 75 degrees)

cos3 = COS(CPI180 * (/ 15.0, 45.0, 75.0 /))

call Common_InitRad_Scalings( xphi1, xphi2, xk, xvlai2, c1, rhoch,      &
                            mp, nrb, Cpi180,cLAI_thresh, veg_mask,             &
                            reducedLAIdue2snow,                &
                            VegXfang, VegTaul, VegRefl)

    WHERE (reducedLAIdue2snow> CLAI_THRESH ) ! vegetated

       ! Extinction coefficient for diffuse radiation for black leaves:
       ExtCoeff_dif = -LOG( SUM(                                                   &
            SPREAD( CGAUSS_W, 1, mp ) * EXP( -xk * xvlai2 ), 2) )       &
            / reducedLAIdue2snow

    ELSEWHERE ! i.e. bare soil
       ExtCoeff_dif = 0.7
    END WHERE

rad%extkd = ExtCoeff_dif

    ! Canopy REFLection of diffuse radiation for black leaves:
    DO ictr=1,nrb

       rad%rhocdf(:,ictr) = rhoch(:,ictr) *  2. *                                &
            ( cGAUSS_W(1) * xk(:,1) / ( xk(:,1) + rad%extkd(:) )&
            + cGAUSS_W(2) * xk(:,2) / ( xk(:,2) + rad%extkd(:) )&
            + cGAUSS_W(3) * xk(:,3) / ( xk(:,3) + rad%extkd(:) ) )

    ENDDO

! Define beam fraction, fbeam:
radfbeam(:,1) = spitter(mp,Cpi, real(metDoY), coszen, SW_down(:,1))
radfbeam(:,2) = spitter(mp,Cpi, real(metDoY), coszen, SW_down(:,2))

! coszen is set during met data read in.
WHERE (coszen <1.0e-2)
   RadFbeam(:,1) = 0.0
   RadFbeam(:,2) = 0.0
END WHERE

rad%fbeam = radfbeam

    ! In gridcells where vegetation exists....

    WHERE (reducedLAIdue2snow> CLAI_THRESH .AND. coszen > 1.e-6 )

       ! SW beam extinction coefficient ("black" leaves, extinction neglects
       ! leaf SW transmittance and REFLectance):
       ExtCoeff_beam= xphi1 / coszen + xphi2

    ELSEWHERE ! i.e. bare soil
       ExtCoeff_beam= 0.5
    END WHERE

    WHERE ( ABS(ExtCoeff_beam - ExtCoeff_dif)  < 0.001 )
       ExtCoeff_beam= ExtCoeff_dif+ 0.001
    END WHERE

    WHERE( coszen < 1.e-6 )
       ! higher value precludes sunlit leaves at night. affects
       ! nighttime evaporation - Ticket #90
       ExtCoeff_beam=1.0e5
    END WHERE
rad%extkb = ExtCoeff_beam

END SUBROUTINE init_radiation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Common_InitRad_Scalings( xphi1, xphi2, xk, xvlai2, c1, rhoch,      &
                            mp, nrb, Cpi180,cLAI_thresh, veg_mask,             &
                            reducedLAIdue2snow,                &
                            VegXfang, VegTaul, VegRefl)
!subrs
USE cbl_rhoch_module,   ONLY : calc_rhoch
implicit none
!re-decl in args
integer :: mp
integer :: nrb
real :: Cpi180
real :: cLAI_thresh
real :: xphi1(mp)    ! leaf angle parmameter 1
real :: xphi2(mp)    ! leaf angle parmameter 2
REAL :: xvlai2(mp,nrb)  ! 2D vlai
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
real :: reducedLAIdue2snow(mp)
real :: VegXfang(mp)
REAL :: VegTaul(mp,nrb)
REAL :: VegRefl(mp,nrb)
logical :: veg_mask(mp)

call common_InitRad_coeffs( xphi1, xphi2, xk, xvlai2, mp, nrb, Cpi180,&
                            cLAI_thresh, veg_mask, VegXfang, reducedLAIdue2snow  )

CALL calc_rhoch( c1,rhoch, mp, nrb, VegTaul, VegRefl )

End subroutine Common_InitRad_Scalings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  common_InitRad_coeffs( xphi1, xphi2, xk, xvlai2, mp, nrb, Cpi180,&
                            cLAI_thresh, veg_mask, VegXfang, reducedLAIdue2snow  )


implicit none
!re-decl in args
integer :: mp
integer :: nrb
real :: Cpi180
real :: xphi1(mp)    ! leaf angle parmameter 1
real :: xphi2(mp)    ! leaf angle parmameter 2
REAL :: xvlai2(mp,nrb)  ! 2D vlai
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves
real :: VegXfang(mp)
real :: reducedLAIdue2snow(mp)
logical :: veg_mask(mp)
real:: cLAI_thresh

!local vars
REAL :: cos3(nrb)      ! cos(15 45 75 degrees)

cos3 = COS(CPI180 * (/ 15.0, 45.0, 75.0 /))

! See Sellers 1985, eq.13 (leaf angle parameters):
WHERE ( veg_mask )
  xphi1 = 0.5 - VegXfang * (0.633 + 0.33 * VegXfang)
  xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
END WHERE

! 2 dimensional LAI
xvlai2 = SPREAD(reducedLAIdue2snow, 2, 3)

! Extinction coefficient for beam radiation and black leaves;
! eq. B6, Wang and Leuning, 1998
WHERE (xvlai2 > cLAI_THRESH) ! vegetated
   xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
ELSEWHERE ! i.e. bare soil
   xk = 0.0
END WHERE

End subroutine common_InitRad_coeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtinctionCoeff( ExtCoeff_beam, ExtCoeff_dif, mp, nrb, CGauss_w, reducedLAIdue2snow, &
                            sunlit_mask, veg_mask, sunlit_veg_mask,  &
                            cLAI_thresh, coszen, xphi1, xphi2, xk, xvlai2)

implicit none
!re-decl in args
integer :: mp
integer :: nrb
real :: ExtCoeff_beam(mp)        !extinction co-eff RETURNED
real :: ExtCoeff_dif(mp)         !extinction co-eff RETURNED
logical:: veg_mask(mp)           !vegetated mask based on a minimum LAI 
logical :: sunlit_mask(mp)       !sunlit mask based on zenith angle
logical :: sunlit_veg_mask(mp)   !BOTH sunlit and vegetated mask 
real :: Cgauss_w(nrb)
real :: cLAI_thresh
real :: coszen(mp)
real :: reducedLAIdue2snow(mp)
real :: xphi1(mp)
real :: xphi2(mp)
REAL :: xvlai2(mp,nrb)  ! 2D vlai
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves

call ExtinctionCoeff_dif( ExtCoeff_dif, mp, nrb, CGauss_w, reducedLAIdue2snow, &
                          veg_mask, cLAI_thresh, xk, xvlai2)

call ExtinctionCoeff_beam( ExtCoeff_beam, mp, nrb, &
                           sunlit_mask, veg_mask, sunlit_veg_mask,  &
                           coszen, xphi1, xphi2, ExtCoeff_dif )

End subroutine ExtinctionCoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtinctionCoeff_beam( ExtCoeff_beam, mp, nrb, sunlit_mask, veg_mask, sunlit_veg_mask,  &
                                 coszen, xphi1, xphi2, ExtCoeff_dif )

implicit none
integer :: mp
integer :: nrb

logical :: sunlit_mask(mp)       !sunlit mask based on zenith angle
logical :: sunlit_veg_mask(mp)   !BOTH sunlit and vegetated mask 
real :: coszen(mp)
real :: ExtCoeff_beam(mp)
real :: xphi1(mp)
real :: xphi2(mp)
real :: ExtCoeff_dif(mp)
logical:: veg_mask(mp)

! retain this initialization for bare soil
ExtCoeff_beam = 0.5

! SW beam extinction coefficient ("black" leaves, extinction neglects
! leaf SW transmittance and REFLectance):
WHERE ( sunlit_veg_mask ) &
  ExtCoeff_beam = xphi1 / Coszen + xphi2

! higher value precludes sunlit leaves at night. affects
! nighttime evaporation - Ticket #90 
WHERE( .NOT. sunlit_mask ) ExtCoeff_beam = 1.0e5 

! Seems to be for stability only
WHERE ( abs(ExtCoeff_beam - ExtCoeff_dif )  < 0.001 ) &
  ExtCoeff_beam = ExtCoeff_dif + 0.001

End subroutine ExtinctionCoeff_beam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtinctionCoeff_dif( ExtCoeff_dif, mp, nrb, Cgauss_w, reducedLAIdue2snow, veg_mask, &
                                 cLAI_thresh, xk, xvlai2)
implicit none
integer :: mp
integer :: nrb
real :: Cgauss_w(nrb)

real :: ExtCoeff_dif(mp)    !return Extinctino Coefficient 
logical :: veg_mask(mp)
real :: reducedLAIdue2snow(mp)
real:: cLAI_thresh
REAL :: xvlai2(mp,nrb)  ! 2D vlai
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves


!local vars
REAL :: cos3(nrb)      ! cos(15 45 75 degrees)

ExtCoeff_dif = 0.7

WHERE ( veg_mask ) ! vegetated
  ! Extinction coefficient for diffuse radiation for black leaves:
  ExtCoeff_dif = -LOG( SUM(                                                   &
                            SPREAD( CGAUSS_W, 1, mp ) &
                            * EXP( -xk * xvlai2 ), 2) &
                     )                                                         &
                     / reducedLAIdue2snow

   END WHERE
End subroutine ExtinctionCoeff_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine EffectiveExtinctCoeffs( EffExtCoeff_beam, EffExtCoeff_dif, mp, nrb, &
                                   sunlit_veg_mask,                        &
                                   ExtCoeff_beam, ExtCoeff_dif, c1 )
implicit none
integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

REAL :: EffExtCoeff_beam(mp,nrb)!Effective Extinction co-efficient for Direct Beam component of SW radiation
REAL :: EffExtCoeff_dif(mp,nrb) !Effective Extinction co-efficient for Diffuse component of SW radiation

REAL :: c1(mp,nrb)
logical :: sunlit_veg_mask(mp)  !combined mask - BOTH sunlit and vegetated
REAL :: ExtCoeff_beam(mp)       !"raw" Extinction co-efficient for Direct Beam component of SW radiation
REAL :: ExtCoeff_dif(mp)        !"raw"Extinction co-efficient for Diffuse component of SW radiation

EffExtCoeff_dif = 0.
call EffectiveExtinctCoeff( EffExtCoeff_dif, mp, ExtCoeff_dif, c1 )

EffExtCoeff_beam = 0.
call EffectiveExtinctCoeff( EffExtCoeff_beam, mp, ExtCoeff_beam, c1, &
                            sunlit_veg_mask )
End Subroutine EffectiveExtinctCoeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! modified k diffuse(6.20)(for leaf scattering)
subroutine EffectiveExtinctCoeff(Eff_ExtCoeff, mp, ExtCoeff, c1, mask )
implicit none
integer :: mp 
real :: Eff_ExtCoeff(mp,2) 
real :: ExtCoeff(mp) 
real :: c1(mp,2) 
logical, optional :: mask(mp) 
integer :: i, b

DO i = 1,mp
  DO b = 1, 2
    !IF mask is present we are doing the beam component then: 
    if( present(mask)) then 
      !then ONLY IF it is sunlit and vegetated -else default 
      if( mask(i) ) Eff_ExtCoeff(i,b) = ExtCoeff(i) * c1(i,b)
    else         
      Eff_ExtCoeff(i,b) = ExtCoeff(i) * c1(i,b)          
    endif

  enddo
enddo

End subroutine EffectiveExtinctCoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BeamFraction( RadFbeam, mp, nrb, Cpi,sunlit_mask, metDoy, coszen, SW_down ) 
USE cbl_spitter_module, ONLY : Spitter

integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: RadFbeam(mp,nrb)        !Beam Fraction of Downward SW radiation [formerly rad%fbeam]

real :: Cpi                     !PI - from cable_math_constants originally
logical :: sunlit_mask(mp)       !sunlit mask based on zenith angle

real :: metDoY(mp)          !Day of the Year [formerly met%doy]
real :: coszen(mp)          !Day of the Year [formerly met%doy]
REAL :: SW_down(mp,nrb)     !Downward SW radiation [formerly met%fsd]


! Define beam fraction, fbeam:
RadFbeam(:,1) = spitter(mp, cpi, real(metDoy), coszen, SW_down(:,1))
RadfBeam(:,2) = spitter(mp, cpi, real(metDoy), coszen, SW_down(:,2))

! coszen is set during met data read in.
WHERE ( .NOT. sunlit_mask )
  RadFbeam(:,1) = 0.0
  RadFbeam(:,2) = 0.0
END WHERE

End subroutine BeamFraction



END MODULE cbl_init_radiation_module
