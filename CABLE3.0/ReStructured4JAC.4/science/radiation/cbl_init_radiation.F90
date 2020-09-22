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

  IMPLICIT NONE

  PUBLIC init_radiation
  PRIVATE

!FUDGED local pars -masks tuned at other times - review conssitency!!
real :: Ccoszen_tols_huge  ! 1e-4 * threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Ccoszen_tols_tiny  ! 1e-4 * threshold cosine of sun's zenith angle, below which considered SUNLIT

CONTAINS
SUBROUTINE init_radiation( ExtCoeff_beam, ExtCoeff_dif,                        &
                        EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,           &
                        c1, rhoch, xk,                                         &
                        mp,nrb,                                                &
                        Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,      &
                        cbl_standalone, jls_standalone, jls_radiation,         &
                        subr_name,                                             &
                        veg_mask, sunlit_mask, sunlit_veg_mask,                &
                        VegXfang, VegTaul, VegRefl,                            &
                        coszen, metDoY, SW_down,                               & 
                        reducedLAIdue2snow )
 !  SUBROUTINE init_radiation( met, rad, veg, canopy, veg_mask, sunlit_veg_mask)

    ! Alternate version of init_radiation that uses only the
    ! zenith angle instead of the fluxes. This means it can be called
    ! before the cable albedo calculation.
    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
         veg_parameter_type!, nrb, mp
    USE cable_common_module

USE cbl_spitter_module, ONLY : Spitter
USE cbl_rhoch_module, ONLY : calc_rhoch
!USE cable_math_constants_mod, ONLY : CPI => PI
!USE cable_math_constants_mod, ONLY : CPI180 => PI180
!USE cable_other_constants_mod, ONLY : CLAI_thresh => LAI_thresh
!USE cable_other_constants_mod, ONLY : CGauss_W => Gauss_W


!re-decl input args
!model dimensions
integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

!returned variables
REAL :: ExtCoeff_beam(mp)       !"raw" Extinction co-efficient for Direct Beam component of SW radiation
REAL :: ExtCoeff_dif(mp)        !"raw"Extinction co-efficient for Diffuse component of SW radiation
REAL :: EffExtCoeff_beam(mp,nrb)!Effective Extinction co-efficient for Direct Beam component of SW radiation
REAL :: EffExtCoeff_dif(mp,nrb) !Effective Extinction co-efficient for Diffuse component of SW radiation
REAL :: RadFbeam(mp,nrb)        !Beam Fraction of Downward SW radiation [formerly rad%fbeam]


!constants
real :: Clai_thresh             !threshold LAI below which considered UN-vegetated
real :: Ccoszen_tols            !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)
real :: Cpi                     !PI - from cable_math_constants originally
real :: Cpi180                  !PI in radians - from cable_math_constants originally
!what model am i in
LOGICAL :: cbl_standalone       !runtime switch defined in cable_*main routines signifying this is cable_standalone
LOGICAL :: jls_standalone       !runtime switch defined in cable_*main routines signifying this is jules_standalone
LOGICAL :: jls_radiation        !runtime switch defined in cable_*main routines signifying this is the radiation pathway 
character(len=*) :: subr_name !where am i called from
!masks
logical :: veg_mask(mp)         !vegetated mask [formed by comparrisson of LAI CLAI_thresh ]
logical :: sunlit_mask(mp)      !sunlit mask [formed by comparrisson of coszen to coszen_tols i.e. is the sun up]
logical :: sunlit_veg_mask(mp)  !combined mask - BOTH sunlit and vegetated

!vegetation parameters input via namelist
REAL :: VegXfang(mp)
REAL :: VegTaul(mp,nrb)
REAL :: VegRefl(mp,nrb)

REAL :: reducedLAIdue2snow(mp)  !Effective LAI given (potential sno coverage)
REAL :: coszen(mp)              ! cosine zenith angle of sun
REAL :: SW_down(mp,nrb)         !Downward SW radiation [formerly met%fsd]
integer :: metDoY(mp)           !Day of the Year [formerly met%doy]

!co-efficients used throughout init_radiation used in albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)              ! extinct. coef.for beam rad. and black leaves

!local_vars - common scaling co-efficients used throughout init_radiation
REAL :: xvlai2(mp,nrb) ! 2D vlai
REAL :: xphi1(mp)      ! leaf angle parmameter 1
REAL :: xphi2(mp)      ! leaf angle parmameter 2

    REAL, DIMENSION(nrb) ::                                                     &
         cos3       ! cos(15 45 75 degrees)

    INTEGER :: ictr

    cos3 = COS(CPI180 * (/ 15.0, 45.0, 75.0 /))

!Null Initializations
ExtCoeff_beam(:) = 0.0
ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0
EffExtCoeff_dif(:,:) = 0.0
RadFbeam(:,:) = 0.0
c1(:,:) = 0.0
rhoch(:,:) = 0.0
xk(:,:) = 0.0

! Compute common scaling co-efficients used throughout init_radiation
call Common_InitRad_Scalings( xphi1, xphi2, xk, xvlai2, c1, rhoch,             &
                            mp, nrb, Cpi180,cLAI_thresh, veg_mask,             &
                            reducedLAIdue2snow, VegXfang, VegTaul, VegRefl)

!Limiting Initializations for stability
Ccoszen_tols_huge = Ccoszen_tols * 1e2 
Ccoszen_tols_tiny = Ccoszen_tols * 1e-2 

    WHERE ( veg_mask )

       ! Extinction coefficient for diffuse radiation for black leaves:
       ExtCoeff_dif = -LOG( SUM(                                                   &
            SPREAD( CGAUSS_W, 1, mp ) * EXP( -xk * xvlai2 ), 2) )       &
            / reducedLAIdue2snow

    ELSEWHERE ! i.e. bare soil
       ExtCoeff_dif = 0.7
    END WHERE


CALL calc_rhoch( c1,rhoch, mp, nrb, VegTaul, VegRefl )

    IF( .NOT. cable_runtime%um) THEN

       ! Define beam fraction, fbeam:
       RadFbeam(:,1) = spitter(mp, cpi, MetDoy, coszen, SW_down(:,1))
       RadFbeam(:,2) = spitter(mp, cpi, MetDoy, coszen, SW_down(:,2))

       ! coszen is set during met data read in.

       WHERE (coszen <1.0e-2)
          RadFbeam(:,1) = 0.0
          RadFbeam(:,2) = 0.0
       END WHERE

    ENDIF

    ! In gridcells where vegetation exists....
    WHERE ( sunlit_veg_mask )

       ! SW beam extinction coefficient ("black" leaves, extinction neglects
       ! leaf SW transmittance and REFLectance):
       ExtCoeff_beam = xphi1 / coszen + xphi2

    ELSEWHERE ! i.e. bare soil
       ExtCoeff_beam = 0.5
    END WHERE

    WHERE ( ABS(ExtCoeff_beam - ExtCoeff_dif)  < 0.001 )
       ExtCoeff_beam = ExtCoeff_dif + 0.001
    END WHERE

    WHERE( coszen < 1.e-6 )
       ! higher value precludes sunlit leaves at night. affects
       ! nighttime evaporation - Ticket #90
       ExtCoeff_beam=1.0e5
    END WHERE

  END SUBROUTINE init_radiation
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

  xphi1 = 0.0
  xphi2 = 0.0
  xvlai2 = 0.0
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


END MODULE cbl_init_radiation_module
