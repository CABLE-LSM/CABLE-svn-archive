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

USE cbl_rhoch_module,   ONLY : calc_rhoch
USE cbl_spitter_module, ONLY : Spitter
USE cable_um_tech_mod, ONLY : canopy, met, rad, veg 
   USE cable_common_module

implicit none

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
   
!logical ::mask(mp)         !vegetated mask [formed by comparrisson of LAI CLAI_thresh ]
!local vars
integer :: ictr
REAL :: cos3(nrb)      ! cos(15 45 75 degrees)

  
!Null Initializations
ExtCoeff_beam(:) = 0.0
ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0
EffExtCoeff_dif(:,:) = 0.0
RadFbeam(:,:) = 0.0
c1(:,:) = 0.0
rhoch(:,:) = 0.0
xk(:,:) = 0.0

!Limiting Initializations for stability
Ccoszen_tols_huge = Ccoszen_tols * 1e2 
Ccoszen_tols_tiny = Ccoszen_tols * 1e-2 

   CALL point2constants( C ) 
   
   cos3 = COS(C%PI180 * (/ 15.0, 45.0, 75.0 /))

   ! See Sellers 1985, eq.13 (leaf angle parameters):
   WHERE (canopy%vlaiw > C%LAI_THRESH)
      xphi1 = 0.5 - veg%xfang * (0.633 + 0.33 * veg%xfang)
      xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
   END WHERE

   ! 2 dimensional LAI
   xvlai2 = SPREAD(canopy%vlaiw, 2, 3)

   ! Extinction coefficient for beam radiation and black leaves;
   ! eq. B6, Wang and Leuning, 1998
   WHERE (xvlai2 > C%LAI_THRESH) ! vegetated
      xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
   ELSEWHERE ! i.e. bare soil
      xk = 0.0          
   END WHERE

   WHERE (canopy%vlaiw > C%LAI_THRESH ) ! vegetated
   
      ! Extinction coefficient for diffuse radiation for black leaves:
      rad%extkd = -LOG( SUM(                                                   &
                  SPREAD( C%GAUSS_W, 1, mp ) * EXP( -xk * xvlai2 ), 2) )       &
                  / canopy%vlaiw
   
   ELSEWHERE ! i.e. bare soil
      rad%extkd = 0.7
   END WHERE

CALL calc_rhoch( c1,rhoch, mp, nrb, veg%taul, veg%refl )

   ! Canopy REFLection of diffuse radiation for black leaves:
   DO ictr=1,nrb
     
     rad%rhocdf(:,ictr) = rhoch(:,ictr) *                                      &
                          ( C%GAUSS_W(1) * xk(:,1) / ( xk(:,1) + rad%extkd(:) )&
                          + C%GAUSS_W(2) * xk(:,2) / ( xk(:,2) + rad%extkd(:) )&
                          + C%GAUSS_W(3) * xk(:,3) / ( xk(:,3) + rad%extkd(:) ) )

   ENDDO
   
   WHERE (canopy%vlaiw > C%LAI_THRESH)    
      
      ! SW beam extinction coefficient ("black" leaves, extinction neglects
      ! leaf SW transmittance and REFLectance):
      rad%extkb = xphi1 / met%coszen + xphi2
   
   ELSEWHERE ! i.e. bare soil
      rad%extkb = 0.5
   END WHERE
   
   WHERE ( abs(rad%extkb - rad%extkd)  < 0.001 )
      rad%extkb = rad%extkd + 0.001
   END WHERE
   
   WHERE(rad%fbeam(:,1) < C%RAD_THRESH )
      rad%extkb=30.0         ! keep cexpkbm within real*4 range (BP jul2010)
   END WHERE
!borrowed from cbl_iniit_radiation IN CABLE3 - bc there it is moved out oof Albedo and into iniit_radiation
!so here we are missing this definition   
rad%extkbm(:,:) = 0.0
rad%extkdm(:,:) = 0.0

call EffectiveExtinctCoeffs( rad%extkbm, rad%extkdm,               &
                             mp, nrb, sunlit_veg_mask,                        &
                             rad%extkb, rad%extkd, c1 )

! Offline/standalone forcing gives us total downward Shortwave. We have
! previosuly, arbitratily split this into NIR/VIS (50/50). We use 
! Spitter function to split these bands into direct beam and diffuse components
IF( cbl_standalone .OR. jls_standalone .AND. .NOT. jls_radiation ) &
  CALL BeamFraction( RadFbeam, mp, nrb, Cpi, Ccoszen_tols_huge, metDoy,  &
                     coszen, SW_down ) 
   
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
real :: Eff_ExtCoeff(mp,3) 
real :: ExtCoeff(mp) 
real :: c1(mp,3) 
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
subroutine BeamFraction( RadFbeam, mp, nrb, Cpi,Ccoszen_tols_huge, metDoy, &
coszen, SW_down ) 
USE cbl_spitter_module, ONLY : Spitter

integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: RadFbeam(mp,nrb)        !Beam Fraction of Downward SW radiation [formerly rad%fbeam]

real :: Cpi !PI - from cable_math_constants originally
real :: Ccoszen_tols_huge !PI - from cable_math_constants originally

integer:: metDoY(mp)          !Day of the Year [formerly met%doy]
real :: coszen(mp)          !Day of the Year [formerly met%doy]
REAL :: SW_down(mp,nrb)     !Downward SW radiation [formerly met%fsd]


! Define beam fraction, fbeam:
RadFbeam(:,1) = spitter(mp, cpi, metDoy, coszen, SW_down(:,1))
RadfBeam(:,2) = spitter(mp, cpi, metDoy, coszen, SW_down(:,2))

! coszen is set during met data read in.
WHERE (coszen < 1.e-2 )
  RadFbeam(:,1) = 0.0
  RadFbeam(:,2) = 0.0
END WHERE

End subroutine BeamFraction



END MODULE cbl_init_radiation_module
