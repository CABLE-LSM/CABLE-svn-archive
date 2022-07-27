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
! Purpose: Computes 4-band (visible/near-infrared, beam/diffuse) reflectances
!          and albedo
!
! IMPORTANT NOTE regarding the masks (Ticket 333)
! Prior to #333, 3 masks were used here - veg_mask, sunlit_mask, veg_mask
! - although passed in sunlit_mask was not used
! For JAC we will not be able to populate the sunlit masks and instead will 
! evaluate the EffExtCoeff outside their bounds of applicability here by 
! using inclusive masks in their place from the calling routines,
!  ie. JAC will use veg_mask in place of veg_mask
!
! To avoid confusion the mask names here are renamed:
! sunlit_mask now called mask1, veg_mask now called mask2
! 
! ==============================================================================

MODULE cbl_albedo_mod

IMPLICIT NONE

PUBLIC albedo
PRIVATE

CONTAINS

SUBROUTINE Albedo( AlbSnow, AlbSoil, mp, nrb, ICE_SoilType, lakes_cable,       & 
                   jls_radiation, veg_mask, Ccoszen_tols, cgauss_w,            &
                   SurfaceType, SoilType, VegRefl, VegTaul,                    &
                   coszen, reducedLAIdue2snow, SnowDepth, SnowDensity,         &
                   SoilTemp, SnowAge, xk, c1, rhoch, RadFbeam, RadAlbedo,      &
                   ExtCoeff_dif, ExtCoeff_beam, EffExtCoeff_dif,               &
                   EffExtCoeff_beam, CanopyRefl_dif,CanopyRefl_beam,           &
                   CanopyTransmit_dif, CanopyTransmit_beam,                    &
                   EffSurfRefl_dif, EffSurfRefl_beam )

!subrs called
USE cbl_snow_albedo_module, ONLY: surface_albedosn

IMPLICIT NONE

!model dimensions
INTEGER :: mp                       ! total number of "tiles"
INTEGER :: nrb                      ! # rad bands: VIS,NIR. 3rd dim was for LW

! Return: Effective Surface Relectance as seen by atmosphere
REAL :: EffSurfRefl_dif(mp,nrb)     ! Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    ! Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)

!--- IN: CABLE specific surface_type indexes 
INTEGER, INTENT(IN) :: ICE_SoilType
INTEGER, INTENT(IN) :: lakes_cable

!constants
REAL :: Ccoszen_tols                !threshold cosine of sun's zenith angle, below which considered SUNLIT
REAL :: Cgauss_w(nrb)
LOGICAL :: jls_radiation            !runtime switch def. in cable_*main routines
                                    !signifying this is the radiation pathway

!masks
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI)

!Vegetation parameters
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
INTEGER:: SurfaceType(mp)           !Integer index of Surface type (veg%iveg)
INTEGER:: SoilType(mp)              !Integer index of Soil    type (soil%isoilm)

REAL :: reducedLAIdue2snow(mp)      !Reduced LAI given snow coverage

! Albedos
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)

!Forcing
REAL :: coszen(mp)                  !cosine zenith angle  (met%coszen)
REAL :: SW_down(mp,nrb)             !Downward shortwave "forced" (met%fsd)

!Prognostics
REAL :: SnowDepth(mp)               !Total Snow depth - water eqivalent - packed from snow_surft (SnowDepth)
REAL :: SnowDensity(mp)             !Total Snow density (assumes 1 layer describes snow cover) (SnowDensity)
REAL :: SoilTemp(mp)                !Soil Temperature of top layer - for lake alebdo (ssnow%tgg)
REAL :: SnowAge(mp)                 !Snow age (assumes 1 layer describes snow cover) (SnowAge)

REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings - computed  in init_radiation()
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)

!Variables shared primarily between radiation and albedo and possibly elsewhere
!Extinction co-efficients computed in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient
                                    !Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for
                                    !Diffuse component of SW radiation (rad%extkd)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-eff
                                    !Direct Beam component of SW radiation (rad%extkbm)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-eff
                                    !Diffuse component of SW radiation (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo()
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance  (rad%rhocdf
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance  (rad%rhocbm)
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (rad%cexpkdm)
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (rad%cexpkbm)

REAL :: SumEffSurfRefl_beam(1)
REAL :: SumEffSurfRefl_dif(1)
INTEGER :: i

    ! END header

AlbSnow(:,:) = 0.0
CanopyTransmit_beam(:,:) = 1.0
CanopyRefl_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0
CanopyTransmit_dif(:,:) = 1.0  ! MPI (at least inits this = 1.0 at dt=0)

!Modify parametrised soil albedo based on snow coverage
CALL surface_albedosn( AlbSnow, AlbSoil, mp, nrb, ICE_SoilType, lakes_cable,      &
                       SurfaceType, SoilType, SnowDepth, SnowDensity,         &
                       SoilTemp, SnowAge, Coszen )

! Update fractional leaf transmittance and reflection
!---1 = visible, 2 = nir radiaition

! Define canopy Reflectance for diffuse/direct radiation
! Formerly rad%rhocbm, rad%rhocdf
CALL CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif,                       &
                        mp, nrb, CGauss_w, veg_mask,                           &
                        AlbSnow, xk, rhoch,                                    &
                        ExtCoeff_beam, ExtCoeff_dif)

! Define canopy diffuse transmittance
! Formerly rad%cexpkbm, rad%cexpkdm
CALL CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,      &
                        veg_mask, reducedLAIdue2snow,                          &
                        EffExtCoeff_dif, EffExtCoeff_beam)

!---1 = visible, 2 = nir radiaition
! Finally compute Effective 4-band albedo for diffuse/direct radiation-
! In the UM this is the required variable to be passed back on the rad call
! Formerly rad%reffbm, rad%reffdf

! Even when there is no vegetation, albedo is at least snow modified soil albedo
EffSurfRefl_dif = AlbSnow
EffSurfRefl_beam = AlbSnow

CALL EffectiveSurfaceReflectance( EffSurfRefl_beam, EffSurfRefl_dif,           &
                                  mp, nrb, veg_mask, CanopyRefl_beam,          &
                                  CanopyRefl_dif, CanopyTransmit_beam,         &
                                  CanopyTransmit_dif, AlbSnow )

! Compute total albedo to SW given the Effective Surface Reflectance
! (considering Canopy/Soil/Snow contributions)
! we dont need to do this on rad call AND may not haveappropriate RadFbeam
RadAlbedo = AlbSnow
IF (.NOT. jls_radiation)                                                       &
  CALL FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam,                 &
                       EffSurfRefl_dif, EffSurfRefl_beam, AlbSnow )

END SUBROUTINE albedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif,                 &
                         mp, nrb, CGauss_w, veg_mask,                   &
                         AlbSnow, xk, rhoch,                                   &
                         ExtCoeff_beam, ExtCoeff_dif)
IMPLICIT NONE
!re-decl in args
INTEGER :: mp                       !total number of "tiles"
INTEGER :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm)
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)
REAL :: Cgauss_w(nrb)
LOGICAL :: veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: xk(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)

! Initialise canopy beam reflectance:
!HACHvstrunk!CanopyRefl_beam  = AlbSnow !Formerly rad%reffbm
!HACHvstrunk!CanopyRefl_dif   = AlbSnow ! Formerly rad%refdfm

CALL CanopyReflectance_beam( CanopyRefl_beam, mp, nrb, veg_mask,        &
                             ExtCoeff_beam, ExtCoeff_dif, rhoch )

CALL CanopyReflectance_dif( CanopyRefl_dif, mp, nrb, CGauss_w,                 &
                             ExtCoeff_dif, xk, rhoch )
END SUBROUTINE CanopyReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyReflectance_beam( CanopyRefl_beam, mp, nrb, veg_mask,  &
              ExtCoeff_beam,ExtCoeff_dif, rhoch )
IMPLICIT NONE
INTEGER :: mp
INTEGER ::nrb
REAL :: CanopyRefl_beam(mp,nrb)
REAL :: ExtCoeff_dif(mp)
REAL :: ExtCoeff_beam(mp)
LOGICAL :: veg_mask(mp)
REAL :: rhoch(mp, nrb)
INTEGER :: i, b

! Canopy reflection (6.21) beam:
DO i = 1,mp
  DO b = 1, 2
    IF ( veg_mask(i) )                                                  &
      CanopyRefl_beam(i,b) = 2.0 * ExtCoeff_beam(i) /                          &
                            ( ExtCoeff_beam(i) + ExtCoeff_dif(i) )             &
                            * rhoch(i,b)
  END DO
END DO

END SUBROUTINE CanopyReflectance_beam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyReflectance_dif( CanopyRefl_dif, mp, nrb, CGauss_w,           &
                                  ExtCoeff_dif, xk, rhoch )

IMPLICIT NONE
INTEGER :: mp
INTEGER :: nrb
REAL :: Cgauss_w(nrb)
REAL :: CanopyRefl_dif(mp,nrb)
REAL :: ExtCoeff_dif(mp)
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves
REAL :: rhoch(mp,nrb)
!local vars
INTEGER :: ictr

! Canopy REFLection of diffuse radiation for black leaves:
DO ictr=1,2

  CanopyRefl_dif(:,ictr) = rhoch(:,ictr) *  2.0 *                              &
                       ( cgauss_w(1) * xk(:,1) / ( xk(:,1) + ExtCoeff_dif(:) ) &
                       + cgauss_w(2) * xk(:,2) / ( xk(:,2) + ExtCoeff_dif(:) ) &
                       + cgauss_w(3) * xk(:,3) / ( xk(:,3) + ExtCoeff_dif(:) ) )

END DO

END SUBROUTINE CanopyReflectance_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,&
                              mask, reducedLAIdue2snow,                        &
                              EffExtCoeff_dif, EffExtCoeff_beam)
IMPLICIT NONE
!re-decl in args
INTEGER :: mp                       !total number of "tiles"
INTEGER :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: CanopyTransmit_dif(mp,nrb)      !Canopy Transmitance (rad%cexpkdm)
REAL :: CanopyTransmit_beam(mp,nrb)     !Canopy Transmitance (rad%cexpkbm)
LOGICAL :: mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated
LOGICAL :: dummyMask(mp)
REAL :: reducedLAIdue2snow(mp)
REAL :: EffExtCoeff_beam(mp,nrb)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: EffExtCoeff_dif(mp,nrb)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)

! For beam, compute canopy trasmitance when sunlit (and vegetated)
CALL CanopyTransmitance_beam( CanopyTransmit_beam, mp, nrb, EffExtCoeff_beam,  &
                              reducedLAIdue2snow, mask )

!'=1.0' initialization remains the calculated value where "mask"=FALSE
dummyMask(:) = .TRUE.

! For diffuse rad, always compute canopy trasmitance
CALL CanopyTransmitance_dif( CanopyTransmit_dif, mp, nrb, EffExtCoeff_dif,     &
                             reducedLAIdue2snow, dummyMask )

END SUBROUTINE CanopyTransmitance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CanopyTransmitance_dif(CanopyTransmit, mp, nrb, ExtinctionCoeff, reducedLAIdue2snow, mask )
IMPLICIT NONE
INTEGER :: mp
INTEGER :: nrb
LOGICAL :: mask(mp)
REAL :: CanopyTransmit(mp,nrb)
REAL :: ExtinctionCoeff(mp,nrb)
REAL :: reducedLAIdue2snow(mp)
REAL :: dummy(mp,nrb)
INTEGER :: i, b

DO i = 1,mp
  DO b = 1, 2
    dummy(i,b) = ExtinctionCoeff(i,b) * reducedLAIdue2snow(i)
    CanopyTransmit(i,b) = EXP( -1.0* dummy(i,b) )
  END DO
END DO

END SUBROUTINE  CanopyTransmitance_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CanopyTransmitance_beam(CanopyTransmit, mp, nrb, ExtinctionCoeff, reducedLAIdue2snow, mask )
IMPLICIT NONE
INTEGER :: mp
INTEGER :: nrb
LOGICAL :: mask(mp)
REAL :: CanopyTransmit(mp,nrb)
REAL :: ExtinctionCoeff(mp,nrb)
REAL :: reducedLAIdue2snow(mp)
REAL :: dummy(mp,nrb)
INTEGER :: i, b

DO i = 1,mp
  DO b = 1, 2
    IF ( mask(i) ) THEN
      dummy(i,b) = MIN( ExtinctionCoeff(i,b) * reducedLAIdue2snow(i), 20.0 )
      CanopyTransmit(i,b) = EXP( -1.0* dummy(i,b) )
    END IF
  END DO
END DO

END SUBROUTINE  CanopyTransmitance_beam


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EffectiveSurfaceReflectance(EffSurfRefl_beam, EffSurfRefl_dif,      &
                                       mp, nrb, veg_mask, CanopyRefl_beam,     &
                                       CanopyRefl_dif, CanopyTransmit_beam,    &
                                       CanopyTransmit_dif, AlbSnow )
IMPLICIT NONE

INTEGER :: mp                       !total number of "tiles"
INTEGER :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
REAL :: CanopyRefl_beam(mp,nrb)
REAL :: CanopyRefl_dif(mp,nrb)
REAL :: CanopyTransmit_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm)
REAL :: CanopyTransmit_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)
REAL :: AlbSnow(mp,nrb)

CALL EffectiveReflectance( EffSurfRefl_dif, mp, nrb, CanopyRefl_dif, AlbSnow,  &
                           CanopyTransmit_dif, veg_mask )

CALL EffectiveReflectance( EffSurfRefl_beam, mp, nrb, CanopyRefl_beam, AlbSnow,&
                           CanopyTransmit_beam, veg_mask )

END SUBROUTINE EffectiveSurfaceReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE EffectiveReflectance( EffRefl, mp, nrb, CanopyRefl, AlbSnow,        &
                                 CanopyTransmit, mask )
IMPLICIT NONE
INTEGER :: mp
INTEGER :: nrb
REAL :: AlbSnow(mp,nrb)
REAL :: CanopyRefl(mp,nrb)
REAL :: CanopyTransmit(mp,nrb)
REAL :: EffRefl(mp,nrb)
LOGICAL :: mask(mp)
INTEGER :: i,b

DO i = 1,mp
  DO b = 1, 2!ithis is fixed as 2  because nrb=3 due to legacy
    IF ( mask(i) ) THEN

       ! Calculate effective beam reflectance (fraction):
      EffRefl(i,b) = CanopyRefl(i,b)                                           &
                         + ( AlbSnow(i,b) - CanopyRefl(i,b) )                  &
                         * CanopyTransmit(i,b)**2

    END IF
  END DO
END DO

END SUBROUTINE EffectiveReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam,             &
                           EffSurfRefl_dif, EffSurfRefl_beam, AlbSnow )
IMPLICIT NONE
!re-decl input args
INTEGER :: mp                       !total number of "tiles"
INTEGER :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)
REAL :: AlbSnow(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI)
REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!local vars
INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
INTEGER :: i

! Initialise total albedo:
RadAlbedo = AlbSnow
DO i = 1,mp
  DO b = 1, 2 !nrb -1 -nrb shouldnt be =3 anyway
    ! Define albedo:
    IF ( veg_mask(i) )                                                         &
       RadAlbedo(i,b) = ( 1.0 - radfbeam(i,b) )*EffSurfRefl_dif(i,b) +         &
                         radfbeam(i,b) * EffSurfRefl_beam(i,b)
  END DO
END DO

END SUBROUTINE FbeamRadAlbedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cbl_albedo_mod
