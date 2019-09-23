MODULE cbl_albedo_mod

  IMPLICIT NONE

  PUBLIC albedo
  PRIVATE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Albedo( &
  ssnow, veg, met, rad, soil, canopy, &
#                 include "cbl_albedo_args.inc" 
                 )
!subrs called
USE cbl_rhoch_module, ONLY : calc_rhoch
USE cbl_snow_albedo_module, ONLY : surface_albedosn

    USE cable_common_module
    USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &
         canopy_type, met_type, radiation_type,      &
         !soil_snow_type, mp, r_2, nrb
         soil_snow_type, r_2

implicit none
    TYPE (canopy_type),INTENT(IN)       :: canopy
    TYPE (met_type),INTENT(INOUT)       :: met
    TYPE (radiation_type),INTENT(INOUT) :: rad
    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow

    TYPE (veg_parameter_type),INTENT(INOUT)  :: veg
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL(r_2), DIMENSION(mp)  ::                                                &
         dummy2, & !
         dummy

!    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: c1, rhoch

    LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation

    INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave

!model dimensions
!-------------------------------------------------------------------------------
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
!-------------------------------------------------------------------------------

!This is what we are returning here
!-------------------------------------------------------------------------------REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!-------------------------------------------------------------------------------

!constants
!-------------------------------------------------------------------------------
real :: Ccoszen_tols                !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)
LOGICAL :: jls_radiation            !runtime switch defined in cable_*main routines signifying this is the radiation pathway 
!-------------------------------------------------------------------------------

!masks
!-------------------------------------------------------------------------------
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL :: sunlit_mask(mp)          ! this "mp" is sunlit (uses zenith angle)
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  
!-------------------------------------------------------------------------------

integer:: surface_type(mp)          ! Integer index of Surface type (veg%iveg)

real :: reducedLAIdue2snow(mp)             ! Reduced LAI given snow coverage

! Albedos
!-------------------------------------------------------------------------------
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)
!-------------------------------------------------------------------------------

!Forcing
!-------------------------------------------------------------------------------
REAL :: MetTk(mp)                   !Air Temperture at surface - atmospheric forcing (met%tk)
REAL :: coszen(mp)                  !cosine zenith angle  (met%coszen)
REAL :: metDoY(mp)                  !Day of the Year - not always available (met%doy)
REAL :: SW_down(mp,nrb)             !Downward shortwave "forced" (met%fsd)
!-------------------------------------------------------------------------------

!Prognostics
!-------------------------------------------------------------------------------
REAL :: SnowDepth(mp)               !Total Snow depth - water eqivalent - packed from snow_surft (SnowDepth)
REAL :: SnowODepth(mp)              !Total Snow depth before any update this timestep (SnowODepth)
REAL :: SnowDensity(mp)             !Total Snow density (assumes 1 layer describes snow cover) (SnowDensity)
REAL :: SoilTemp(mp)                !Soil Temperature of top layer (soil%tgg)
REAL :: SnowAge(mp)                 !Snow age (assumes 1 layer describes snow cover) (SnowAge)
integer:: SnowFlag_3L(mp)           !Flag to treat snow as 3 layer  - if enough present. Updated depending on total depth (SnowFlag_3L)
!-------------------------------------------------------------------------------

REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings [computed in albedo() ]
!-------------------------------------------------------------------------------
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
!-------------------------------------------------------------------------------

!Variables shared primarily between radiation and albedo and possibly elsewhere
!-------------------------------------------------------------------------------
!Extinction co-efficients compued in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-efficient for Direct Beam component of SW radiation (rad%extkbm)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-efficient for Diffuse component of SW radiation (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo() 
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance  (rad%rhodf   
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance  (rad%rhobm)   
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (rad%cexpkdm)   
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (rad%cexpkbm)
!-------------------------------------------------------------------------------

!Vegetation parameters
!-------------------------------------------------------------------------------
REAL :: VegXfang(mp)                !leaf angle PARAMETER (veg%xfang)
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
!-------------------------------------------------------------------------------


call surface_albedosn( AlbSnow, AlbSoil, mp, surface_type, &
                       SnowDepth, SnowODepth, SnowFlag_3L, &
                       SnowDensity, &
                       SoilTemp, SnowAge, &
                       metTk, coszen )
   
CanopyTransmit_beam = 0.0
    rad%extkbm  = 0.0
    rad%rhocbm  = 0.0

    ! Initialise effective conopy beam reflectance:
    EffSurfRefl_beam = ssnow%albsoilsn
    EffSurfRefl_dif = ssnow%albsoilsn
    rad%albedo = ssnow%albsoilsn

    ! Define vegetation mask:
    mask = canopy%vlaiw > 0.001 .AND.                                    &
         ( met%fsd(:,1) + met%fsd(:,2) ) > 0.001

call calc_rhoch( c1,rhoch, mp, nrb, veg%taul, veg%refl )

    ! Update extinction coefficients and fractional transmittance for
    ! leaf transmittance and reflection (ie. NOT black leaves):
    !---1 = visible, 2 = nir radiaition
    DO b = 1, 2

       rad%extkdm(:,b) = rad%extkd * c1(:,b)

       !--Define canopy diffuse transmittance (fraction):
CanopyTransmit_dif(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)

       !---Calculate effective diffuse reflectance (fraction):
       WHERE( canopy%vlaiw > 1e-2 )                                             &
            EffSurfRefl_dif(:,b) = rad%rhocdf(:,b) + (ssnow%albsoilsn(:,b)             &
            - rad%rhocdf(:,b)) * canopytransmit_dif(:,b)**2

       !---where vegetated and sunlit
       WHERE (mask)

          rad%extkbm(:,b) = rad%extkb * c1(:,b)

          ! Canopy reflection (6.21) beam:
          rad%rhocbm(:,b) = 2. * rad%extkb / ( rad%extkb + rad%extkd )          &
               * rhoch(:,b)

          ! Canopy beam transmittance (fraction):
          dummy2 = MIN(rad%extkbm(:,b)*canopy%vlaiw, 20.)
          dummy  = EXP(-dummy2)
CanopyTransmit_beam(:,b) = REAL(dummy)

          ! Calculate effective beam reflectance (fraction):
          EffSurfRefl_beam(:,b) = rad%rhocbm(:,b) + (ssnow%albsoilsn(:,b)             &
               - rad%rhocbm(:,b))*canopytransmit_beam(:,b)**2

       END WHERE

       ! Define albedo:
       WHERE( canopy%vlaiw> .001 )                                      &
            rad%albedo(:,b) = ( 1. - rad%fbeam(:,b) )*EffSurfRefl_dif(:,b) +           &
            rad%fbeam(:,b) * EffSurfRefl_beam(:,b)

    END DO




   
!H!!Modify albedo based on snow coverage 
!H!call surface_albedosn( AlbSnow, AlbSoil, mp, surface_type, &
!H!                       SnowDepth, SnowODepth, SnowFlag_3L, &
!H!                       SnowDensity, &
!H!                       SoilTemp, SnowAge, &
!H!                       metTk, coszen )
!H!   
!H!! Define canopy Reflectance for diffuse/direct radiation
!H!! Formerly rad%rhocbm, rad%rhocdf
!H!call CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif, &
!H!                        mp, nrb, CGauss_w, sunlit_veg_mask, &
!H!                        AlbSnow, xk, rhoch,                  &
!H!                        ExtCoeff_beam, ExtCoeff_dif)
!H!
!H!! Define canopy diffuse transmittance 
!H!! Formerly rad%cexpkbm, rad%cexpkdm
!H!call CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,&
!H!                              sunlit_veg_mask, reducedLAIdue2snow, &
!H!                              EffExtCoeff_dif, EffExtCoeff_beam)
!H!  
!H!! Finally compute Effective 4-band albedo for diffuse/direct radiation- 
!H!! In the UM this is the required variable to be passed back on the rad call
!H!! Formerly rad%reffbm, rad%reffdf
!H!call EffectiveSurfaceReflectance( EffSurfRefl_beam, EffSurfRefl_dif,           &
!H!                                  mp, nrb, veg_mask, sunlit_veg_mask,          &
!H!                                  CanopyRefl_beam, CanopyRefl_dif,             &
!H!                                  CanopyTransmit_beam,CanopyTransmit_dif,      &
!H!                                  AlbSnow )
!H!
!H!! Compute total albedo to SW given the Effective Surface Reflectance 
!H!! (considering Canopy/Soil/Snow contributions) 
!H!! we dont need to do this on rad call AND may not haveappropriate RadFbeam
!H!if(.NOT. jls_radiation) &
!H!  call FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam, &
!H!                       EffSurfRefl_dif, EffSurfRefl_beam, AlbSnow )

END SUBROUTINE albedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyReflectance( CanopyRefl_beam, CanopyRefl_dif, &
                         mp, nrb, CGauss_w, sunlit_veg_mask, &
                         AlbSnow, xk, rhoch,                  &
                         ExtCoeff_beam, ExtCoeff_dif)
implicit none 
!re-decl in args
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
real :: Cgauss_w(nrb)
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: xk(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)

! Initialise canopy beam reflectance:
CanopyRefl_beam  = AlbSnow !Formerly rad%reffbm
CanopyRefl_dif   = AlbSnow ! Formerly rad%refdfm

call CanopyReflectance_beam( CanopyRefl_beam, mp, nrb, sunlit_veg_mask, &
                             ExtCoeff_beam, ExtCoeff_dif, rhoch )

call CanopyReflectance_dif( CanopyRefl_dif, mp, nrb, CGauss_w, &
                             ExtCoeff_dif, xk, rhoch )
End subroutine CanopyReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyReflectance_beam( CanopyRefl_beam, mp, nrb, sunlit_veg_mask, &
              ExtCoeff_beam,ExtCoeff_dif, rhoch )
implicit none
integer :: mp
integer ::nrb 
real :: CanopyRefl_beam(mp,nrb)
real :: ExtCoeff_dif(mp) 
real :: ExtCoeff_beam(mp) 
LOGICAL :: sunlit_veg_mask(mp) 
REAL :: rhoch(mp, nrb)
integer :: i, b

! Canopy reflection (6.21) beam:
DO i = 1,mp
  DO b = 1, 2
    IF( sunlit_veg_mask(i) ) &
      CanopyRefl_beam(i,b) = 2. * ExtCoeff_beam(i) / &
                            ( ExtCoeff_beam(i) + ExtCoeff_dif(i) )          & 
                            * rhoch(i,b)
  END DO
END DO

End subroutine CanopyReflectance_beam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyReflectance_dif( CanopyRefl_dif, mp, nrb, CGauss_w,  &
                                  ExtCoeff_dif, xk, rhoch )

implicit none
INTEGER :: mp
integer :: nrb
real :: Cgauss_w(nrb)
REAL :: CanopyRefl_dif(mp,nrb)  
real :: ExtCoeff_dif(mp)    
REAL :: xk(mp,nrb)      ! extinct. coef.for beam rad. and black leaves
REAL :: rhoch(mp,nrb)      
!local vars
INTEGER :: ictr

! Canopy REFLection of diffuse radiation for black leaves:
DO ictr=1,nrb

  CanopyRefl_dif(:,ictr) = rhoch(:,ictr) *  2. *                                &
                       ( CGAUSS_W(1) * xk(:,1) / ( xk(:,1) + ExtCoeff_dif(:) )&
                       + CGAUSS_W(2) * xk(:,2) / ( xk(:,2) + ExtCoeff_dif(:) )&
                       + CGAUSS_W(3) * xk(:,3) / ( xk(:,3) + ExtCoeff_dif(:) ) )

ENDDO

End subroutine CanopyReflectance_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyTransmitance(CanopyTransmit_beam, CanopyTransmit_dif, mp, nrb,&
                              sunlit_veg_mask, reducedLAIdue2snow, &
                              EffExtCoeff_dif, EffExtCoeff_beam)
implicit none
!re-decl in args
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
REAL :: CanopyTransmit_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyTransmit_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  
real :: reducedLAIdue2snow(mp)
REAL :: EffExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: EffExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)

!Zero initialization remains the calculated value where "mask"=FALSE
CanopyTransmit_beam = 0.0 ! Formerly rad%rhocbm

call CanopyTransmitance_dif(CanopyTransmit_dif, mp, nrb, EffExtCoeff_dif, reducedLAIdue2snow)

call CanopyTransmitance_beam( CanopyTransmit_beam, mp, nrb, EffExtCoeff_beam,  &
                              reducedLAIdue2snow, sunlit_veg_mask )

End subroutine CanopyTransmitance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyTransmitance_dif(CanopyTransmit, mp, nrb, Eff_Transmit, reducedLAIdue2snow )
implicit none
integer :: mp 
integer :: nrb
real :: CanopyTransmit(mp,nrb) 
real :: Eff_Transmit(mp,nrb) 
real :: reducedLAIdue2snow(mp)
integer :: i, b
 
DO i = 1,mp
  DO b = 1, 2
    CanopyTransmit(i,b) = EXP( -1.* Eff_Transmit(i,b) * reducedLAIdue2snow(i) )
  enddo
enddo

End subroutine  CanopyTransmitance_dif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CanopyTransmitance_beam(CanopyTransmit, mp, nrb, Eff_Transmit, reducedLAIdue2snow, mask )
implicit none
integer :: mp 
integer :: nrb
real :: CanopyTransmit(mp,nrb) 
real :: Eff_Transmit(mp,nrb) 
real :: reducedLAIdue2snow(mp)
logical :: mask(mp) 
real :: dummy(mp,nrb) 
integer :: i, b
 
DO i = 1,mp
  DO b = 1, 2 !ithis is fixed as 2  because nrb=3 due to legacy  
    dummy(i,b) = min( Eff_Transmit(i,b) * reducedLAIdue2snow(i), 20. )
    CanopyTransmit(i,b) = EXP( -1.* dummy(i,b) )
  enddo
enddo

End subroutine  CanopyTransmitance_beam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EffectiveSurfaceReflectance(EffSurfRefl_beam, EffSurfRefl_dif,      &
                                       mp, nrb, veg_mask, sunlit_veg_mask,     &
                                       CanopyRefl_beam, CanopyRefl_dif,        &
                                       CanopyTransmit_beam,CanopyTransmit_dif, & 
                                       AlbSnow )
implicit none
!re-decl input args 
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is vegetated (uses minimum LAI) 
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
REAL :: CanopyRefl_beam(mp,nrb)  
REAL :: CanopyRefl_dif(mp,nrb)  
REAL :: CanopyTransmit_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyTransmit_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
real :: AlbSnow(mp,nrb)



call EffectiveReflectance( EffSurfRefl_dif, mp, nrb, CanopyRefl_dif, AlbSnow, &
                           CanopyTransmit_dif, veg_mask )

call EffectiveReflectance( EffSurfRefl_beam, mp, nrb, CanopyRefl_beam, AlbSnow,&
                           CanopyTransmit_beam, sunlit_veg_mask )

End subroutine EffectiveSurfaceReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EffectiveReflectance( EffRefl, mp, nrb, CanopyRefl, AlbSnow, &
          CanopyTransmit, mask )
implicit none
integer :: mp
integer :: nrb
real :: AlbSnow(mp,nrb)
real :: CanopyRefl(mp,nrb)
real :: CanopyTransmit(mp,nrb) 
real :: EffRefl(mp,nrb) 
logical :: mask(mp) 
integer :: i,b  

DO i = 1,mp
  DO b = 1, 2!ithis is fixed as 2  because nrb=3 due to legacy  
      IF( mask(i) ) then 
      
         ! Calculate effective beam reflectance (fraction):
         EffRefl(i,b) = CanopyRefl(i,b) &
                            + ( AlbSnow(i,b) - CanopyRefl(i,b) ) &
                            * CanopyTransmit(i,b)**2

    endif
  END DO
END DO

End subroutine EffectiveReflectance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FbeamRadAlbedo( RadAlbedo, mp, nrb, veg_mask, radfbeam, &
                           EffSurfRefl_dif, EffSurfRefl_beam, AlbSnow )
implicit none
!re-decl input args 
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
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
    IF( veg_mask(i) )                                      &
       RadAlbedo(i,b) = ( 1. - radfbeam(i,b) )*EffSurfRefl_dif(i,b) +           &
                         radfbeam(i,b) * EffSurfRefl_beam(i,b)
  END DO
END DO

End subroutine FbeamRadAlbedo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cbl_albedo_mod
