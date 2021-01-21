MODULE cable_albedo_module

  IMPLICIT NONE

  PUBLIC albedo
  PRIVATE

CONTAINS

SUBROUTINE Albedo( &
AlbSnow, AlbSoil,                                 & 
mp, nrb,                                          &
jls_radiation ,                                   &
veg_mask, sunlit_mask, sunlit_veg_mask,           &  
Ccoszen_tols, CGAUSS_W,                           & 
surface_type, soil_type, VegRefl, VegTaul,        &
metTk, coszen,                                    & 
reducedLAIdue2snow,                               &
SnowDepth, SnowODepth, SnowFlag_3L,               & 
SnowDensity, SoilTemp, SnowTemp, SnowAge,                   &
xk, c1, rhoch,                                    & 
RadFbeam, RadAlbedo,                              &
ExtCoeff_dif, ExtCoeff_beam,                      &
EffExtCoeff_dif, EffExtCoeff_beam,                &
CanopyRefl_dif,CanopyRefl_beam,                   &
CanopyTransmit_dif, CanopyTransmit_beam,          &
EffSurfRefl_dif, EffSurfRefl_beam                 )

!subrs called
USE cbl_rhoch_module, ONLY : calc_rhoch
    USE cable_common_module
    USE cable_def_types_mod, ONLY : r_2
USE cbl_snow_albedo_module, ONLY : surface_albedosn
USE cbl_rhoch_module, ONLY : calc_rhoch

implicit none

!model dimensions
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

!This is what we are returning here
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (EffSurfRefl_dif)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (EffSurfRefl_beam)

!constants
real :: Ccoszen_tols                !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb)
LOGICAL :: jls_radiation            !runtime switch def. in cable_*main routines 
                                    !signifying this is the radiation pathway 

!masks
LOGICAL :: veg_mask(mp)             ! this "mp" is vegetated (uses minimum LAI) 
LOGICAL :: sunlit_mask(mp)          ! this "mp" is sunlit (uses zenith angle)
LOGICAL :: sunlit_veg_mask(mp)      ! this "mp" is BOTH sunlit AND  vegetated  

!Vegetation parameters
REAL :: VegTaul(mp,nrb)             !PARAMETER leaf transmisivity (veg%taul)
REAL :: VegRefl(mp,nrb)             !PARAMETER leaf reflectivity (veg%refl)
integer:: surface_type(mp)          !Integer index of Surface type (veg%iveg)
integer:: soil_type(mp)          !Integer index of Soil    type (soil%isoilm)

real :: reducedLAIdue2snow(mp)      !Reduced LAI given snow coverage

! Albedos
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)

!Forcing
REAL :: MetTk(mp)                   !Air Temperture at surface - atmospheric forcing (met%tk)
REAL :: coszen(mp)                  !cosine zenith angle  (met%coszen)
REAL :: metDoY(mp)                  !Day of the Year - not always available (met%doy)
REAL :: SW_down(mp,nrb)             !Downward shortwave "forced" (met%fsd)

!Prognostics
REAL :: SnowDepth(mp)               !Total Snow depth - water eqivalent - packed from snow_surft (SnowDepth)
REAL :: SnowODepth(mp)              !Total Snow depth before any update this timestep (SnowODepth)
REAL :: SnowDensity(mp)             !Total Snow density (assumes 1 layer describes snow cover) (SnowDensity)
REAL :: SoilTemp(mp)                !Soil Temperature of top layer (soil%tgg)
REAL :: SnowTemp(mp)                !Soil Temperature of top layer (soil%tgg)
REAL :: SnowAge(mp)                 !Snow age (assumes 1 layer describes snow cover) (SnowAge)
integer:: SnowFlag_3L(mp)           !Flag to treat snow as 3 layer - if enough present 
                                    !Updated depending on total depth (SnowFlag_3L)

REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings - computed  in init_radiation()
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)

!Variables shared primarily between radiation and albedo and possibly elsewhere
!Extinction co-efficients computed in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient 
                                    !Direct Beam component of SW radiation (ExtCoeff_beam)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for 
                                    !Diffuse component of SW radiation (ExtCoeff_dif)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-eff 
                                    !Direct Beam component of SW radiation (ExtCoeff_beam)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-eff 
                                    !Diffuse component of SW radiation (ExtCoeff_difm)

!Canopy reflectance/transmitance compued in albedo() 
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance  (CanopyRefl_dif   
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance  (CanopyRefl_beam)   
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (CanopyTransmit_dif)   
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (CanopyTransmit_beam)

real :: SumEffSurfRefl_beam(1)
real :: SumEffSurfRefl_dif(1)
integer :: i

    REAL(r_2), DIMENSION(mp)  ::                                                &
         dummy2, & !
         dummy

    INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave

    ! END header

!CanopyTransmit_dif(:,:) = 0.0 !open in 7589
!CanopyTransmit_beam(:,:) = 0.0 !open in 7589
!!CanopyRefl_dif(:,:) = 0.0
!CanopyRefl_beam(:,:) = 0.0 !open in 7589
!AlbSnow(:,:) = 0.0 !open in 7589

!Modify parametrised soil albedo based on snow coverage 
call surface_albedosn( AlbSnow, AlbSoil, mp, nrb, jls_radiation, surface_type, soil_type, &
                       SnowDepth, SnowODepth, SnowFlag_3L,                      & 
                       SnowDensity, SoilTemp, SnowTemp, SnowAge,                     & 
                       MetTk, Coszen )

! Initialise effective conopy beam reflectance:
EffSurfRefl_beam = AlbSnow
EffSurfRefl_dif = AlbSnow
RadAlbedo = AlbSnow

CALL calc_rhoch( c1,rhoch, mp, nrb, VegTaul, VegRefl )

    ! Update extinction coefficients and fractional transmittance for
    ! leaf transmittance and reflection (ie. NOT black leaves):
    !---1 = visible, 2 = nir radiaition
    DO b = 1, 2
       EffExtCoeff_dif(:,b) = ExtCoeff_dif(:) * c1(:,b)
    END DO

    DO b = 1, 2
         !---where vegetated and sunlit
       WHERE (sunlit_veg_mask)
          EffExtCoeff_beam(:,b) = ExtCoeff_beam(:) * c1(:,b)
       END WHERE
    END DO


    DO b = 1, 2
       !--Define canopy diffuse transmittance (fraction):
       CanopyTransmit_dif(:,b) = EXP(-EffExtCoeff_dif(:,b) * reducedLAIdue2snow)
    END DO
    DO b = 1, 2
       !---Calculate effective diffuse reflectance (fraction):
       WHERE( veg_mask )                                      &
            EffSurfRefl_dif(:,b) = CanopyRefl_dif(:,b) + (AlbSnow(:,b)             &
            - CanopyRefl_dif(:,b)) * CanopyTransmit_dif(:,b)**2
    END DO

    ! Canopy beam transmittance (fraction):
    DO b = 1, 2
       WHERE (sunlit_veg_mask)

          ! Canopy reflection (6.21) beam:
          CanopyRefl_beam(:,b) = 2. * ExtCoeff_beam / ( ExtCoeff_beam + ExtCoeff_dif )          &
               * rhoch(:,b)

          ! Canopy beam transmittance (fraction):
          dummy2 = MIN(EffExtCoeff_beam(:,b)*reducedLAIdue2snow, 20.)
          dummy  = EXP(-dummy2)
          CanopyTransmit_beam(:,b) = REAL(dummy)
       END WHERE
    END DO

    DO b = 1, 2
       !---where vegetated and sunlit
       WHERE (sunlit_veg_mask)
          ! Calculate effective beam reflectance (fraction):
          EffSurfRefl_beam(:,b) = CanopyRefl_beam(:,b) + (AlbSnow(:,b)             &
               - CanopyRefl_beam(:,b))*CanopyTransmit_beam(:,b)**2
       END WHERE
    END DO
 
    DO b = 1, 2
       ! Define albedo:
       WHERE( veg_mask )                                      &
            RadAlbedo(:,b) = ( 1. - RadFbeam(:,b) )*EffSurfRefl_dif(:,b) +           &
            RadFbeam(:,b) * EffSurfRefl_beam(:,b)
    END DO
  
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
!HACHvstrunk!CanopyRefl_beam  = AlbSnow !Formerly rad%reffbm
!HACHvstrunk!CanopyRefl_dif   = AlbSnow ! Formerly rad%refdfm

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


End MODULE cable_albedo_module

