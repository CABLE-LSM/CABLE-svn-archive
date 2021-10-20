MODULE cbl_albedo_mod

   USE cable_data_module, ONLY : ialbedo_type, point2constants 
   
   IMPLICIT NONE
   
  PUBLIC albedo
   PRIVATE

   TYPE(ialbedo_type) :: C


CONTAINS

SUBROUTINE Albedo( AlbSnow, AlbSoil,              & 
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
   
!SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy)
   USE cable_common_module   
   USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &     
                                   canopy_type, met_type, radiation_type,      &
                                   soil_snow_type, r_2
   
USE cable_um_tech_mod, ONLY : ssnow, veg, met, rad, soil, canopy
!subrs called
USE cbl_snow_albedo_module, ONLY : surface_albedosn
   
implicit none

!model dimensions
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

!This is what we are returning here
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)

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

real :: SumEffSurfRefl_beam(1)
real :: SumEffSurfRefl_dif(1)
integer :: i

    ! END header

   REAL(r_2), DIMENSION(mp)  ::                                                &
      dummy2, & !
      dummy

   LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation

   INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
      
   CALL point2constants(C) 
   
!Modify parametrised soil albedo based on snow coverage 
call surface_albedosn( AlbSnow, AlbSoil, mp, nrb, jls_radiation, surface_type, soil_type, &
                       SnowDepth, SnowODepth, SnowFlag_3L,                      & 
                       SnowDensity, SoilTemp, SnowTemp, SnowAge,                     & 
                       MetTk, Coszen )


   rad%cexpkbm = 0.0
   rad%extkbm  = 0.0
   rad%rhocbm  = 0.0

   ! Initialise effective conopy beam reflectance:
   rad%reffbm = ssnow%albsoilsn
   rad%reffdf = ssnow%albsoilsn
   rad%albedo = ssnow%albsoilsn

   ! Define vegetation mask:
   mask = canopy%vlaiw > C%LAI_THRESH .AND.                                    &
          ( met%fsd(:,1) + met%fsd(:,2) ) > C%RAD_THRESH     

   CALL calc_rhoch( veg, c1, rhoch )

   ! Update extinction coefficients and fractional transmittance for 
   ! leaf transmittance and reflection (ie. NOT black leaves):
   !---1 = visible, 2 = nir radiaition
   DO b = 1, 2        
      
      rad%extkdm(:,b) = rad%extkd * c1(:,b)
   
      !--Define canopy diffuse transmittance (fraction):
      rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)

      !---Calculate effective diffuse reflectance (fraction):
      WHERE( canopy%vlaiw > C%lai_thresh )                                             &
         rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssnow%albsoilsn(:,b)             &
                           - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2
      
      !---where vegetated and sunlit 
      WHERE (mask)                
      
         rad%extkbm(:,b) = rad%extkb * c1(:,b)
      
      ! Canopy reflection (6.21) beam:
         rad%rhocbm(:,b) = 2. * rad%extkb / ( rad%extkb + rad%extkd )          &
                        * rhoch(:,b)

         ! Canopy beam transmittance (fraction):
         dummy2 = -rad%extkbm(:,b)*canopy%vlaiw
         dummy  = EXP(dummy2)

         rad%cexpkbm(:,b) = REAL(dummy)

         ! Calculate effective beam reflectance (fraction):
         rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssnow%albsoilsn(:,b)             &
               - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2

      END WHERE

      ! Define albedo:
      WHERE( canopy%vlaiw> C%LAI_THRESH )                                      &
         rad%albedo(:,b) = ( 1. - rad%fbeam(:,b) )*rad%reffdf(:,b) +           &
                           rad%fbeam(:,b) * rad%reffbm(:,b)
       
   END DO

END SUBROUTINE Albedo 

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

!jhan:subr was reintroduced here to temporarily resolve issue when 
!creating libcable.a  (repeated in cable_radiation.F90)
SUBROUTINE calc_rhoch(veg,c1,rhoch) 

   USE cable_def_types_mod, ONLY : veg_parameter_type
   TYPE (veg_parameter_type), INTENT(INOUT) :: veg
   REAL, INTENT(INOUT), DIMENSION(:,:) :: c1, rhoch

   c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
   c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
   c1(:,3) = 1.
    
   ! Canopy reflection black horiz leaves 
   ! (eq. 6.19 in Goudriaan and van Laar, 1994):
   rhoch = (1.0 - c1) / (1.0 + c1)

END SUBROUTINE calc_rhoch 


END MODULE cbl_albedo_mod
