MODULE cable_albedo_module

   USE cable_data_module, ONLY : ialbedo_type, point2constants 
   
   IMPLICIT NONE
   
   PUBLIC surface_albedo
   PRIVATE

   TYPE(ialbedo_type) :: C


CONTAINS

  
!d1!SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy, &
!d1!AlbSnow, AlbSoil,              & 
!d1!mp, nrb,                                          &
!d1!jls_radiation ,                                   &
!d1!veg_mask, sunlit_mask, sunlit_veg_mask,           &  
!d1!Ccoszen_tols, CGAUSS_W,                           & 
!d1!surface_type, soil_type, VegRefl, VegTaul,        &
!d1!metTk, coszen,                                    & 
!d1!reducedLAIdue2snow,                               &
!d1!SnowDepth, SnowODepth, SnowFlag_3L,               & 
!d1!SnowDensity, SoilTemp, SnowTemp, SnowAge,                   &
!d1!xk, c1, rhoch,                                    & 
!d1!RadFbeam, RadAlbedo,                              &
!d1!ExtCoeff_dif, ExtCoeff_beam,                      &
!d1!EffExtCoeff_dif, EffExtCoeff_beam,                &
!d1!CanopyRefl_dif,CanopyRefl_beam,                   &
!d1!CanopyTransmit_dif, CanopyTransmit_beam,          &
!d1!EffSurfRefl_dif, EffSurfRefl_beam                 )
   
SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy)
   USE cable_common_module   
   USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &     
                                   canopy_type, met_type, radiation_type,      &
                                   soil_snow_type, r_2, mp, nrb
   
!subrs called
!CBL3!USE cbl_snow_albedo_module, ONLY : surface_albedosn
   
implicit none

!model dimensions
!d1!integer :: mp                       !total number of "tiles"  
!d1!integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

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

!SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy)
   TYPE (canopy_type),INTENT(IN)       :: canopy
   TYPE (met_type),INTENT(INOUT)       :: met
   TYPE (radiation_type),INTENT(INOUT) :: rad
   TYPE (soil_snow_type),INTENT(INOUT) :: ssnow

   TYPE (veg_parameter_type),INTENT(INOUT)  :: veg
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil   

   REAL(r_2), DIMENSION(mp)  ::                                                &
      dummy2, & !
      dummy

   LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation

   INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
      
   CALL point2constants(C) 
   
   CALL surface_albedosn(ssnow, veg, met, soil)

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
      WHERE( canopy%vlaiw > 1e-2 )                                             &
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

END SUBROUTINE surface_albedo 

! ------------------------------------------------------------------------------

SUBROUTINE surface_albedosn(ssnow, veg, met, soil)
   
   USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &     
                                   met_type, soil_snow_type, mp 
   USE cable_common_module
   
   TYPE (soil_snow_type),INTENT(INOUT) :: ssnow
   TYPE (met_type),INTENT(INOUT)       :: met
   
   TYPE (veg_parameter_type),INTENT(INout)  :: veg
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil   

   REAL, DIMENSION(mp) ::                                                      &
      alv,     &  ! Snow albedo for visible
      alir,    &  ! Snow albedo for near infra-red
      ar1,     &  ! crystal growth  (-ve)
      ar2,     &  ! freezing of melt water
      ar3,     &  !
      dnsnow,  &  ! new snow albedo
      dtau,    &  !
      fage,    &  ! age factor
      fzenm,   &  !    
      sfact,   &  !
      snr,     &  ! 
      snrat,   &  !  
      talb,    &  ! snow albedo
      tmp         ! temporary value
   
   REAL, PARAMETER ::                                                          &
      alvo  = 0.95,  &  ! albedo for vis. on a new snow
      aliro = 0.70      ! albedo for near-infr. on a new snow
   
   INTEGER :: k,i,j,l,l1,l2

   soil%albsoilf = soil%albsoil(:,1)

   ! lakes: hard-wired number to be removed in future
   WHERE( veg%iveg == 16 )                                                     &
      soil%albsoilf = -0.022*( MIN( 275., MAX( 260., met%tk ) ) - 260. ) + 0.45

   WHERE(ssnow%snowd > 1. .and. veg%iveg == 16 ) soil%albsoilf = 0.85

   sfact = 0.68
  
   WHERE (soil%albsoilf <= 0.14)
      sfact = 0.5
   ELSEWHERE (soil%albsoilf > 0.14 .and. soil%albsoilf <= 0.20)
      sfact = 0.62
   END WHERE

   ssnow%albsoilsn(:,2) = 2. * soil%albsoilf / (1. + sfact)
   ssnow%albsoilsn(:,1) = sfact * ssnow%albsoilsn(:,2)
  
   snrat=0.
   alir =0.
   alv  =0.

   WHERE ( ssnow%snowd > 1. .AND. .NOT. cable_runtime%um_radiation ) 
       
      ! new snow (cm H2O)
      dnsnow = MIN ( 1., .1 * MAX( 0., ssnow%snowd - ssnow%osnowd ) ) 
      
      ! Snow age depends on snow crystal growth, freezing of melt water,
      ! accumulation of dirt and amount of new snow.
      tmp = ssnow%isflag * ssnow%tggsn(:,1) + ( 1 - ssnow%isflag )            &
            * ssnow%tgg(:,1)
      tmp = MIN( tmp, C%TFRZ )
      ar1 = 5000. * (1. / (C%TFRZ-0.01) - 1. / tmp) ! crystal growth  (-ve)
      ar2 = 10. * ar1 ! freezing of melt water
      snr = ssnow%snowd / max (ssnow%ssdnn, 200.)
      
      WHERE (soil%isoilm == 9)
         ! permanent ice: hard-wired number to be removed in future version
         ar3 = .0000001
         !  NB. dsnow =1,assumes pristine snow; ignores soot etc. ALTERNATIVELY,
         !dnsnow = max (dnsnow, .5) !increase refreshing of snow in Antarctic
         dnsnow = 1.0
         snrat = 1.
      
      ELSEWHERE

         ! accumulation of dirt
         ar3 = .1
         ! snow covered fraction of the grid
         snrat = min (1., snr / (snr + .1) )    
      
      END WHERE

      dtau = 1.e-6 * (EXP( ar1 ) + EXP( ar2 ) + ar3 ) * kwidth_gl 
      
      WHERE (ssnow%snowd <= 1.0)
         ssnow%snage = 0.
      ELSEWHERE
         ssnow%snage = max (0.,(ssnow%snage+dtau)*(1.-dnsnow))
      END WHERE
      
      fage = 1. - 1. / (1. + ssnow%snage ) !age factor

      tmp = MAX( .17365, met%coszen )
      fzenm = MAX( 0.0, MERGE( 0.0,                                           &
              ( 1. + 1./2. ) / ( 1. + 2.*2. * tmp ) - 1./2., tmp > 0.5 ) )

      tmp = alvo * (1.0 - 0.2 * fage)
      alv = .4 * fzenm * (1. - tmp) + tmp
      tmp = aliro * (1. - .5 * fage)

      ! use dry snow albedo for pernament land ice: hard-wired no to be removed
      WHERE (soil%isoilm == 9)             
         
         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
         tmp = 0.75 * (1. - .5 * fage)
      
      END WHERE
      
      alir = .4 * fzenm * (1.0 - tmp) + tmp
      talb = .5 * (alv + alir) ! snow albedo
    
   ENDWHERE        ! snowd > 0
   
   ! when it is called from cable_rad_driver (UM) 
   ! no need to recalculate snage 
   WHERE (ssnow%snowd > 1 .and. cable_runtime%um_radiation )
      
      snr = ssnow%snowd / MAX (ssnow%ssdnn, 200.)
      
      WHERE (soil%isoilm == 9)
         ! permanent ice: hard-wired number to be removed
         snrat = 1.
      ELSEWHERE
         snrat = MIN (1., snr / (snr + .1) )
      END WHERE
      
      fage = 1. - 1. / (1. + ssnow%snage ) !age factor
      tmp = MAX (.17365, met%coszen )
      fzenm = MAX( 0., MERGE( 0.0,                                             &
              ( 1. + 1./2. ) / ( 1. + 2. * 2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)
      
      ! use dry snow albedo
      WHERE (soil%isoilm == 9)          
         ! permanent ice: hard-wired number to be removed

         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
         tmp = 0.75 * (1. - .5 * fage)
      
      END WHERE
      
      alir = .4 * fzenm * (1.0 - tmp) + tmp
      talb = .5 * (alv + alir) ! snow albedo
    
   ENDWHERE        ! snowd > 0
   

   ssnow%albsoilsn(:,2) = MIN( aliro,                                          &
                          ( 1. - snrat ) * ssnow%albsoilsn(:,2) + snrat * alir)
   
   ssnow%albsoilsn(:,1) = MIN( alvo,                                           &
                          ( 1. - snrat ) * ssnow%albsoilsn(:,1) + snrat * alv )

   WHERE (soil%isoilm == 9)          ! use dry snow albedo
     ssnow%albsoilsn(:,2) = 0.82
     ssnow%albsoilsn(:,1) = 0.82
   END WHERE
   
END SUBROUTINE surface_albedosn

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


END MODULE cable_albedo_module
