
MODULE albedo_module

   IMPLICIT NONE
   PRIVATE
   PUBLIC surface_albedo


   CONTAINS


   SUBROUTINE surface_albedo(ssoil, veg, met, rad, soil, canopy)
   
   USE cable_common_module
   USE define_types
   USE define_dimensions
   USE other_constants, ONLY : LAI_THRESH, RAD_THRESH 
   USE cable_diag_module, ONLY : cable_stat
   
   TYPE (canopy_type),INTENT(IN)          :: canopy
   TYPE (met_type),INTENT(INOUT)       :: met
   TYPE (radiation_type),INTENT(INOUT) :: rad
   TYPE (soil_snow_type),INTENT(INOUT) :: ssoil

   TYPE (veg_parameter_type),INTENT(INOUT)  :: veg
      TYPE(soil_parameter_type), INTENT(INOUT) :: soil   

   REAL(r_2), DIMENSION(mp)  ::                                                &
      dummy2, & !
      dummy

   REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: c1, rhoch
   
      LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation

   INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave
      
   ! END header

   IF (.NOT. allocated(c1)) &
      ALLOCATE( c1(mp,nrb), rhoch(mp,nrb) )

!jhan:this changes Bondville offline
   !CALL surface_albedosn(ssoil, veg, met, soil)

   WHERE (soil%isoilm == 9)          ! use dry snow albedo
     ssoil%albsoilsn(:,2) = 0.82
     ssoil%albsoilsn(:,1) = 0.82
   END WHERE
   
   rad%cexpkbm = 0.0
   rad%extkbm  = 0.0
   rad%rhocbm  = 0.0

   ! Initialise effective conopy beam reflectance:
   rad%reffbm = ssoil%albsoilsn
   rad%reffdf = ssoil%albsoilsn
   rad%albedo = ssoil%albsoilsn

   ! Define vegetation mask:
   mask = canopy%vlaiw > LAI_THRESH .AND.                                      &
          ( met%fsd(:,1) + met%fsd(:,2) ) > RAD_THRESH     

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
         rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssoil%albsoilsn(:,b)             &
                           - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2
      
      !---where vegetated and sunlit 
      WHERE (mask)                
      
         rad%extkbm(:,b) = rad%extkb * c1(:,b)
      
      ! Canopy reflection (6.21) beam:
      rad%rhocbm(:,b) = 2. * rad%extkb / ( rad%extkb + rad%extkd )             &
                        * rhoch(:,b)

         ! Canopy beam transmittance (fraction):
         dummy2 = -rad%extkbm(:,b)*canopy%vlaiw
         dummy  = EXP(dummy2)

         rad%cexpkbm(:,b) = REAL(dummy, r_1)

         ! Calculate effective beam reflectance (fraction):
         rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssoil%albsoilsn(:,b)             &
               - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2

      END WHERE

      ! Define albedo:
      WHERE( canopy%vlaiw> LAI_THRESH )                                        &
         rad%albedo(:,b) = ( 1. - rad%fbeam(:,b) )*rad%reffdf(:,b) +           &
                           rad%fbeam(:,b) * rad%reffbm(:,b)
       
   END DO

END SUBROUTINE surface_albedo 




SUBROUTINE surface_albedosn(ssoil, veg, met, soil)
   
   USE define_types
   USE define_dimensions
   USE cable_common_module
   USE physical_constants
   
   TYPE (soil_snow_type),INTENT(INOUT) :: ssoil
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
    
   WHERE( veg%iveg == 16 )                                                     &
      soil%albsoilf = -0.022*( MIN( 275., MAX( 260., met%tk ) ) - 260. ) + 0.45

   WHERE(ssoil%snowd > 1. .and. veg%iveg == 16 ) soil%albsoilf = 0.85

   sfact = 0.68
  
   WHERE (soil%albsoilf <= 0.14)
      sfact = 0.5
   ELSEWHERE (soil%albsoilf > 0.14 .and. soil%albsoilf <= 0.20)
      sfact = 0.62
   END WHERE

   ssoil%albsoilsn(:,2) = 2. * soil%albsoilf / (1. + sfact)
   ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
  
   snrat=0.
   alir =0.
   alv  =0.

   WHERE ( ssoil%snowd > 1. .AND. .NOT. cable_runtime%um_radiation ) 
       
      ! new snow (cm H2O)
      dnsnow = MIN ( 1., .1 * MAX( 0., ssoil%snowd - ssoil%osnowd ) ) 
      
      ! Snow age depends on snow crystal growth, freezing of melt water,
      ! accumulation of dirt and amount of new snow.
      tmp = ssoil%isflag * ssoil%tggsn(:,1) + ( 1 - ssoil%isflag )            &
            * ssoil%tgg(:,1)
      tmp = MIN( tmp, tfrz )
      ar1 = 5000. * (1. / 273.15 - 1. / tmp) ! crystal growth  (-ve)
      ar2 = 10. * ar1 ! freezing of melt water
      snr = ssoil%snowd / max (ssoil%ssdnn, 200.)
      
      WHERE (soil%isoilm == 9)
         
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
      
      WHERE (ssoil%snowd <= 1.0)
         ssoil%snage = 0.
      ELSEWHERE
         ssoil%snage = max (0.,(ssoil%snage+dtau)*(1.-dnsnow))
      END WHERE
      
      fage = 1. - 1. / (1. + ssoil%snage ) !age factor

      tmp = MAX( .17365, met%coszen )
      fzenm = MAX( 0.0, MERGE( 0.0,                                           &
              ( 1. + 1./2. ) / ( 1. + 2.*2. * tmp ) - 1./2., tmp > 0.5 ) )

      tmp = alvo * (1.0 - 0.2 * fage)
      alv = .4 * fzenm * (1. - tmp) + tmp
      tmp = aliro * (1. - .5 * fage)

      ! use dry snow albedo for pernament land ice
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
   WHERE (ssoil%snowd > 1 .and. cable_runtime%um_radiation )
      
      snr = ssoil%snowd / MAX (ssoil%ssdnn, 200.)
      
      WHERE (soil%isoilm == 9)
         snrat = 1.
      ELSEWHERE
         snrat = MIN (1., snr / (snr + .1) )
      END WHERE
      
      fage = 1. - 1. / (1. + ssoil%snage ) !age factor
      tmp = MAX (.17365, met%coszen )
      fzenm = MAX( 0., MERGE( 0.0,                                             &
              ( 1. + 1./2. ) / ( 1. + 2. * 2. * tmp ) - 1./2., tmp > 0.5 ) )

       tmp = alvo * (1.0 - 0.2 * fage)
       alv = .4 * fzenm * (1. - tmp) + tmp
       tmp = aliro * (1. - .5 * fage)
      
      ! use dry snow albedo
      WHERE (soil%isoilm == 9)          

         tmp = 0.95 * (1.0 - 0.2 * fage)
         alv = .4 * fzenm * (1. - tmp) + tmp
         tmp = 0.75 * (1. - .5 * fage)
      
      END WHERE
      
      alir = .4 * fzenm * (1.0 - tmp) + tmp
      talb = .5 * (alv + alir) ! snow albedo
    
   ENDWHERE        ! snowd > 0
   

   ssoil%albsoilsn(:,2) = MIN( aliro,                                          &
                          ( 1. - snrat ) * ssoil%albsoilsn(:,2) + snrat * alir)
   
   ssoil%albsoilsn(:,1) = MIN( alvo,                                           &
                          ( 1. - snrat ) * ssoil%albsoilsn(:,1) + snrat * alv )

END SUBROUTINE surface_albedosn


!jhan:subr was reintroduced here to temporarily resolve issue when 
!creating libcable.a 
SUBROUTINE calc_rhoch(veg,c1,rhoch) 

   USE define_types
   USE other_constants
   TYPE (veg_parameter_type), INTENT(INOUT) :: veg
   REAL, INTENT(INOUT), DIMENSION(:,:) :: c1, rhoch

   c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
   c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
   c1(:,3) = 1.
    
   ! Canopy reflection black horiz leaves 
   ! (eq. 6.19 in Goudriaan and van Laar, 1994):
   rhoch = (1.0 - c1) / (1.0 + c1)

END SUBROUTINE calc_rhoch 


END MODULE albedo_module
