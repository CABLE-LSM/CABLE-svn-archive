MODULE cbl_snow_albedo_module

   USE cable_data_module, ONLY : ialbedo_type, point2constants 
   
  IMPLICIT NONE

  PUBLIC surface_albedosn
  PRIVATE

   TYPE(ialbedo_type) :: C

CONTAINS

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

   CALL point2constants(C) 
   
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



END MODULE cbl_snow_albedo_module
