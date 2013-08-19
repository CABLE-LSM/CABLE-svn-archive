#define NO_CASA_YET 1

MODULE cable_cbm_module
   
   USE cable_canopy_module
   USE cable_albedo_module
  
   IMPLICIT NONE
  
   PRIVATE
   PUBLIC cbm 

CONTAINS

   SUBROUTINE cbm( dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg )
    
   USE cable_common_module
   USE cable_carbon_module
   USE cable_soil_snow_module
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cable_radiation_module
   USE cable_air_module
#ifndef NO_CASA_YET
   USE casadimension,     only : icycle ! used in casa_cnp
#endif
   USE cable_data_module, ONLY : icbm_type, point2constants 
   
   !ptrs to local constants 
   TYPE( icbm_type ) :: C
   ! CABLE model variables
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (bgc_pool_type),  INTENT(INOUT) :: bgc
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (balances_type),  INTENT(INOUT) :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (sum_flux_type),  INTENT(INOUT) :: sum_flux
    
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg  

   REAL, INTENT(IN)               :: dels ! time setp size (s)
    
   INTEGER :: k,kk,j  
   !INTEGER(i_d) :: idjd,idjd1,idjd2,idjd3,idjd4
   INTEGER :: idjd,idjd1,idjd2,idjd3,idjd4


#ifdef NO_CASA_YET
   INTEGER :: ICYCLE
   ICYCLE = 0
#endif

       idjd1 = 1036
       idjd2 = 9682
       idjd3 = 17386
       idjd4 = 21140

   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

   IF( cable_runtime%um ) THEN
      
      cable_runtime%um_radiation = .FALSE.
      
      IF( cable_runtime%um_explicit ) THEN
         CALL ruff_resist(veg, rough, ssnow, canopy, met)
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF
      
      CALL define_air (met, air)
   
   ELSE
      call ruff_resist(veg, rough, ssnow, canopy, met)
   ENDIF


   CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt

   IF( cable_runtime%um ) THEN
      
      IF( cable_runtime%um_explicit ) THEN
         CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
      ENDIF
   
   ELSE
      CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
   ENDIf
    
!      IF( L_EXPLICIT )  then
!        print 101,ktau_gl,rad%albedo(idjd1,1),rad%albedo(idjd1,2),met%coszen(idjd1), &
!                        met%fsd(idjd1,1)+met%fsd(idjd1,2),ssnow%albsoilsn(idjd1,:),ssnow%snage(idjd1), &
!                        veg%vlai(idjd1),rough%hruff(idjd1),veg%hc(idjd1),ssnow%dfh_dtg(idjd1), &
!                        ssnow%cls(idjd1)*ssnow%dfe_ddq(idjd1)*ssnow%ddq_dtg(idjd1),soil%albsoil(idjd1,1), &
!                        ssnow%tss(idjd1),met%tk(idjd1),met%tvair(idjd1),canopy%tv(idjd1)
!    101 format(1x,'idjd1ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!        print 102,ktau_gl,rad%albedo(idjd2,1),rad%albedo(idjd2,2),met%coszen(idjd2), &
!                        met%fsd(idjd2,1)+met%fsd(idjd2,2),ssnow%albsoilsn(idjd2,:),ssnow%snage(idjd2), &
!                        veg%vlai(idjd2),rough%hruff(idjd2),veg%hc(idjd2),ssnow%dfh_dtg(idjd2), &
!                        ssnow%cls(idjd2)*ssnow%dfe_ddq(idjd2)*ssnow%ddq_dtg(idjd2),soil%albsoil(idjd2,1), &
!                        ssnow%tss(idjd2),met%tk(idjd2),met%tvair(idjd2),canopy%tv(idjd2) 
!    102 format(1x,'idjd2ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!        print 103,ktau_gl,rad%albedo(idjd3,1),rad%albedo(idjd3,2),met%coszen(idjd3), &
!                        met%fsd(idjd3,1)+met%fsd(idjd3,2),ssnow%albsoilsn(idjd3,:),ssnow%snage(idjd3), &
!                        veg%vlai(idjd3),rough%hruff(idjd3),veg%hc(idjd3),ssnow%dfh_dtg(idjd3), &
!                        ssnow%cls(idjd3)*ssnow%dfe_ddq(idjd3)*ssnow%ddq_dtg(idjd3),soil%albsoil(idjd3,1), &
!                        ssnow%tss(idjd3),met%tk(idjd3),met%tvair(idjd3),canopy%tv(idjd3)
!    103 format(1x,'idjd3ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!        print 104,ktau_gl,rad%albedo(idjd4,1),rad%albedo(idjd4,2),met%coszen(idjd4), &
!                        met%fsd(idjd4,1)+met%fsd(idjd4,2),ssnow%albsoilsn(idjd4,:),ssnow%snage(idjd4), &
!                        veg%vlai(idjd4),rough%hruff(idjd4),veg%hc(idjd4),ssnow%dfh_dtg(idjd4), &
!                        ssnow%cls(idjd4)*ssnow%dfe_ddq(idjd4)*ssnow%ddq_dtg(idjd4),soil%albsoil(idjd4,1), &
!                        ssnow%tss(idjd4),met%tk(idjd4),met%tvair(idjd4),canopy%tv(idjd4)
!    104 format(1x,'idjd4ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!       ENDIF
!
!     print 201,ktau_gl,rad%albedo(idjd1,1),rad%albedo(idjd1,2),rad%reffdf(idjd1,1),rad%reffdf(idjd1,2), &
!                         rad%reffbm(idjd1,1),rad%reffbm(idjd1,2),soil%albsoil(idjd1,1),soil%albsoil(idjd1,2), &
!                          ssnow%albsoilsn(idjd1,:),ssnow%snowd(idjd1)
!201 format(1x,'atm_ph1albedo1',i6,11f6.3,2f7.1)
!     print 202,ktau_gl,rad%albedo(idjd2,1),rad%albedo(idjd2,2),rad%reffdf(idjd2,1),rad%reffdf(idjd2,2), &
!                         rad%reffbm(idjd2,1),rad%reffbm(idjd2,2),soil%albsoil(idjd2,1),soil%albsoil(idjd2,2), &
!                          ssnow%albsoilsn(idjd2,:),ssnow%snowd(idjd2)
!202 format(1x,'atm_ph1albedo2',i6,11f6.3,2f7.1)
!     print 203,ktau_gl,rad%albedo(idjd3,1),rad%albedo(idjd3,2),rad%reffdf(idjd3,1),rad%reffdf(idjd3,2), &
!                         rad%reffbm(idjd3,1),rad%reffbm(idjd3,2),soil%albsoil(idjd3,1),soil%albsoil(idjd3,2), &
!                          ssnow%albsoilsn(idjd3,:),ssnow%snowd(idjd3)
!203 format(1x,'atm_ph1albedo3',i6,11f6.3,2f7.1)
!     print 204,rad%albedo(idjd4,1),rad%albedo(idjd4,2),rad%reffdf(idjd4,1),rad%reffdf(idjd4,2), &
!                         rad%reffbm(idjd4,1),rad%reffbm(idjd4,2),soil%albsoil(idjd4,1),soil%albsoil(idjd4,2), &
!                          ssnow%albsoilsn(idjd4,:),ssnow%snowd(idjd4)
!204 format(1x,'atm_ph1albedo4',i6,11f6.3,2f7.1)
!
!



   ! Calculate canopy variables:
   CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy)

   ssnow%otss_0 = ssnow%otss
   ssnow%otss = ssnow%tss

   IF( cable_runtime%um ) THEN
      
      ssnow%owetfac = ssnow%wetfac
      IF( cable_runtime%um_hydrology ) THEN
         CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
      ENDIF

   ELSE
      call soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
   ENDIF

   ssnow%deltss = ssnow%tss-ssnow%otss

   canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
   
   canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
   
    canopy%fh = canopy%fhv + canopy%fhs

   canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) *                        &
                ( ssnow%dfe_ddq * ssnow%ddq_dtg )
                !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
   
   canopy%fes_cor = canopy%fes_cor + ( ssnow%tss-ssnow%otss ) *                &
                    ( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )

   ! need to adjust fe after soilsnow
   canopy%fev  = canopy%fevc + canopy%fevw
  
   ! Calculate total latent heat flux:
   canopy%fe = canopy%fev + canopy%fes

   ! Calculate net radiation absorbed by soil + veg
   canopy%rnet = canopy%fns + canopy%fnv

   ! Calculate radiative/skin temperature:
   rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
              rad%transd * ssnow%tss**4 )**0.25

   ! rml 17/1/11 move all plant resp and soil resp calculations here            
   ! from canopy. in UM only call on implicit step.
   ! put old and new soil resp calculations into soilcarb subroutine
   ! make new plantcarb subroutine
   IF (.not.cable_runtime%um_explicit .AND. icycle == 0) THEN

      !calculate canopy%frp
      CALL plantcarb(veg,bgc,met,canopy)
     
      !calculate canopy%frs
      CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

      CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

      canopy%fnpp = -1.0* canopy%fpn - canopy%frp
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

   ENDIF

  
END SUBROUTINE cbm

END MODULE cable_cbm_module


