
MODULE cbm_module
  USE canopy_module
  USE albedo_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC cbm 
CONTAINS
  SUBROUTINE cbm( dels, air, bgc, canopy, met, bal, &
             rad, rough, soil, ssoil, sum_flux, veg )
   use define_dimensions
    USE carbon_module
    USE soil_snow_module
    USE define_types
    USE physical_constants
    USE roughness_module
    USE radiation_module
    USE air_module
      use cable_common_module
      use cable_diag_module, only : cable_stat
      use cable_data_module, only : model_type, const 
      !USE casadimension, only : icycle !declared in casadimension which is used in casa_cnp
      implicit none
    INTEGER(i_d) :: k,kk,j  
    INTEGER(i_d) :: idjd,idjd1,idjd2,idjd3  
    REAL(r_1), INTENT(IN)               :: dels ! time setp size (s)
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (bgc_pool_type), INTENT(INOUT) :: bgc
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (balances_type), INTENT(INOUT)         :: bal
    TYPE (radiation_type), INTENT(INOUT)        :: rad
    TYPE (roughness_type), INTENT(INOUT)        :: rough
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
    TYPE (soil_snow_type), INTENT(INOUT)        :: ssoil
    TYPE (sum_flux_type), INTENT(INOUT) :: sum_flux
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg  
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints

      !jhan: we can put all this data strucutre stuff into cable_data
      !alos we can put define dimensions into constants% type
      TYPE (model_type) :: cable 
      

      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('cbm')

      !jhan:do this in iniotialization   
       ! Fix in-canopy turbulence scheme globally:
       veg%meth = 1

      call cable_data_alloc( cable, air, bgc, canopy, met, bal, &
                  rad, rough, soil, ssoil, sum_flux, veg )
      

      if( cable_runtime%um ) then
         cable_runtime%um_radiation = .false.
         if( cable_runtime%um_explicit ) then
            !jhan: commented as arg list tested in else
            !call ruff_resist(veg, rough, ssoil, soil, met, canopy)
            met%tk = met%tk + grav/capp*(rough%zref_tq + 0.9*rough%z0m)
            !jhan:consult group b4 rm
            met%tc = met%tk - tfrz
         endif
         CALL define_air ( cable% a% air, & 
                           cable% i% met% met_tvair, &
                           cable% i% met% met_pmb, &
                           const%phys )
         !CALL define_air (met, air)
      else
         call ruff_resist( cable% i% rough, &
                           cable% o% rough, &
                           cable% a% rough, &
                           cable% i% ssnow% ssnow_snowd, &
                           cable% i% ssnow% ssnow_ssdnn, & 
                           cable% i% veg% veg_hc, & 
                           cable% i% veg% veg_vlai, & 
                           cable% i% veg% veg_iveg, & 
                           cable% a% canopy% canopy_vlaiw, &
                           cable% a% canopy% canopy_rghlai, &
                           const%phys & 
                         )
      endif


      CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt

      if( cable_runtime%um ) then
         if( cable_runtime%um_explicit ) then
            !call surface_albedo(ssoil, veg, met, rad, soil, canopy)
            !jhan: commented as arg list tested in else
            !call surface_albedo(ssoil, veg, met, rad, soil, canopy)
         endif
      else
         call surface_albedo(ssoil, veg, met, rad, soil, canopy)
      endif
!      print *,'cable_runtime%um,cable_runtime%um_explicit',ktau_gl,cable_runtime%um,cable_runtime%um_explicit, &
!      cable_runtime%um_hydrology
    
!      IF( L_EXPLICIT )  then
!        print 101,ktau_gl,rad%albedo(idjd1,1),rad%albedo(idjd1,2),met%coszen(idjd1), &
!                        met%fsd(idjd1,1)+met%fsd(idjd1,2),ssoil%albsoilsn(idjd1,:),ssoil%snage(idjd1), &
!                        veg%vlai(idjd1),rough%hruff(idjd1),veg%hc(idjd1),ssoil%dfh_dtg(idjd1), &
!                        ssoil%cls(idjd1)*ssoil%dfe_ddq(idjd1)*ssoil%ddq_dtg(idjd1),soil%albsoil(idjd1,1), &
!                        ssoil%tss(idjd1),met%tk(idjd1),met%tvair(idjd1),canopy%tv(idjd1)
!    101 format(1x,'idjd1ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!        print 102,ktau_gl,rad%albedo(idjd2,1),rad%albedo(idjd2,2),met%coszen(idjd2), &
!                        met%fsd(idjd2,1)+met%fsd(idjd2,2),ssoil%albsoilsn(idjd2,:),ssoil%snage(idjd2), &
!                        veg%vlai(idjd2),rough%hruff(idjd2),veg%hc(idjd2),ssoil%dfh_dtg(idjd2), &
!                        ssoil%cls(idjd2)*ssoil%dfe_ddq(idjd2)*ssoil%ddq_dtg(idjd2),soil%albsoil(idjd2,1), &
!                        ssoil%tss(idjd2),met%tk(idjd2),met%tvair(idjd2),canopy%tv(idjd2)
!    102 format(1x,'idjd2ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!        print 103,ktau_gl,rad%albedo(idjd3,1),rad%albedo(idjd3,2),met%coszen(idjd3), &
!                        met%fsd(idjd3,1)+met%fsd(idjd3,2),ssoil%albsoilsn(idjd3,:),ssoil%snage(idjd3), &
!                        veg%vlai(idjd3),rough%hruff(idjd3),veg%hc(idjd3),ssoil%dfh_dtg(idjd3), &
!                        ssoil%cls(idjd3)*ssoil%dfe_ddq(idjd3)*ssoil%ddq_dtg(idjd3),soil%albsoil(idjd3,1), &
!                        ssoil%tss(idjd3),met%tk(idjd3),met%tvair(idjd3),canopy%tv(idjd3)
!    103 format(1x,'idjd3ALBEDO',i5,3f6.3,1x,f6.0,3f5.2,x,f6.3,x,f5.2,f10.4,f7.3,2f7.3,f6.2,4f6.1)
!       ENDIF

    ! Calculate canopy variables:
    CALL define_canopy(bal,rad,rough,air,met,dels,ssoil,soil,veg, canopy, &
         cable%a%air, & 
                           cable%i%met% met_tvair, &
                           cable%i%met% met_pmb )

   ssoil%otss = ssoil%tss
   if( cable_runtime%um ) then
      if( cable_runtime%um_hydrology ) then
         ssoil%owetfac = ssoil%wetfac
         call soil_snow(dels, soil, ssoil, canopy, met, bal,veg)
      endif
    else
      call soil_snow(dels, soil, ssoil, canopy, met, bal,veg)
    endif

    !jhan:Eva addded this
    canopy%fhs = canopy%fhs + (ssoil%tss-ssoil%otss)*ssoil%dfh_dtg
    canopy%fes = canopy%fes + (ssoil%tss-ssoil%otss)*(ssoil%cls*ssoil%dfe_ddq*ssoil%ddq_dtg)

    canopy%fh = canopy%fhv + canopy%fhs

    !   need to adjust fe after soilsnow
    canopy%fev  = canopy%fevc + canopy%fevw
    ! Calculate total latent heat flux:
    canopy%fe = canopy%fev + canopy%fes
    ! Calculate net radiation absorbed by soil + veg
    canopy%rnet = canopy%fns + canopy%fnv

    ! Calculate radiative/skin temperature:
    rad%trad = ( (1.-rad%transd)*canopy%tv**4 + &
          rad%transd * ssoil%tss**4 )**0.25

! rml 17/1/11 move all plant resp and soil resp calculations here from canopy
! and only call on implicit step
! put old and new soil resp calculations into soilcarb subroutine
! make new plantcarb subroutine
   if( cable_runtime%um ) then
     if (.not.cable_runtime%um_explicit) then
!      calculate canopy%frp
       CALL plantcarb(veg,bgc,met,canopy)
!      calculate canopy%frs
       CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)

! rml - why was this here, is it needed??
!        rad%flws = sboltz*emsoil* ssoil%tss **4
       CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)

       canopy%fnpp = -1.0* canopy%fpn - canopy%frp
       canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

     endif
   endif

  
END SUBROUTINE cbm


   subroutine cable_data_alloc( cable, air, bgc, canopy, met, bal, &
                  rad, rough, soil, ssoil, sum_flux, veg )
      use define_dimensions
      use define_types
      use cable_data_module, only : model_type 
      !USE casadimension, only : icycle !declared in casadimension which is used in casa_cnp
      implicit none

      !alos we can put define dimensions into constants% type
      TYPE (model_type) :: cable 
      TYPE (air_type), INTENT(INOUT)      :: air
      TYPE (bgc_pool_type), INTENT(INOUT) :: bgc
      TYPE (canopy_type), INTENT(INOUT)   :: canopy
      TYPE (met_type), INTENT(INOUT)      :: met
      TYPE (balances_type), INTENT(INOUT)         :: bal
      TYPE (radiation_type), INTENT(INOUT)        :: rad
      TYPE (roughness_type), INTENT(INOUT)        :: rough
      TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
      TYPE (soil_snow_type), INTENT(INOUT)        :: ssoil
      TYPE (sum_flux_type), INTENT(INOUT) :: sum_flux
      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg  
      
      allocate( cable% a% air% air_rho(mp) ) 
      allocate( cable% a% air% air_rlam(mp) ) 
      allocate( cable% a% air% air_epsi(mp) )
      allocate( cable% a% air% air_visc(mp) )
      allocate( cable% a% air% air_psyc(mp) )
      allocate( cable% a% air% air_dsatdk(mp) ) 
      allocate( cable% a% air% air_cmolar(mp) ) 
 
      allocate( cable% i% met% met_tvair(mp) ) 
      allocate( cable% i% met% met_pmb(mp) ) 
    
      cable% a% air% air_rho => air%rho
      cable% a% air% air_rlam => air%rlam ! latent heat for water (j/kg)
      cable% a% air% air_epsi => air%epsi! d(qsat)/dT ((kg/kg)/K)
      cable% a% air% air_visc => air%visc ! air kinematic viscosity (m2/s)
      cable% a% air% air_psyc => air%psyc ! psychrometric constant
      cable% a% air% air_dsatdk => air%dsatdk! d(es)/dT (mb/K)
      cable% a% air% air_cmolar => air%cmolar ! conv. from m/s to mol/m2/s
 
      cable% i% met% met_tvair => met%tvair ! conv. from m/s to mol/m2/s
      cable% i% met% met_pmb => met%pmb ! conv. from m/s to mol/m2/s
    
      allocate( cable% i% rough% rough_hruff_grmx(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% i% rough% rough_za_uv(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% i% rough% rough_za_tq(mp) ) ! max ht of canopy from tiles on same grid 

      allocate( cable% o% rough% rough_z0m(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% o% rough% rough_zref_tq(mp) ) ! max ht of canopy from tiles on same grid 
   
      allocate( cable% a% rough% rough_hruff(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_coexp(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_disp(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_rt0us(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_rt1usa(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_rt1usb(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_rt1(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_usuh(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_zref_uv(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_zruffs(mp) ) ! max ht of canopy from tiles on same grid 
      allocate( cable% a% rough% rough_z0soilsn(mp) ) ! max ht of canopy from tiles on same grid 


      cable% i% rough% rough_hruff_grmx => rough%hruff_grmx
      cable% i% rough% rough_za_uv => rough%za_uv
      cable% i% rough% rough_za_tq => rough%za_tq

      cable% o% rough% rough_z0m =>rough%z0m
      cable% o% rough% rough_zref_tq => rough%zref_tq

      cable% a% rough% rough_hruff => rough%hruff
      cable% a% rough% rough_coexp => rough%coexp
      cable% a% rough% rough_disp => rough%disp
      cable% a% rough% rough_rt0us => rough%rt0us
      cable% a% rough% rough_rt1usa => rough%rt1usa
      cable% a% rough% rough_rt1usb => rough%rt1usb
      cable% a% rough% rough_rt1 => rough%rt1
      cable% a% rough% rough_usuh => rough%usuh
      cable% a% rough% rough_zref_uv => rough%zref_uv
      cable% a% rough% rough_zruffs => rough%zruffs
      cable% a% rough% rough_z0soilsn => rough%z0soilsn


      allocate( cable% i% ssnow% ssnow_snowd(mp) ) 
      allocate( cable% i% ssnow% ssnow_ssdnn(mp) ) 

      cable% i% ssnow% ssnow_snowd => ssoil%snowd 
      cable% i% ssnow% ssnow_ssdnn => ssoil%ssdnn 


      allocate( cable% i% veg% veg_hc(mp) ) 
      allocate( cable% i% veg% veg_vlai(mp) ) 
      allocate( cable% i% veg% veg_iveg(mp) ) 

      cable% i% veg% veg_hc => veg%hc 
      cable% i% veg% veg_vlai => veg%vlai
      cable% i% veg% veg_iveg  => veg%iveg


      allocate( cable% a% canopy% canopy_vlaiw (mp) ) 
      allocate( cable% a% canopy% canopy_rghlai (mp) ) 

      cable% a% canopy% canopy_vlaiw => canopy%vlaiw 
      cable% a% canopy% canopy_rghlai => canopy%rghlai 

 
      !allocate( cable% a% % (mp) ) 
      !cable% a% % => 
      return
   end subroutine cable_data_alloc
 


END MODULE cbm_module
