
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
    USE cable_common_module
    USE carbon_module
    USE soil_snow_module
    USE define_types
    USE physical_constants
    USE roughness_module
    USE radiation_module
    USE air_module
      use cable_diag_module, only : cable_stat
      USE casadimension, only : icycle !declared in casadimension which is used in casa_cnp
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
    REAL(r_1), DIMENSION(mp)            :: xx1


      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('cbm')

      !jhan:do this in iniotialization   
       ! Fix in-canopy turbulence scheme globally:
       veg%meth = 1

      if( cable_runtime%um ) then
         cable_runtime%um_radiation = .false.
         if( cable_runtime%um_explicit ) then
            call ruff_resist(veg, rough, ssoil, soil, met, canopy)
            met%tk = met%tk + grav/capp*(rough%zref_tq + 0.9*rough%z0m)
            !jhan:consult group b4 rm
            met%tc = met%tk - tfrz
         endif
         CALL define_air (met, air)
      else
         call ruff_resist(veg, rough, ssoil, soil, met, canopy)
      endif


      CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt

      if( cable_runtime%um ) then
         if( cable_runtime%um_explicit ) then
            !call surface_albedo(ssoil, veg, met, rad, soil, canopy)
            call surface_albedo(ssoil, veg, met, rad, soil, canopy)
         endif
      else
         call surface_albedo(ssoil, veg, met, rad, soil, canopy)
      endif
    
    ! Calculate canopy variables:
    CALL define_canopy(bal,rad,rough,air,met,dels,ssoil,soil,veg, canopy)

   !canopy%fes_cor = 0.
   !canopy%fhs_cor = 0.
   ssoil%otss_0 = ssoil%otss
   ssoil%otss = ssoil%tss
   if( cable_runtime%um ) then
      ssoil%owetfac = ssoil%wetfac
      if( cable_runtime%um_hydrology ) then
         call soil_snow(dels, soil, ssoil, canopy, met, bal,veg)
      endif
    else
      call soil_snow(dels, soil, ssoil, canopy, met, bal,veg)
    endif

    !jhan:Eva addded this
    ssoil%deltss = ssoil%tss-ssoil%otss
    canopy%fhs = canopy%fhs + (ssoil%tss-ssoil%otss)*ssoil%dfh_dtg
    canopy%fhs_cor = canopy%fhs_cor + (ssoil%tss-ssoil%otss)*ssoil%dfh_dtg
    canopy%fes = canopy%fes + (ssoil%tss-ssoil%otss)*(ssoil%cls*ssoil%dfe_ddq*ssoil%ddq_dtg)
    canopy%fes_cor = canopy%fes_cor + (ssoil%tss-ssoil%otss)*(ssoil%cls*ssoil%dfe_ddq*ssoil%ddq_dtg)

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
    xx1 =  rad%trad
    rad%trad = (1.-rad%transd)*canopy%tv + rad%transd * ssoil%tss
 
    ! rml 17/1/11 move all plant resp and soil resp calculations here from canopy
    ! and only call on implicit step
    ! put old and new soil resp calculations into soilcarb subroutine
    !make new plantcarb subroutine
    if (.not.cable_runtime%um_explicit .AND. icycle == 0) then
        !calculate canopy%frp
        CALL plantcarb(veg,bgc,met,canopy)
        !calculate canopy%frs
        CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)

        ! rml - why was this here, is it needed??
        !rad%flws = sboltz*emsoil* ssoil%tss **4
        CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)

        canopy%fnpp = -1.0* canopy%fpn - canopy%frp
        canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

    endif

  
END SUBROUTINE cbm

END MODULE cbm_module
