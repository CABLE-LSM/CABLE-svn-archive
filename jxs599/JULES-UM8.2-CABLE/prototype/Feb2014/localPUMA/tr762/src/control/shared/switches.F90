! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Module containing switches imported directly by JULES subroutines
! (as opposed to receiving the switches as arguments)
!

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE switches

IMPLICIT NONE

  LOGICAL ::                                                      &
   l_point_data   = .FALSE.                                       &
                     ! Switch for using point rainfall data
  ,l_spec_albedo  = .FALSE.                                       & 
                     ! Switch spectrally varying land albedo 
  ,l_spec_alb_bs  = .FALSE.                                       & 
                     ! Switch to have only a bluesky albedo 
                     ! when using spectrally varying albedo
  ,l_snow_albedo  = .FALSE.                                       &
                     ! Switch for prognostic snow albedo (on land)
  ,l_phenol       = .FALSE.                                       &
                     ! Switch for phenology
  ,l_triffid      = .FALSE.                                       &
                     ! Switch for TRIFFID
  ,l_vg_soil      = .FALSE.                                       &
                     ! Switch for using Van Genuchten soil scheme
  ,l_aggregate    = .FALSE.                                       &
                     ! Switch for setting an aggregate surface
                     ! scheme
  ,l_360          = .FALSE.                                       &
                     ! Switch for setting a 360 day year
  ,l_um_jules     = .FALSE.                                       &
                     ! Switch for using JULES in the UM
  ,l_flake_model  = .FALSE.                                       &
                     ! Switch for using the FLake lake model
  ,l_o3_damage    = .FALSE.                                       &
                     ! Switch for ozone damage
  ,l_veg_compete  = .TRUE.                                        &
                     ! Switch for competing vegetation
                     ! Setting l_triffid = .TRUE. and this as
                     ! .FALSE. means that the carbon pools evolve
                     ! but the PFT distribution does not change
                     ! The default of .TRUE. means that enabling
                     ! TRIFFID has competing veg on by default
#if defined(UM_JULES)
  ,l_epot_corr    = .FALSE.                                       &
                     ! Switch for using a correction to the
                     ! calculation of potential evaporation
                     ! Default for UM is FALSE
#else
  ,l_epot_corr    = .TRUE.                                        &
                     ! Switch for using a correction to the
                     ! calculation of potential evaporation
                     ! Default for standalone is TRUE
#endif
  ,l_snowdep_surf = .FALSE.                                       &
                     ! use equivalent canopy snow depth for surface
                     ! calculations on tiles with a snow canopy
  ,l_land_ice_imp = .FALSE.                                       &
                     ! use implicit numerics to update land ice
                     ! temperatures
  ,l_tstar_sice_new = .FALSE.                                     &
                     ! calculate sea ice surface temperature in scheme
                     ! compatible with multi-categories. (This logical
                     ! is only used in a single category run.)
  ,l_rho_snow_corr  = .TRUE.                                      &
                     ! Switch for using a correction to the density
                     ! of the snow pack when nsnow=0 when
                     ! relayering in the new snow scheme
                     ! Has no effect for nsmax < 1
  ,l_baseflow_corr  = .TRUE.                                      &
                     ! Switch for using a correction to the
                     ! calculation of baseflow
                     ! Only used if l_top = T
  ,l_dpsids_dsdz    = .FALSE.                                     &
                     ! Switch to calculate vertical gradient of
                     ! soil suction with the assumption of
                     ! linearity only for fractional saturation
                     ! (consistent with the calculation of hydraulic
                     ! conductivity)
  ,l_albedo_obs      = .FALSE.                                    &
                     ! scale the albedo on tiles to agree with obs 
  ,l_dolr_land_black = .TRUE.                                     &
                     ! Do not use the surface emissivity in
                     ! adjusting the OLR at land points. 
                     ! This flag is introduced for historical 
                     !compatibility only. There is no
                     ! equivalent choice at sea points.
  ,l_bvoc_emis       = .FALSE.                                    &
                     ! Switch to enable calculation of BVOC emissions
  ,l_spec_sea_alb    = .FALSE.                                    &
                     ! Switch spectrally varying open sea albedo
  ,l_sea_alb_var_chl = .FALSE.                                    &
                     ! Switch varying chlorophyll in open sea albedo
  ,l_ssice_albedo    = .FALSE.                                    &
                     ! Switch for including the effect of
                     ! snow on the sea-ice albedo
  ,l_sice_scattering = .FALSE.                                    &
                     ! Switch for seaice albedo internal scatter
  ,l_sice_meltponds  = .FALSE.                                    &
                     ! Sea-ice albedo affected by meltponds
  ,l_sice_hadgem1a   = .FALSE.                                    &
                     ! HadGEM1 sea-ice albedo bug corrected
  ,l_sice_multilayers= .FALSE.                                    &
                     ! True if coupled to sea ice multilayer model
  ,l_cice_alb        = .FALSE.                                    &
                     ! T = use sea ice albedo scheme from the CICE model
                     ! The sea ice radiation code in control.F90 assumes
                     ! this is always FALSE (standalone JULES only)
  ,l_sice_heatflux   = .FALSE.                                    &
                     ! T: semi-implicit update of TI
  ,l_soil_sat_down     =.FALSE.                                   &
                     ! Direction of super_saturated soil moisture
                     ! TRUE excess water is pushed down
                     ! FALSE excess water is pushed up (as in JULES2.0)
  ,l_top         = .FALSE.                                        &
                     ! Switch for TOPMODEL-based hydrology
  ,l_pdm         = .FALSE.                                        &
                     ! Switch for PDM hydrology
  ,l_anthrop_heat_src  = .FALSE.                                  &
                     ! Switch for anthropogenic heat source on urban
                     ! tile
  ,l_ctile             = .FALSE.                                  &
                     ! True if coastal tiling
  ,l_hydrology         = .FALSE. 
                     ! Turns off hydrology code for UM_JULES

#if defined(UM_JULES)
  LOGICAL :: l_neg_tstar = .FALSE.  ! Test for negative surface temperature.
#endif

  INTEGER ::                                                      &
   can_model           = 4                                        &
!                            switch for thermal vegetation
  ,soilhc_method       = 1                                        &
!                            switch for the calculation method
!                            of soil thermal conductivity
!        SOILHC_METHOD=1: Method of Cox et al (1999).
!        SOILHC_METHOD=2: Simplified Johansen (1975).
  ,i_modiscopt         = 0                                        &
!                            Method of discretization
!                            in the surface layer
  ,frac_snow_subl_melt = 0                                        &
!                            switch for use of snow-cover
!                            fraction in the calculation of
!                            sublimation and melting
!                            0 = off
!                            1 = on
  ,all_tiles           = 0                                        &
!                            switch for doing calculations
!                            of tile properties on all tiles
!                            for all gridpoints even when the
!                            tile fraction is zero
!                            (except for land ice).
  ,can_rad_mod         = 4                                        &
!                            Canopy radiation model
  ,cor_mo_iter         = 1                                        &
!                            Switch for MO iteration correction
  ,iseaz0t             = 0                                        &
!                            Switch for the definition of
!                            the thermal roughness length over the sea.
  ,buddy_sea           = 0                                        &
!                            Switch to use the wind speed from
!                            adjacent sea points for the
!                            sea part of coastal grid points
  ,iscrntdiag          = 0                                        &
!                            Method of diagnosing the screen temperature
  ,i_aggregate_opt     = 0                                        &
!                            Method of aggregating tiled properties 
!                            ! i_aggregate_opt=0: Original option 
!                            ! i_aggregate_opt=1: Separate aggregation 
!                            !                    of z0h 
  ,i_sea_alb_method    = 1
!                            Method of diagnosing the Ocean Surface Albedo
!                            1 - Briegleb and Ramanathan, 1982, J. Appl. Met.
!                 (doi:10.1175/1520-0450(1982)021<1160:SADVIC>2.0.CO;2)
!                            2 - Modified Barker and Li, 1995, J. Climate,
!                 (doi:10.1175/1520-0442(1995)008<2213:ISOCSS>2.0.CO;2)
!                            3 - Jin et al. 2011, Optics Express
!                            (doi:10.1364/OE.19.026429)

  REAL ::                                                         &
    dz_pdm = 1.0                                                  &
!                            Soil layer thickness for PDM (m):
   ,b_pdm  = 1.0
!                            Shape factor for PDM:

!----------------------------------------------------------------
! Switch for IMOGEN (never changed from default in the UM)
!----------------------------------------------------------------
  LOGICAL ::                                                     &
   l_imogen = .FALSE.

  INTEGER ::                                                      &
   ilayers = 10      ! Number of layers for canopy radiation model

!----------------------------------------------------------------
! Options moved here from run_bl / bl_option_mod
!----------------------------------------------------------------
! Options for form drag
  INTEGER, PARAMETER :: no_drag         = 0
  INTEGER, PARAMETER :: effective_z0    = 1
  INTEGER, PARAMETER :: explicit_stress = 2
! Switch for orographic form drag
  INTEGER :: formdrag = no_drag

! Drag coefficient for orographic form drag
  REAL :: orog_drag_param = 0.3

! Switch to implement stability dependence of orographic form drag
  INTEGER :: fd_stab_dep = 0

! Switch to include the effect of convective downdraughts on surface exchange
  INTEGER :: ISrfExCnvGust = 0
! OFF (=0) => not used: only boundary-layer gustiness is considered
! (original version)
! IP_SrfExWithCnv (=1) the impact of gustiness due to boundary layer eddies
! is reduced relative to the above, but eddies driven by convective
! downdraughts are included
  INTEGER, PARAMETER :: IP_SrfExWithCnv = 1

  INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
  INTEGER, PARAMETER :: Limit_ObukhovL = 3
!       Option under the COR_MO_ITER switch for imposing
!       lower limit on L in very stable conditions.

! Options for iseaz0t
! Standard flixed value of thermal roughness length over sea
  INTEGER, PARAMETER :: Fixed_Z0T = 0
! Thermal roughness length over sea defined from surface divergence theory
  INTEGER, PARAMETER :: SurfDivZ0T = 1

#if !defined(UM_JULES)
!-----------------------------------------------------------------------
! START OF NON-UM SECTION
!
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------

  LOGICAL ::                                                      &
   route         = .FALSE.                                        &
                !  Switch for runoff routing.
  ,routeonly     = .FALSE.                                        &
                !  Switch to only do runoff routing, nothing else
  ,l_trif_eq     = .FALSE.                                        &
                ! Switch for vegetation equilibrium
  ,l_neg_tstar   = .TRUE.                                         &
                ! Switch for -ve TSTAR error check
  ,l_cosz        = .TRUE.
                ! Switch for turning on calculations of cosz

!-----------------------------------------------------------------------
! Switches introduced during reconciliation with UM that can be set
! via the control file but are passed down in the UM
!-----------------------------------------------------------------------
  LOGICAL ::                                                      &
   l_q10               = .TRUE.                                  
                ! True if using Q10 for soil resp
  

!-----------------------------------------------------------------------
! Switches that are always .FALSE. standalone
!-----------------------------------------------------------------------
  LOGICAL ::                                                      &
   lq_mix_bl           = .FALSE.                                  &
                ! TRUE if mixing ratios used in
                ! boundary layer code
   ,l_spec_z0           = .FALSE.                                  &
                ! T if using prescribed sea surface
                ! roughness lengths
   ,l_co2_interactive   = .FALSE.                                  &
                ! Switch for 3D CO2 field
  ,l_dust              = .FALSE.                                  &
                ! Switch for mineral dust
  ,l_z0_orog           = .FALSE.                                  &
                ! T to use orog.roughness in surface calcs
  ,l_inland            = .FALSE.
                ! True if re-routing inland basin flow
                ! to soil moisture

!-----------------------------------------------------------------------
! END OF NON-UM SECTION
!-----------------------------------------------------------------------
#endif

  LOGICAL :: l_cable = .FALSE.  ! runtime switch to call CABLE and switch
                                ! components of JULES on/off  
!-----------------------------------------------------------------------
! Set up a namelist to allow switches to be set
! UM and standalone JULES currently have different requirements since
! some JULES options are defined higher up in the UM, and some options
! available in the UM are not subject to change standalone
!-----------------------------------------------------------------------
#if defined(UM_JULES)
  NAMELIST /jules_switches/ l_point_data, l_spec_albedo,          &
                            l_snow_albedo,l_phenol, l_triffid,    &
                            l_vg_soil, l_epot_corr,               & 
                            l_aggregate, i_aggregate_opt,         & 
                            l_snowdep_surf,                       &
                            l_land_ice_imp,                       &
                            l_flake_model,                        &
                            l_tstar_sice_new,                     &
                            can_model,soilhc_method,              &
                            i_modiscopt,                          &
                            frac_snow_subl_melt,all_tiles,        &
                            can_rad_mod, cor_mo_iter,             &
                            iseaz0t, buddy_sea,                   &
                            iscrntdiag, b_pdm, dz_pdm,            &
                            l_rho_snow_corr,l_baseflow_corr,      &
                            l_dpsids_dsdz, l_albedo_obs,          &
                            l_dolr_land_black, l_bvoc_emis,       &
                            l_ssice_albedo, l_sice_meltponds,     &
                            l_sice_scattering, l_sice_hadgem1a,   &
                            l_sice_multilayers, l_cice_alb,       &
                            l_sice_heatflux, i_sea_alb_method,    &
                            l_top, l_pdm, l_soil_sat_down,        &
                            l_anthrop_heat_src, l_ctile, ilayers, &
                            l_spec_alb_bs,                        &
                            l_spec_sea_alb, l_sea_alb_var_chl,    &
                            formdrag, orog_drag_param,            &
                            fd_stab_dep, ISrfExCnvGust,           &
                            l_hydrology, l_cable
#else
  NAMELIST /jules_switches/ l_aggregate,i_aggregate_opt,          & 
                            can_model,can_rad_mod,                & 
                            ilayers,l_cosz,l_spec_albedo,         &
                            l_snow_albedo,l_phenol,l_triffid,     &
                            l_veg_compete,l_trif_eq,l_top,l_pdm,  &
                            l_anthrop_heat_src,l_o3_damage,       &
                            l_imogen,l_epot_corr,l_snowdep_surf,  &
                            l_tstar_sice_new, l_land_ice_imp,     &
                            iscrntdiag,l_360,l_q10,               &
                            l_soil_sat_down,l_vg_soil,            &
                            soilhc_method,l_point_data,           &
                            l_rho_snow_corr,l_baseflow_corr,      &
                            l_dpsids_dsdz, l_albedo_obs,          &
                            l_dolr_land_black, l_bvoc_emis,       &
                            l_spec_alb_bs, i_sea_alb_method,      &
                            l_spec_sea_alb, l_sea_alb_var_chl,    &
                            l_cable
#endif

CONTAINS

  SUBROUTINE print_nlist_jules_switches()
    USE jules_print_mgr, ONLY : jules_print
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer

    CALL jules_print('switches', &
        'Contents of namelist jules_switches')

#if defined(UM_JULES)
    WRITE(lineBuffer,*)' l_point_data = ',l_point_data
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_snow_albedo = ',l_snow_albedo
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_phenol = ',l_phenol
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_triffid = ',l_triffid
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_vg_soil = ',l_vg_soil
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_aggregate = ',l_aggregate
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_epot_corr = ',l_epot_corr
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_snowdep_surf = ',l_snowdep_surf
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_flake_model = ',l_flake_model
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' can_model = ',can_model
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' soilhc_method = ',soilhc_method
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' i_modiscopt = ',i_modiscopt
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' frac_snow_subl_melt = ',frac_snow_subl_melt
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' all_tiles = ',all_tiles
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' can_rad_mod = ',can_rad_mod
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' cor_mo_iter = ',cor_mo_iter
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' iseaz0t = ',iseaz0t
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' buddy_sea = ',buddy_sea
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' iscrntdiag = ',iscrntdiag
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' b_pdm = ',b_pdm
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' dz_pdm = ',dz_pdm
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_rho_snow_corr = ',l_rho_snow_corr
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_baseflow_corr = ',l_baseflow_corr
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' i_sea_alb_method = ',i_sea_alb_method
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)'l_spec_sea_alb  = ',l_spec_sea_alb
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_sea_alb_var_chl = ',l_sea_alb_var_chl
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' formdrag = ',formdrag
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' orog_drag_param = ',orog_drag_param
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' fd_stab_dep = ',fd_stab_dep
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' isrfexcnvgust = ',isrfexcnvgust
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_hydrology = ',l_hydrology
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_cable = ',l_cable
    CALL jules_print('switches',lineBuffer)
#else
    WRITE(lineBuffer,*)' l_aggregate = ',l_aggregate
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' can_model = ',can_model
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' can_rad_mod = ',can_rad_mod
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' ilayers = ',ilayers
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_cosz = ',l_cosz
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_spec_albedo = ',l_spec_albedo
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_phenol = ',l_phenol
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_triffid = ',l_triffid
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_veg_compete = ',l_veg_compete
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_trif_eq = ',l_trif_eq
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_top = ',l_top
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_pdm = ',l_pdm
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_anthrop_heat_src = ',l_anthrop_heat_src
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_o3_damage = ',l_o3_damage
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_imogen = ',l_imogen
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_epot_corr = ',l_epot_corr
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_snowdep_surf = ',l_snowdep_surf
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' iscrntdiag = ',iscrntdiag
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_360 = ',l_360
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_q10 = ',l_q10
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_soil_sat_down = ',l_soil_sat_down
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_vg_soil = ',l_vg_soil
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' soilhc_method = ',soilhc_method
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_point_data = ',l_point_data
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_rho_snow_corr = ',l_rho_snow_corr
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_baseflow_corr = ',l_baseflow_corr
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' i_sea_alb_method = ',i_sea_alb_method
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)'l_spec_sea_alb  = ',l_spec_sea_alb
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_sea_alb_var_chl = ',l_sea_alb_var_chl
    CALL jules_print('switches',lineBuffer)
    WRITE(lineBuffer,*)' l_cable = ',l_cable
    CALL jules_print('switches',lineBuffer)
#endif
    CALL jules_print('switches', &
        '- - - - - - end of namelist - - - - - -')

  END SUBROUTINE print_nlist_jules_switches

END MODULE switches
