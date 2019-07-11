#if !defined(UM_JULES)
  SUBROUTINE control ( a_step )

!-------------------------------------------------------------------------------
! Control level routine, to call the main parts of the model.
!-------------------------------------------------------------------------------
  USE datetime_mod, ONLY : SECS_IN_DAY, day_of_year

  USE model_time_mod, ONLY : timestep_len, current_time

  USE update_mod, ONLY : io_rad_type, use_diff_rad, diff_frac_const

  USE aero, ONLY :  &
!  imported scalars with intent(in)
     co2_mmr  &
!  imported arrays with intent(in)
    ,co2_3d  &
!  imported arrays with intent(out)
    ,aresist,aresist_tile,cd_std_dust,r_b_dust,resist_b         &
    ,resist_b_tile,rho_aresist,rho_aresist_tile,RHO_CD_MODV1,u_s_std_tile

  USE ancil_info, ONLY :  &
              co2_dim_len,co2_dim_row,frac,halo_i,halo_j                     &
             ,ice_fract,ice_fract_ncat,land_index,land_pts                   &
             ,land_mask,land_pts_trif,lice_index,lice_pts,n_rows,nice        &
             ,npft_trif,ntiles,off_x,off_y,row_length,rows                   &
             ,sm_levels,soil_index,soil_pts,tile_index,tile_pts              &
             ,z1_tq,z1_uv,nsmax,dim_cs1,dim_cs2,nice_use,sice_pts_ncat       &
             ,sice_index_ncat,ssi_index

  USE nstypes, ONLY : ntype

  USE bl_option_mod, ONLY : on

! Surface elevation
  USE c_elevate, ONLY :  &
!  imported arrays with intent(in)
     z_land

  USE c_epslon, ONLY :  &
!  imported scalar parameters
     c_virtual

  USE c_gamma, ONLY :  &
!  imported arrays with intent(in)
     gamma

  USE coastal, ONLY :  &
!  imported arrays with intent(in)
     fland,flandg  &
!  imported arrays with intent(out)
    ,surf_ht_flux_land,surf_ht_flux_sice,taux_land,taux_ssi,tauy_land  &
    ,tauy_ssi                                                          &
    ,tstar_land,tstar_sea,tstar_sice,tstar_sice_ncat,tstar_ssi         &
    ,vshr_land,vshr_ssi

  USE csigma, ONLY :  &
!  imported scalar parameters
      sbcon

  USE diag_swchs, ONLY : sfme,simlt,slh,smlt,sq1p5,st1p5  &
            ,stf_hf_snow_melt,stf_sub_surf_roff,su10,sv10,SZ0HEFF

  USE fluxes, ONLY :  &
!  imported arrays with intent(out)
     alb_tile,e_sea,ecan,ecan_tile,ei,ei_tile,esoil,esoil_tile  &
    ,ext,fqw_1,fqw_tile,FQW_ICE,fsmc,ftl_1,FTL_ICE,ftl_tile,H_SEA  &
    ,HF_SNOW_MELT,LAND_ALBEDO,LATENT_HEAT,LE_TILE,MELT_TILE     &
    ,SEA_ICE_HTF,SICE_MLT_HTF,SNOMLT_SUB_HTF,SNOMLT_SURF_HTF    &
    ,SNOW_MELT,SNOWMELT,SUB_SURF_ROFF,SURF_HT_FLUX,surf_htf_tile  &
    ,SURF_ROFF,RADNET_TILE,TAUX_1,TAUY_1                          &
    ,TOT_TFALL,TSTAR  &
    ,sw_tile,emis_tile,alb_sice

  USE forcing, ONLY :  &
     con_rain,con_snow,ls_rain,ls_snow,lw_down  &
    ,pstar,qw_1,sw_down,tl_1,u_0,u_1,v_0,v_1,diff_rad

  USE orog, ONLY :  &
!  imported arrays
     h_blend_orog,ho2r2_orog,sil_orog_land,z0m_eff

  USE p_s_parms, ONLY : albsoil,albobs_sw, albobs_vis, albobs_nir,b     &
       ,catch,catch_snow,cosz,hcap,hcon,infil_tile,satcon,sathh,smvccl  &
       ,smvcst,smvcwt,sthf,sthu,z0_tile,z0h_tile_bare,soil_clay

  USE prognostics, ONLY :  &
!  imported arrays with intent(inout)
     canht_ft,canopy,canopy_gb,cs                                  &
    ,di,di_ncat,k_sice,gc,gs,lai,rgrain,rho_snow_grnd,smc,smcl     &
    ,snow_mass,snow_mass_sea_ncat,snow_grnd,snow_tile,soot,t_soil  &
    ,ti,tstar_tile,z0msea,nsnow,snowdepth,rgrainl                  &
    ,sice,sliq,tsnow

!  USE route_mod, ONLY :  &
!  imported araays with intent(inout)
!     roffAccumLand

  USE screen, ONLY : q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m


  USE sea_ice, ONLY :  &
!  imported scalars
    beta,dalb_bare_wet,dt_bare,pen_rad_frac            &
   ,sw_alphab,sw_alphac,sw_alpham,sw_dtice

  USE snow_param, ONLY :  &
!  imported arrays with intent(inout)
     ds

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     can_model,l_co2_interactive,l_cosz,l_dust          &
    ,l_pdm,l_phenol,l_snow_albedo                       &
    ,l_top,l_trif_eq,l_triffid                          &
    ,route,routeOnly,l_aggregate,can_rad_mod            &
    ,ilayers,L_Q10,lq_mix_bl,L_CTILE                    &
    ,L_spec_z0,L_SICE_HEATFLUX,L_INLAND                 &
    ,l_sice_multilayers, i_sea_alb_method               &
    ,L_SOIL_SAT_DOWN,l_anthrop_heat_src,buddy_sea       &
!CABLE{
	! set .TRUE. to use CABLE LSM and switch corresponding
	! JULES components on/off [default=.FALSE.]
    ,L_cable 

  !jhan: revise HACKs here
  USE cable_data_mod, ONLY : cable_control, cable_control2, &
                        cable_control3, cable_control4, cable_control5, &
                        cable_control6, cable_control7
!jhan: check this mype is even what we want here
!AND if we can get it another way. ATM disabled
!#ifdef CABLE-UM  
!  USE UM_ParCore, ONLY : mype
!#else
!  USE parallel_mod, ONLY : mype => task_id
  USE soil_param, ONLY : dzsoil 
  USE nstypes, ONLY : npft 
  USE model_grid_mod, ONLY : latitude, longitude 
!#endif  
!CABLE}
  
  USE top_pdm, ONLY :  &
!  imported arrays with intent(in)
     a_fsat,a_fwet,c_fsat,c_fwet,fexp,gamtot,ti_mean,ti_sig  &
!  imported arrays with intent(inout)
    ,fsat,sthzw,zw  &
!  imported arrays with intent(out)
   ,drain,dun_roff,fch4_wetl,fwetl,qbase,qbase_zw,inlandout_atm

  USE trifctl, ONLY :  &
!  imported scalars with intent(in)
     PHENOL_PERIOD,TRIFFID_PERIOD  &
!  imported scalars with intent(inout)
    ,ASTEPS_SINCE_TRIFFID  &
!  imported arrays with intent(inout)
    ,G_LEAF_ACC,NPP_FT_ACC,G_LEAF_PHEN_ACC,RESP_W_FT_ACC,RESP_S_ACC  &
    ,GPP,NPP,RESP_P,G_LEAF,G_LEAF_PHEN,GPP_FT,NPP_FT,RESP_P_FT       &
    ,RESP_S,RESP_W_FT,LAI_PHEN,C_VEG,CV,G_LEAF_DAY,G_LEAF_DR_OUT     &
    ,LIT_C,LIT_C_MN,NPP_DR_OUT,RESP_W_DR_OUT,RESP_S_DR_OUT,FRAC_AGR

  USE u_v_grid, ONLY :  &
!  imported arrays
     dtrdz_charney_grid_1,u_0_p,u_1_p,v_0_p,v_1_p

  USE zenith_mod, ONLY :  &
!  imported procedures
     zenith

  USE ozone_vars, ONLY :  &
     o3

  USE atm_fields_bounds_mod, ONLY: tdims, udims, vdims, pdims,  &
                                   pdims_s, udims_s, vdims_s

  USE surf_param, ONLY : diff_frac

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: a_step     ! IN Atmospheric timestep number.

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER :: error        ! OUT 0 - AOK;
!                         !     1 to 7  - bad grid definition detected
  INTEGER :: phenol_call  ! indicates whether phenology is to be called
  INTEGER :: triffid_call ! indicates whether TRIFFID is to be called
  INTEGER :: nstep_trif   ! Number of atmospheric timesteps between calls to
!                         ! TRIFFID vegetation model
  INTEGER :: i,j,k,l,ll,n      ! Loop counters

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  REAL :: CON_RAIN_LAND(LAND_PTS)     ! Convective rain on land points only (kg/m2/s).
  REAL :: LS_RAIN_LAND(LAND_PTS)      ! Large-scale rain on land points only (kg/m2/s).
  REAL :: CON_SNOW_LAND(LAND_PTS)     ! Convective snow on land points only (kg/m2/s).
  REAL :: LS_SNOW_LAND(LAND_PTS)      ! Large-scale snow on land points only (kg/m2/s).
  REAL :: SURF_HT_FLUX_LD(LAND_PTS)   ! Surface heat flux on land (W/m2).
!                                     ! This is the heat flux into the
!                                     ! uppermost subsurface layer (soil or
!                                     ! snow/soil composite) on land.
  REAL :: LYING_SNOW(LAND_PTS)        ! Gridbox snowmass (kg/m2)
  REAL :: albobs_sc(t_i_length*t_j_length,ntype,2)
                                      ! albedo scaling factors to obs

!-------------------------------------------------------------------------------
!  SCREEN (additional)
!-------------------------------------------------------------------------------
  REAL :: T1_SD(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Standard deviation of turbulent
!                                     ! fluctuations of layer 1 temp;
!                                     ! used in initiating convection.
  REAL :: Q1_SD(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Standard deviation of turbulent
!                                     ! flux of layer 1 humidity;
!                                     ! used in initiating convection.
  REAL :: cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)
                                      ! Ratio of CD's reqd for
!                                     ! calculation of 10 m wind. On
!                                     ! U-grid; comments as per RHOKM.
  REAL :: cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                      ! Ratio of CD's reqd for
!                                     ! calculation of 10 m wind. On
!                                     ! V-grid; comments as per RHOKM.
  REAL :: CHR1P5M(LAND_PTS,NTILES)    ! Ratio of coefffs for
!                                     ! calculation of 1.5m temp for
!                                     ! land tiles.
  REAL :: CHR1P5M_SICE(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                                     ! CHR1P5M for sea and sea-ice
!                                (leads ignored).

!-------------------------------------------------------------------------------
! Local radiation variables.
!-------------------------------------------------------------------------------
  REAL :: PHOTOSYNTH_ACT_RAD(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                                     ! Net downward shortwave radiation
!                                     !  in band 1 (w/m2).

!-------------------------------------------------------------------------------
!  LOCAL SURF
!-------------------------------------------------------------------------------
  REAL :: TILE_FRAC(LAND_PTS,NTILES)  ! Tile fractions including
!                                     ! snow cover in the ice tile.
  REAL :: RAD_SICE(tdims%i_start:tdims%i_end,                       &
                   tdims%j_start:tdims%j_end,nice_use)
                                      ! Surface net shortwave and
!                                     ! downward LW radiation for
!                                     ! sea-ice (W/sq m).
  REAL :: ZH(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Height above surface of top of
!                                     ! boundary layer (metres).
  REAL :: RHO_SNOW(LAND_PTS,NTILES,NSMAX)
!                                     ! Snow layer densities (m)
  REAL :: HCONS(LAND_PTS)             ! Thermal conductivity of top
!                                     ! soil layer, including water and
!                                     ! ice (W/m/K)
  REAL :: RADNET_SICE(tdims%i_start:tdims%i_end,                   &
                      tdims%j_start:tdims%j_end,nice_use)
                                      ! Surface net radiation on
!                                     ! sea-ice (W/m2)
  REAL :: EI_SICE(tdims%i_start:tdims%i_end,                   &
                  tdims%j_start:tdims%j_end,nice_use)
!                                     ! Sea ice sublimation
  REAL :: RHOKM_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
!                                     ! Exchange coefficients for
!                                     ! momentum on P-grid
  REAL :: rhokm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)
                                   ! Exchange coefficients for
!                                  ! momentum (on U-grid, with 1st
!                                  ! and last rows undefined or, at
!                                  ! present, set to "missing data")
  REAL :: rhokm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                   ! Exchange coefficients for
!                                  ! momentum (on V-grid, with 1st
!                                  ! and last rows undefined or, at
!                                  ! present, set to "missing data")
  REAL :: RIB(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Mean bulk Richardson number for
!                                     ! lowest layer.
  REAL :: RIB_TILE(LAND_PTS,NTILES)   ! RIB for land tiles.
  REAL :: FME(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Wind mixing "power" (W/m2).
  REAL :: FB_SURF(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Surface flux buoyancy over
!                                     ! density (m^2/s^3)
  REAL :: U_S(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Surface friction velocity (m/s)
  REAL :: ALPHA1(LAND_PTS,NTILES)     ! Mean gradient of saturated
!                                     ! specific humidity with respect
!                                     ! to temperature between the
!                                     ! bottom model layer and tile
!                                     ! surfaces
  REAL :: ALPHA1_SICE(tdims%i_start:tdims%i_end,                &
                      tdims%j_start:tdims%j_end,nice_use)
                                      ! ALPHA1 for sea-ice.
  REAL :: ASHTF_PRIME(tdims%i_start:tdims%i_end,                &
                      tdims%j_start:tdims%j_end,nice_use)
                                      ! Adjusted SEB coefficient for sea-ice
  REAL :: ASHTF_PRIME_TILE(LAND_PTS,NTILES)
!                                     ! Adjusted SEB coefficient for land tiles
  REAL :: FRACA(LAND_PTS,NTILES)      ! Fraction of surface moisture
!                                     ! flux with only aerodynamic
!                                     ! resistance for snow-free land tiles.
  REAL :: RHOSTAR(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Surface air density
  REAL :: RESFS(LAND_PTS,NTILES)      ! Combined soil, stomatal
!                                     ! and aerodynamic resistance
!                                     ! factor for fraction (1-FRACA)
!                                     ! of snow-free land tiles.
  REAL :: RESFT(LAND_PTS,NTILES)      ! Total resistance factor.
!                                     ! FRACA+(1-FRACA)*RESFS for
!                                     ! snow-free land, 1 for snow.
  REAL :: RHOKH(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Grid-box surface exchange
!                                     ! coefficients
  REAL :: RHOKH_TILE(LAND_PTS,NTILES) ! Surface exchange coefficients
!                                     ! for land tiles
  REAL :: RHOKH_SICE(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Surface exchange coefficients
!                                     ! for sea and sea-ice
  REAL :: DTSTAR_TILE(LAND_PTS,NTILES)! Change in TSTAR over timestep
!                                     ! for land tiles
  REAL :: DTSTAR(tdims%i_start:tdims%i_end,                &
                 tdims%j_start:tdims%j_end,nice_use)
                                      ! Change is TSTAR over timestep
!                                     ! for sea-ice
  REAL :: FLANDG_U(udims%i_start:udims%i_end,udims%j_start:udims%j_end)
                                      ! Land frac (on U-grid, with 1st
!                                     ! and last rows undefined or, at
!                                     ! present, set to "missing data")
  REAL :: FLANDG_V(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                      ! Land frac (on V-grid, with 1st
!                                     ! and last rows undefined or, at
!                                     ! present, set to "missing data")

  REAL :: Z0HSSI(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Roughness length for heat and
!                                     ! moisture over sea (m).
  REAL :: Z0H_TILE(LAND_PTS,NTILES)   ! Tile roughness lengths for heat
!                                     ! and moisture (m).
  REAL :: Z0MSSI(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Roughness length for
!                                     ! momentum over sea (m).
  REAL :: Z0M_TILE(LAND_PTS,NTILES)   ! Tile roughness lengths for
!                                     ! momentum.
  REAL :: VSHR(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Magnitude of surface-to-lowest
!                                     ! atm level wind shear (m per s).
  REAL :: CANHC_TILE(LAND_PTS,NTILES) ! Areal heat capacity of canopy
!                                     ! for land tiles (J/K/m2).
  REAL :: WT_EXT_TILE(LAND_PTS,SM_LEVELS,NTILES)
!                                     ! Fraction of evapotranspiration
!                                     ! which is extracted from each
!                                     ! soil layer by each tile.
  REAL :: FLAKE(LAND_PTS,NTILES)      ! Lake fraction.
  REAL :: CT_CTQ_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Coefficient in T and q
!                                     ! tri-diagonal implicit matrix
  REAL :: CQ_CM_U_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)
                                      ! Coefficient in U tri-diagonal
!                                     ! implicit matrix
  REAL :: CQ_CM_V_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                      ! Coefficient in V tri-diagonal
!                                     ! implicit matrix
  REAL :: DQW_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Level 1 increment to q field
  REAL :: DTL_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Level 1 increment to T field
  REAL :: DU_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end)
!                                     ! Level 1 increment to u wind field
  REAL :: DV_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end)
!                                     ! Level 1 increment to v wind field
  REAL :: TI_GB(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! GBM ice surface temperature (K)
  REAL :: OLR(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      !    TOA - surface upward LW on
!                                     !    last radiation timestep
!                                     !    Corrected TOA outward LW
  REAL :: RHOKH_MIX(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Exchange coeffs for moisture.
  REAL :: BQ_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! A buoyancy parameter (beta q tilde).
  REAL :: BT_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! A buoyancy parameter (beta T tilde).
  REAL :: EMIS_SOIL(LAND_PTS)         ! Emissivity of underlying soil

  REAL :: SEA_ICE_ALBEDO(t_i_length,t_j_length,4)  ! Sea ice albedo
  REAL :: OPEN_SEA_ALBEDO(t_i_length,t_j_length,2) ! Open ocean albedo

!-----------------------------------------------------------------------
! Variables used for semi-implicit, semi-Lagrangian scheme in UM
! Not used standalone
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: NumCycles = 1  !  Number of cycles (iterations) for iterative SISL.
  INTEGER, PARAMETER :: CycleNo = 1    !  Iteration no

!-----------------------------------------------------------------------
! These variables are required for prescribed roughness lengths in
! SCM mode in UM - not used standalone
!-----------------------------------------------------------------------
  REAL :: Z0M_SCM(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Fixed Sea-surface roughness
!                                     ! length for momentum (m).(SCM)
  REAL :: Z0H_SCM(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Fixed Sea-surface roughness
!                                     ! length for heat (m). (SCM)

!-----------------------------------------------------------------------
! These variables are INTENT(OUT) in sf_expl, so just define them here.
!-----------------------------------------------------------------------
  REAL :: RECIP_L_MO_SEA(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                                     ! Reciprocal of the surface
!                                     ! Obukhov  length at sea
!                                     ! points. (m-1).
  REAL :: EPOT_TILE(LAND_PTS,NTILES)  ! Local EPOT for land tiles.
  REAL :: RHOKPM(LAND_PTS,NTILES)     ! Land surface exchange coeff.
!                                     ! (Note used with JULES)
  REAL :: RHOKPM_POT(LAND_PTS,NTILES) ! Potential evaporation
!                                     ! exchange coeff.
!                                     ! (Dummy - not used with JULES)
  REAL :: RHOKPM_SICE(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Sea-ice surface exchange coeff.
!                                     ! (Dummy - not used with JULES)
  REAL :: Z0H_EFF(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Effective grid-box roughness
!                                     ! length for heat, moisture (m)
  REAL :: Z0M_GB(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                      ! Gridbox mean roughness length
!                                     ! for momentum (m).
  REAL :: RESP_S_TOT(DIM_CS2)         ! Total soil respiration (kg C/m2/s).
  REAL :: WT_EXT(LAND_PTS,SM_LEVELS)  ! cumulative fraction of transp'n
  REAL :: RA(LAND_PTS)                ! Aerodynamic resistance (s/m).

  LOGICAL, PARAMETER :: l_ukca = .FALSE.   ! switch for UKCA scheme - NEVER USED!!!
  LOGICAL, PARAMETER :: l_aero_classic = .FALSE. ! switch for CLASSIC aerosol scheme - NEVER USED!!!
  LOGICAL, PARAMETER :: L_FLUX_BC = .FALSE.
                      ! SCM logical for prescribed
                      ! surface flux forcing - why is this in
                      ! the surface scheme?!!!

  LOGICAL :: l_dust_diag
                      ! In the standalone JULES, thise switch essentially
                      ! does the same job as l_dust

!-------------------------------------------------------------------------
! These variables are INTENT(IN) to sf_expl, but not used with the
! current configuration of standalone JULES (initialised to 0 below)
!-------------------------------------------------------------------------
REAL :: z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Height of top of lowest uv-layer
REAL :: z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Height of top of lowest Tq-layer
REAL :: ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                            ! Convective downdraught
!                            ! mass-flux at cloud base

!-----------------------------------------------------------------------
! These variables are INTENT(IN) to sf_impl2 but not used with the
! current configuration of JULES
!-----------------------------------------------------------------------
REAL :: rho1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Density on lowest level
REAL :: f3_at_p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Coriolis parameter
REAL :: uStargbm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! GBM surface friction velocity

!-----------------------------------------------------------------------
! These variables are INTENT(INOUT) to sf_impl2 but not used with the
! current configuration of JULES
!-----------------------------------------------------------------------
REAL :: TScrnDcl_SSI(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                            !    Decoupled screen-level temperature
                            !    over sea or sea-ice
REAL :: TScrnDcl_TILE(land_pts,ntiles)
                            !    Decoupled screen-level temperature
                            !    over land tiles
REAL :: tStbTrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                            !    Time since the transition

!-----------------------------------------------------------------------
! New message passing variables for sf_expl
! These were previously calculated inside sf_expl
!-----------------------------------------------------------------------
REAL ::                                                                       &
  u_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  v_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

! In standalone JULES, usages of flandfac and fseafac are preprocessed out
! Hence, these are dummy arrays, initialised below
REAL ::                                                                       &
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

REAL ::                                                                       &
  rhokm_land(pdims_s%i_start:pdims_s%i_end,                                   &
             pdims_s%j_start:pdims_s%j_end),                                  &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,                                    &
            pdims_s%j_start:pdims_s%j_end),                                   &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

! Versions of the above fields on u and v grids - previously calculated
! inside sf_expl
! In standalone JULES, these grids are the same
REAL ::                                                                       &
  rhokm_u_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),          &
  rhokm_v_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),          &
  rhokm_u_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),           &
  rhokm_v_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

! Other variables previously calculated inside sf_expl
REAL ::                                                                       &
  flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),            &
  flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),            &
  fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
  fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

!-----------------------------------------------------------------------------
! Weights for the implicit solver
! In JULES, these are configured so that the coupling is explicit
!-----------------------------------------------------------------------------
REAL ::                                                                       &
  gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                             ! weights for new BL solver
  gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Additional variables required by the new unconditionally stable implicit
! solver (these are set to ensure explicit coupling standalone)
REAL :: ctctq1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: dqw1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: dtl1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: du_star1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end)
REAL :: dv_star1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end)

! With the configuration of gamma2 and du_1, these are always set to 0
! by im_sf_pt2
REAL :: taux_land_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end)
REAL :: taux_ssi_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end)
REAL :: tauy_land_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
REAL :: tauy_ssi_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)


! Variables required for open sea albedo parameterisation:
REAL :: ws10m(t_i_length,t_j_length) ! 10m wind speed (m s-1)
REAL :: chloro(t_i_length,t_j_length) ! ocean nr surface chlorophyll (kg m-3)
INTEGER, PARAMETER :: n_swbands = 1  !  Nnumber of SW bands
REAL, PARAMETER :: spec_band_bb(1) = -1.0 ! spectral band boundary (bb=-1)

!CABLE{
!jhan: declaring vars which are otherwise NA 

! Soil respiration on tiles (kg C/m2/s). kdcorbin, 10/10
real :: resp_s_tile(land_pts,ntiles)

! local dec in sf_exch - but we need here if filling cable%
real, dimension(land_pts,ntiles) ::                                 &
  cd_tile,                                         &
  ch_tile,                                         &
  recip_l_mo_tile
  
! In UM 7.3 this is recieved in cable and passed back to where
! available in glue_rad. NA in JULES. we need this passed to
! tile_albedo online and declared offline
REAL :: land_alb(t_i_length,t_j_length) 
!CABLE}
           
!-----------------------------------------------------------------------
! Initialise these to zero as they are never used
!-----------------------------------------------------------------------
  z0m_scm(:,:)       = 0.0
  z0h_scm(:,:)       = 0.0
  z1_uv_top(:,:)     = 0.0
  z1_tq_top(:,:)     = 0.0
  ddmfx(:,:)         = 0.0
  rho1(:,:)          = 0.0
  f3_at_p(:,:)       = 0.0
  uStargbm(:,:)      = 0.0
  TScrnDcl_SSI(:,:)  = 0.0
  TScrnDcl_TILE(:,:) = 0.0
  tStbTrans(:,:)     = 0.0
  flandfac(:,:)      = 0.0
  fseafac(:,:)       = 0.0

  l_dust_diag = l_dust

! Initialise the olr diagnostic
  olr(:,:) = 0.0

!CABLE{
  cd_tile = 0.0
  ch_tile = 0.0
  recip_l_mo_tile = 0.0

  ! there are four albedos passed to cable_rad_driver: 
  ! 	  albsoil - PACKS cable var and is used
  !    alb_tile - unpacked from CABLE's rad%reffXX
  !    land_albedo_cable - fn(alb_tile) 
  !    land_alb_cable    - fn(land_albedo_cable) 
  ! so in the first instance we dont need this for offline?
  ! online we will need to pass it back though?
  ! Here canopy_gb == canopy_water in atm_Step(7.3) 
  call cable_control( L_cable, a_step, timestep_len, row_length,         &
         rows, land_pts, ntiles, sm_levels, dim_cs1, dim_cs2,                  & 
         latitude, longitude,                                                  &
         land_index, b, hcon, satcon, sathh, smvcst, smvcwt,                   &
         smvccl, albsoil, lw_down, cosz, ls_rain, ls_snow, pstar,              &
         co2_mmr, sthu, smcl, sthf, gs, canopy_gb, land_alb )

  ! canopy here == canopy_tile, called canopy in atmos_physics2
  ! comment implies snow free only?
  call cable_control2( npft, tile_frac, snow_tile, vshr_land,                  &
         canopy, canht_ft, lai, con_rain, con_snow, NPP, NPP_FT,               &  
         GPP, GPP_FT, RESP_S, RESP_S_TOT, RESP_S_TILE,                         &  
         RESP_P, RESP_P_FT, G_LEAF, Radnet_TILE, Lying_snow,                   &  
         surf_roff, sub_surf_roff, tot_tfall ) 

! in ACCESS tl,qw are passed and recieved as 1st layer tl_1, qw_1
! tl_1, qw_1 are available in JULES forcing module
  call cable_control3( tl_1, qw_1 )  

!  !this will be tricky as there is no surf_down_Sw in JULES
!  !it seems in 8.2 that thisdec locally in glue)rad for the 8A bdylayer scheme
!  ! and maybe set in r2_swrad
! in the UM this has to be passed from glue_rad

! spitter routine in cable_radiation doess this, but simplistically

! CABLE needs diffue, direct, vis, NIR, breakdown of SW radiation
! offline CABLE uses simple spitter FUNCTION. online uses UM's rad
! scheme to compute surf_down_sw. CABLE pulls this upto a higher level
! and then passes down to CABLE
! jhan: HACK. just divide sw_down for now
  !call cable_control4( surf_down_sw )        
  call cable_control4( sw_down )        

  call cable_control5( alb_tile, land_albedo,         &
                  TILE_PTS, TILE_INDEX  )        


  call cable_control6( z1_tq, z1_uv, Fland, dzsoil,                            &  
         LAND_MASK, FTL_TILE, FQW_TILE, TSTAR_TILE,                            &
         U_S, U_S_STD_TILE,                                                    & 
         CD_TILE, CH_TILE, & ! local dec in sf_exch
         FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,                              &
         RECIP_L_MO_TILE, & ! local dec in sf_exch
         EPOT_TILE )


call cable_control7( dtl_1, dqw_1, T_SOIL, FTL_1, FQW_1,                       &
       SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE, EI_TILE,                      &
       T1P5M_TILE, Q1P5M_TILE, MELT_TILE )

!CABLE}

!------------------------------------------------------------------------------
! If we're only doing river routing, most routines need not be called.
!------------------------------------------------------------------------------
  IF ( .NOT. routeOnly ) THEN

!-------------------------------------------------------------------------------
!   Calculate the cosine of the zenith angle
!-------------------------------------------------------------------------------
    IF ( l_cosz ) THEN
      CALL zenith( t_i_length*t_j_length,cosz )
    ELSE
!     Set cosz to a default of 1.0
      cosz(:) = 1.0
    ENDIF

    ! put something in ocean nr. surface chlorophyll - a global mean value:
    chloro(:,:) = 5.0E-7
    ! and, for now, the 10m wind speed
    IF (a_step == 1) THEN
      ! the u/v at 10m not calculated yet, so make a first guess:
      ws10m(:,:) = 2.0
    ELSE
      DO j=1,t_j_length
        DO i=1,t_i_length
          ws10m(i,j) = SQRT(u10m(i,j)**2 + v10m(i,j)**2)
        END DO
      END DO
    END IF
    CALL ftsa (flandg, ice_fract, ice_fract_ncat, tstar, tstar_sice_ncat,      &
           cosz, ws10m, chloro,                                                &
           snow_mass, snow_mass_sea_ncat, di_ncat,                             &
           sw_alpham, sw_alphac, sw_alphab, sw_dtice,                          &
           dt_bare, dalb_bare_wet, pen_rad_frac, beta,                         &
           row_length*rows, n_swbands, n_swbands, nice_use,                    &
           spec_band_bb, spec_band_bb,                                         &
           sea_ice_albedo,open_sea_albedo )

!-------------------------------------------------------------------------------
!   Calculate albedo on land tiles.
!-------------------------------------------------------------------------------
    CALL TILE_ALBEDO (                                                       &
           t_i_length*t_j_length,land_pts,land_index,ntiles,tile_pts         &
          ,tile_index,l_aggregate,l_snow_albedo,albsoil,albobs_sw,albobs_vis &
          ,albobs_nir,cosz,frac,lai,rgrain,snow_tile,soot,tstar_tile,z0_tile &
          ,alb_tile,land_albedo,albobs_sc,can_rad_mod )

!-------------------------------------------------------------------------------
!   Change radiation to be downward components if not using io_rad_type=1
!   NOTE this assumes that the point is either 100% land or 100% sea-ice.
!   1 downward fluxes provided
!   2 net (all wavelength) downward (in downward longwave variable) and downward
!      shortwave fluxes are provided
!   3 net downward fluxes are provided (in downward variables)
!   One day we would probably like to do this in driveUpdate or the likes, but
!   at present we don't
!   have all the required variables/masks there (I think).
!-------------------------------------------------------------------------------

    IF ( io_rad_type == 2) THEN
!-------------------------------------------------------------------------------
!     Convert net total radiation to net longwave. Net total is currently stored
!     in lw_down. Use the average of diffuse albedos in visible and NIR on land.
!-------------------------------------------------------------------------------
      DO j=1,t_j_length
        DO i=1,t_i_length
          IF ( land_mask(i,j) ) THEN
            lw_down(i,j) = lw_down(i,j) - sw_down(i,j) * ( 1.0 - 0.5 *  &
                                ( land_albedo(i,j,2) + land_albedo(i,j,4) ) )
          ELSE
! Assume not using CICE scheme so albedo same for all rad bands
            lw_down(i,j) = lw_down(i,j)                                       &
                         - sw_down(i,j) * ( 1.0 - sea_ice_albedo(i,j,1) )
          ENDIF
        ENDDO
      ENDDO
    ENDIF   !  io_rad_type

    IF ( io_rad_type == 3 ) THEN
!-------------------------------------------------------------------------------
!     Convert shortwave from net to downward.
!     Net flux is currently stored in sw_down.
!     Use the average of diffuse albedos in visible and NIR on land.
!-------------------------------------------------------------------------------
      DO j=1,t_j_length
        DO i=1,t_i_length
          IF ( land_mask(i,j) ) THEN
            sw_down(i,j) = sw_down(i,j) / ( 1.0 - 0.5 *  &
                                ( land_albedo(i,j,2) + land_albedo(i,j,4) ) )
          ELSE
! Assume not using CICE scheme so albedo same for all rad bands
            sw_down(i,j) = sw_down(i,j) / ( 1.0 - sea_ice_albedo(i,j,1) )
          ENDIF
        ENDDO
      ENDDO
    ENDIF   !  io_rad_type

    IF ( io_rad_type == 2 .OR. io_rad_type == 3) THEN
!-------------------------------------------------------------------------------
!     Convert longwave from net to downward. Net longwave is currently stored in lw_down.
!-------------------------------------------------------------------------------
      DO j=1,t_j_length
        DO i=1,t_i_length
          IF ( .NOT. land_mask(I,J) ) THEN
            DO n=1,nice_use
              lw_down(i,j) = lw_down(i,j) + ice_fract_ncat(i,j,n) * sbcon  &
                                          * tstar_sice_ncat(i,j,n)**4.0
            END DO
          ENDIF
        ENDDO
      ENDDO
      DO n=1,ntiles
        DO l=1,land_pts
          j = ( land_index(l)-1 ) / t_i_length + 1
          i = land_index(l) - ( j-1 ) * t_i_length
          lw_down(i,j) = lw_down(i,j) + frac(l,n) * sbcon * tstar_tile(l,n)**4.0
        ENDDO
      ENDDO
    ENDIF   !  io_rad_type

!-----------------------------------------------------------------------
! Now we know sw_down, we can update diff_frac if required
!-----------------------------------------------------------------------
    IF ( use_diff_rad ) THEN
      k = 0
      DO j = 1,t_j_length
        DO i = 1,t_i_length
          k = k + 1
          IF ( sw_down(i,j) > 1.0 ) THEN
            diff_frac(k) = diff_rad(i,j) / sw_down(i,j)
            diff_frac(k) = MIN( 1.0, diff_frac(k) )
          ELSE
            diff_frac(k) = 0.0
          END IF
        END DO
      END DO
    ELSE
      diff_frac(:) = diff_frac_const
    END IF

!-----------------------------------------------------------------------
!   Calculate radiation for sea ice.
!-----------------------------------------------------------------------
    rad_sice(:,:,:) = 0.0
    DO n = 1, nice_use
      DO l = 1, sice_pts_ncat(n)
        ll =  sice_index_ncat(l,n)
        j = ( ssi_index(l) - 1 ) / t_i_length + 1
        i =  ssi_index(l) - (j - 1) * t_i_length
        rad_sice(i,j,n) = ( 1. - alb_sice(ll,n,1) ) * sw_down(i,j)    &
                        + lw_down(i,j)
      ENDDO
    ENDDO

!-----------------------------------------------------------------------
!   Calculate net SW radiation on tiles.
!   Use the average of diffuse albedos in visible and NIR.
!-----------------------------------------------------------------------
    DO n=1,ntiles
      DO l=1,land_pts
        j = ( land_index(l)-1 ) / t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length
        sw_tile(l,n) = ( 1. - 0.5 * ( alb_tile(l,n,2) + alb_tile(l,n,4)  &
                       ) ) * sw_down(i,j)
      ENDDO
    ENDDO

!-----------------------------------------------------------------------
!   Calculate photosynthetically active radiation (PAR).
!-----------------------------------------------------------------------
    DO j=1,t_j_length
      DO i=1,t_i_length
        photosynth_act_rad(i,j) = 0.5 * sw_down(i,j)
      ENDDO
    ENDDO

!-----------------------------------------------------------------------
!   Calculate buoyancy parameters bt and bq. set ct_ctq_1, cq_cm_u_1,
!   cq_cm_v_1, dtl_1, dqw_1, du_1 , dv_1 all to zero (explicit coupling)
!   and boundary-layer depth.
!-----------------------------------------------------------------------
    DO j=1,t_j_length
      DO i=1,t_i_length
        bt_1(i,j) = 1. / tl_1(i,j)
        bq_1(i,j) = c_virtual / (1. + c_virtual*qw_1(i,j))
        zh(i,j) = 1000.
      ENDDO
    ENDDO

!-----------------------------------------------------------------------
!   Generate the anthropogenic heat for surface calculations
!   dummy variables given for yr, hr, min, sec as they are not used
!-----------------------------------------------------------------------
    CALL generate_anthropogenic_heat(                                         &
      0,                                                                      &
      day_of_year(current_time%year, current_time%month, current_time%day),   &
      0, 0, 0, ntiles, land_pts, frac, l_anthrop_heat_src                     &
    )


!-----------------------------------------------------------------------
! Set message passing variables for call to sf_expl
!-----------------------------------------------------------------------
! Copy wind fields into halo versions (note that halos are 0)
    u_0_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = u_0_p(:,:)
    v_0_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = v_0_p(:,:)
    u_1_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = u_1_p(:,:)
    v_1_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end) = v_1_p(:,:)


!-----------------------------------------------------------------------
!   Explicit calculations.
!-----------------------------------------------------------------------
    CALL sf_expl_l_cable (                                                  &
!     IN values defining field dimensions and subset to be processed :
        land_pts,nice_use                                             &
!     IN  parameters for iterative SISL scheme
       ,numcycles,cycleno                                             &
!     IN parameters required from boundary-layer scheme :
       ,bq_1,bt_1,z1_uv,z1_uv_top,z1_tq,z1_tq_top,qw_1,tl_1           &
!     IN soil/vegetation/land surface data :
       ,land_index,ntiles,sm_levels,canopy,catch,catch_snow,hcon      &
       ,ho2r2_orog,fland,flandg                                       &
       ,snow_tile,sil_orog_land,smvccl,smvcst,smvcwt,sthf,sthu        &
       ,z0_tile,z0h_tile_bare                                         &
!     IN sea/sea-ice data :
       ,ice_fract_ncat,k_sice                                         &
!     IN everything not covered so far :
       ,pstar,lw_down,rad_sice,sw_tile,zh,ddmfx                       &
       ,co2_mmr,co2_3d(1:co2_dim_len,1:co2_dim_row)                   &
       ,l_co2_interactive,l_phenol,l_triffid                          &
       ,asteps_since_triffid,cs,frac,canht_ft,photosynth_act_rad,lai  &
       ,lq_mix_bl,t_soil,ti,tstar,tstar_sea,tstar_sice_ncat           &
       ,tstar_tile,z_land,albsoil,cosz                                &
       ,l_aero_classic,l_dust,l_dust_diag,soil_clay,o3                &
! IN idealised and SCM things
       ,l_spec_z0, z0m_scm, z0h_scm                                   &
!     IN variables for message passing
       ,u_1_px, v_1_px, u_0_px, v_0_px                                &
!     IN STASH flags :-
       ,sfme,sq1p5,st1p5,su10,sv10,sz0heff                            &
!     INOUT data :
       ,l_q10,z0msea,gs,g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc&
!     OUT diagnostic not requiring STASH flags :
       ,recip_l_mo_sea,fqw_1,ftl_1,ftl_tile                           &
       ,radnet_sice,rhokm_1,rib,rib_tile                              &
! OUT variables for message passing
       ,flandfac, fseafac,rhokm_land, rhokm_ssi, cdr10m               &
!     OUT diagnostic requiring STASH flags :
       ,fme                                                           &
!     OUT diagnostics required for soil moisture nudging scheme :
       ,wt_ext,ra                                                     &
!     OUT data required for tracer mixing :
       ,rho_aresist,aresist,resist_b                                  &
       ,rho_aresist_tile,aresist_tile,resist_b_tile                   &
!     OUT data required for mineral dust scheme
       ,r_b_dust,cd_std_dust,u_s_std_tile                             &
!     OUT data required elsewhere in UM system
       ,fb_surf,u_s,t1_sd,q1_sd                                       &
!     OUT data required elsewhere in boundary layer or surface code
       ,alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile,fqw_tile      &
       ,epot_tile,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft           &
       ,rhokh,rhokh_tile,rhokh_sice,dtstar_tile,dtstar                &
       ,h_blend_orog,z0hssi,z0h_tile,z0h_eff,z0m_gb,z0mssi,z0m_tile   &
       ,z0m_eff,chr1p5m,chr1p5m_sice,smc,hcons,vshr,vshr_land,vshr_ssi&
       ,gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft                           &
       ,resp_p_ft,resp_s,resp_s_tot,resp_w_ft                         &
       ,gc,canhc_tile,wt_ext_tile,flake                               &
       ,tile_index,tile_pts,tile_frac,fsmc,emis_tile,emis_soil        &
       )


!-----------------------------------------------------------------------------
! Calculate variables required by sf_impl2
! In the UM, these are message passing variables and variables required by
! the implicit solver
!-----------------------------------------------------------------------------

! Popluate the u and v grid versions of variables defined on the p grid
! that are calculated by sf_expl_l
! Since theta, u, v and p grids are the same and halos are 0 in standalone
! model, this is just a straight copy
    rhokm_u_1(:,:) = rhokm_1(:,:)
    rhokm_u_land(:,:) = rhokm_land(:,:)
    rhokm_v_land(:,:) = rhokm_land(:,:)
    rhokm_u_ssi(:,:) = rhokm_ssi(:,:)
    rhokm_v_ssi(:,:) = rhokm_ssi(:,:)
    rhokm_v_1(:,:) = rhokm_1(:,:)
    flandg_u(:,:) = flandg(:,:)
    flandg_v(:,:) = flandg(:,:)

! In standalone JULES, land points are all land and sea points are all sea
! (i.e. no points are part land/part sea)
    flandfac_u(:,:) = 1.0
    flandfac_v(:,:) = 1.0
    fseafac_u(:,:)  = 1.0
    fseafac_v(:,:)  = 1.0

    IF (l_ctile .AND. buddy_sea == on) THEN

      DO j=1,t_j_length
        DO i=1,t_i_length
          taux_land(i,j) = rhokm_u_land(i,j)* u_1(i,j) * flandfac_u(i,j)
          taux_ssi(i,j)  = rhokm_u_ssi(i,j) * ( u_1(i,j) - u_0(i,j) )   &
                                            * fseafac_u(i,j)
          taux_1(i,j) = flandg_u(i,j)*taux_land(i,j)                    &
                        + (1.-flandg_u(i,j))*taux_ssi(i,j)
        END DO
      END DO

      DO j=1,t_j_length
        DO i=1,t_i_length
          tauy_land(i,j) = rhokm_v_land(i,j)* v_1(i,j) * flandfac_v(i,j)
          tauy_ssi(i,j)  = rhokm_v_ssi(i,j) * ( v_1(i,j) - v_0(i,j) )   &
                                            * fseafac_v(i,j)
          tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                    &
                        + (1.-flandg_v(i,j))*tauy_ssi(i,j)
        END DO
      END DO

    ELSE   ! Standard code

      DO j=1,t_j_length
        DO i=1,t_i_length
          taux_land(i,j) = rhokm_u_land(i,j) * u_1(i,j)
          taux_ssi(i,j) = rhokm_u_ssi(i,j) * ( u_1(i,j) - u_0(i,j) )

          taux_1(i,j) = flandg_u(i,j)*taux_land(i,j)                    &
                        + (1.-flandg_u(i,j))*taux_ssi(i,j)
        END DO
      END DO

      DO j=1,t_j_length
        DO i=1,t_i_length
          tauy_land(i,j) = rhokm_v_land(i,j) * v_1(i,j)
          tauy_ssi(i,j) = rhokm_v_ssi(i,j) * ( v_1(i,j) - v_0(i,j) )

          tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                     &
                      + (1.-flandg_v(i,j))*tauy_ssi(i,j)
        END DO
      END DO

    END IF

    cdr10m_u(:,:) = cdr10m(:,:)
    cdr10m_v(:,:) = cdr10m(:,:)


!-----------------------------------------------------------------------
! Implicit calculations.
!
! In the new boundary layer implicit solver in the UM, sf_impl2 is
! called twice - once with l_correct = .FALSE. and once with
! l_correct = .TRUE.
! The variables gamma1 and gamma2 that determine weights for this new
! solver are set so that the new scheme is the same as the old scheme
! (i.e. fully explicit coupling)
!-----------------------------------------------------------------------

! Set values of inputs for the implicit solver so that we get
! explicit coupling
    ct_ctq_1(:,:) = 0.0
    ctctq1(:,:) = 0.0
    cq_cm_u_1(:,:) = 0.0
    cq_cm_v_1(:,:) = 0.0
    dtl_1(:,:) = 0.0
    dtl1_1(:,:) = 0.0
    dqw_1(:,:) = 0.0
    dqw1_1(:,:) = 0.0
    du_1(:,:) = 0.0
    dv_1(:,:) = 0.0
    du_star1(:,:) = 0.0
    dv_star1(:,:) = 0.0
    gamma1(:,:) = gamma(1)

! To get explicit coupling with the new scheme, gamma2 is 0 everywhere
! for the 1st call, and 1 everywhere for the 2nd call
    gamma2(:,:) = 0.0

! Call sf_impl2 with l_correct = .FALSE.
    CALL sf_impl2_cable (                                                     &
!     IN values defining field dimensions and subset to be processed :
         land_pts,land_index,nice,nice_use,ntiles,tile_index,tile_pts,  &
         sm_levels,canhc_tile,canopy,flake,smc,tile_frac,wt_ext_tile,   &
         fland,flandg,lq_mix_bl,l_flux_bc,                              &
!     IN sea/sea-ice data :
         ice_fract,ice_fract_ncat,k_sice,u_0,v_0                        &
!     IN everything not covered so far :
        ,pstar,lw_down,sw_tile                                          &
        ,t_soil,qw_1,tl_1,u_1,v_1,rhokm_u_1,rhokm_v_1,gamma(1)          &
        ,gamma1,gamma2,alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile  &
        ,dtrdz_charney_grid_1,du_1,dv_1                                 &
        ,fraca,resfs,resft,rhokh,rhokh_tile,rhokh_sice,z1_tq            &
        ,z0hssi,z0mssi,z0h_tile,z0m_tile,cdr10m_u,cdr10m_v              &
        ,chr1p5m,chr1p5m_sice,ctctq1                                    &
        ,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1            &
        ,.FALSE.,flandg_u,flandg_v                                      &
        ,emis_tile,ti,tstar_sea,snow_tile                               &
! IN variables used to calculate cooling at the screen level
        ,l_co2_interactive,co2_mmr,co2_3d,rho1,f3_at_p,ustargbm         &
!     IN STASH flags :-
        ,simlt,smlt,slh,sq1p5,st1p5,su10,sv10                           &
!     INOUT data :
        ,epot_tile,fqw_ice,ftl_ice,dtstar_tile,dtstar,rad_sice          &
        ,tstar_sice_ncat,tstar_ssi,tstar_tile,radnet_sice, fqw_tile     &
        ,fqw_1,ftl_1,ftl_tile,olr,taux_land,taux_ssi,tauy_land,tauy_ssi &
        ,tscrndcl_ssi,tscrndcl_tile,tstbtrans                           &
!     OUT diagnostic not requiring STASH flags :
        ,ecan,ei_tile,esoil_tile,sea_ice_htf,surf_ht_flux               &
        ,surf_ht_flux_land,surf_ht_flux_sice,surf_htf_tile              &
!     OUT diagnostic requiring STASH flags :
        ,sice_mlt_htf,snomlt_surf_htf,latent_heat                       &
        ,q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m                    &
!     OUT data required elsewhere in UM system :
        ,tstar,tstar_land,tstar_sice,le_tile,radnet_tile,e_sea,h_sea    &
        ,taux_1,tauy_1,taux_land_star,tauy_land_star,taux_ssi_star      &
        ,tauy_ssi_star,ecan_tile,ei,ei_sice,esoil,ext,snowmelt          &
        ,melt_tile,rhokh_mix,error                                      &
       )

! Adjust the value of gamma2 to ensure explicit coupling
    gamma2(:,:) = 1.0

! Call sf_impl2 again with l_correct = .TRUE.
    CALL sf_impl2_cable (                                                     &
!     IN values defining field dimensions and subset to be processed :
         land_pts,land_index,nice,nice_use,ntiles,tile_index,tile_pts,  &
         sm_levels,canhc_tile,canopy,flake,smc,tile_frac,wt_ext_tile,   &
         fland,flandg,lq_mix_bl,l_flux_bc,                              &
!     IN sea/sea-ice data :
        ice_fract,ice_fract_ncat,k_sice,u_0,v_0                         &
!     IN everything not covered so far :
        ,pstar,lw_down,sw_tile                                          &
        ,t_soil,qw_1,tl_1,u_1,v_1,rhokm_u_1,rhokm_v_1,gamma(1)          &
        ,gamma1,gamma2,alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile  &
        ,dtrdz_charney_grid_1,du_1,dv_1                                 &
        ,fraca,resfs,resft,rhokh,rhokh_tile,rhokh_sice,z1_tq            &
        ,z0hssi,z0mssi,z0h_tile,z0m_tile,cdr10m_u,cdr10m_v              &
        ,chr1p5m,chr1p5m_sice,ctctq1                                    &
        ,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1            &
        ,.TRUE.,flandg_u,flandg_v                                       &
        ,emis_tile,ti,tstar_sea,snow_tile                               &
! IN variables used to calculate cooling at the screen level
        ,l_co2_interactive,co2_mmr,co2_3d,rho1,f3_at_p,ustargbm         &
!     IN STASH flags :-
        ,simlt,smlt,slh,sq1p5,st1p5,su10,sv10                           &
!     INOUT data :
        ,epot_tile,fqw_ice,ftl_ice,dtstar_tile,dtstar,rad_sice          &
        ,tstar_sice_ncat,tstar_ssi,tstar_tile,radnet_sice,fqw_tile      &
        ,fqw_1,ftl_1,ftl_tile,olr,taux_land,taux_ssi,tauy_land,tauy_ssi &
        ,tscrndcl_ssi,tscrndcl_tile,tstbtrans                           &
!     OUT diagnostic not requiring STASH flags :
        ,ecan,ei_tile,esoil_tile,sea_ice_htf,surf_ht_flux               &
        ,surf_ht_flux_land,surf_ht_flux_sice,surf_htf_tile              &
!     OUT diagnostic requiring STASH flags :
        ,sice_mlt_htf,snomlt_surf_htf,latent_heat                       &
        ,q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m                    &
!     OUT data required elsewhere in UM system :
        ,tstar,tstar_land,tstar_sice,le_tile,radnet_tile,e_sea,h_sea    &
        ,taux_1,tauy_1,taux_land_star,tauy_land_star,taux_ssi_star      &
        ,tauy_ssi_star,ecan_tile,ei,ei_sice,esoil,ext,snowmelt          &
        ,melt_tile,rhokh_mix,error                                      &
        )


!-------------------------------------------------------------------------------
!   Compress precip fields to land only for snow and hydrology calculations
!-------------------------------------------------------------------------------
    DO l = 1,land_pts
      j = ( land_index(l) - 1 ) / t_i_length + 1
      i = land_index(l) - (j - 1) * t_i_length
      ls_rain_land(l) = ls_rain(i,j)
      con_rain_land(l) = con_rain(i,j)
      ls_snow_land(l) = ls_snow(i,j)
      con_snow_land(l) = con_snow(i,j)
    ENDDO

!-------------------------------------------------------------------------------
!   Snow processes.
!-------------------------------------------------------------------------------
    CALL snow ( land_pts,REAL(timestep_len),stf_hf_snow_melt,ntiles,tile_pts,  &
                tile_index,catch_snow,con_snow_land,tile_frac,ls_snow_land,    &
                ei_tile,hcap(:,1),hcons,melt_tile,smcl(:,1),sthf(:,1),         &
                surf_htf_tile,t_soil(:,1),tstar_tile,smvcst(:,1),rgrain,       &
                rgrainl,rho_snow_grnd,sice,sliq,snow_grnd,snow_tile,snowdepth, &
                tsnow,nsnow,ds,hf_snow_melt,lying_snow,rho_snow,snomlt_sub_htf,&
                snow_melt,surf_ht_flux_ld )

!-------------------------------------------------------------------------------
!   Reset snowmelt over land points.
!-------------------------------------------------------------------------------
    DO l=1,land_pts
      j = ( land_index(l)-1 ) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      snowmelt(i,j) = snow_melt(l)
    ENDDO

!-------------------------------------------------------------------------------
!   Land hydrology.
!-------------------------------------------------------------------------------
    CALL hydrol (                                                             &
              lice_pts,lice_index,soil_pts,soil_index,nsnow,                  &
              land_pts,sm_levels,b,catch,con_rain_land,                       &
              ecan_tile,ext,hcap,hcon,ls_rain_land,                           &
              satcon,sathh,snowdepth,surf_ht_flux_ld,REAL(timestep_len),      &
              smvcst,smvcwt,canopy,stf_hf_snow_melt,                          &
              stf_sub_surf_roff,smcl,sthf,sthu,t_soil,                        &
              canopy_gb,hf_snow_melt,smc,snow_melt,                           &
              sub_surf_roff,surf_roff,tot_tfall,                              &
              inlandout_atm,l_inland,ntiles,tile_pts,tile_index,              &
              infil_tile,melt_tile,tile_frac,                                 &
              l_top,l_pdm,fexp,gamtot,ti_mean,ti_sig,cs,                      &
              dun_roff,drain,fsat,fwetl,qbase,qbase_zw,                       &
              zw,sthzw,a_fsat,c_fsat,a_fwet,c_fwet,                           &
              fch4_wetl,dim_cs1,l_soil_sat_down,l_triffid)

  ENDIF   !  not routeonly

!-------------------------------------------------------------------------------
! Call runoff routing.
!-------------------------------------------------------------------------------
!  IF ( route ) CALL route_drive( land_pts,sub_surf_roff  &
!                                ,surf_roff,roffaccumland )


  IF ( .NOT. routeonly ) THEN

!-------------------------------------------------------------------------------
!   Copy land points output back to full fields array.
!-------------------------------------------------------------------------------
    DO l = 1,land_pts
      j = ( land_index(l)-1 ) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      snow_mass(i,j) = lying_snow(l)
    ENDDO

    IF ( .NOT. l_sice_multilayers ) THEN
!-----------------------------------------------------------------------
!   Update sea-ice surface layer temperature.
!-----------------------------------------------------------------------
      CALL sice_htf( flandg,simlt,nice,                                  &
                     di_ncat,ice_fract,ice_fract_ncat,surf_ht_flux_sice, &
                     ti,sice_mlt_htf,ti_gb,sea_ice_htf                   &
                   )
    ELSE
      sea_ice_htf(:,:,:) = 0.0
    END IF

!-------------------------------------------------------------------------------
!   Convert sea and sea-ice fluxes to be fraction of grid-box
!   (as required by sea and sea-ice modellers)
!-------------------------------------------------------------------------------
    DO j=1,t_j_length
      DO i=1,t_i_length
        DO n=1,nice
          surf_ht_flux_sice(i,j,n) = ice_fract_ncat(i,j,n) *                  &
                                                      surf_ht_flux_sice(i,j,n)
          sice_mlt_htf(i,j,n) = ice_fract_ncat(i,j,n) * sice_mlt_htf(i,j,n)
          sea_ice_htf(i,j,n) = ice_fract_ncat(i,j,n) * sea_ice_htf(i,j,n)
        ENDDO
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------
!   If leaf phenolgy is activated, check whether the surface model has run an
!   integer number of phenology calling periods
!-------------------------------------------------------------------------------
    phenol_call=1
    triffid_call=1
    IF ( l_phenol ) phenol_call = MOD ( FLOAT(a_step),                        &
                REAL(phenol_period) * REAL(SECS_IN_DAY) / REAL(timestep_len) )

    IF ( l_triffid ) THEN
      nstep_trif = INT( REAL(SECS_IN_DAY) * REAL(triffid_period) / REAL(timestep_len) )
      IF ( asteps_since_triffid == nstep_trif ) triffid_call = 0
    ENDIF

    IF ( triffid_call == 0 ) THEN
!-------------------------------------------------------------------------------
!     Run includes dynamic vegetation
!-------------------------------------------------------------------------------
      CALL veg2( land_pts,land_index,ntiles,can_model                         &
                ,a_step,asteps_since_triffid                                  &
                ,phenol_period,triffid_period,l_phenol,l_triffid,l_trif_eq    &
                ,REAL(timestep_len),frac_agr,satcon                           &
                ,g_leaf_acc,g_leaf_phen_acc,npp_ft_acc                        &
                ,resp_s_acc,resp_w_ft_acc                                     &
                ,cs,frac,lai,soil_clay,canht_ft                               &
                ,catch_snow,catch,infil_tile,z0_tile,z0h_tile_bare            &
                ,c_veg,cv,lit_c,lit_c_mn,g_leaf_day,g_leaf_phen               &
                ,lai_phen,g_leaf_dr_out,npp_dr_out,resp_w_dr_out              &
                ,resp_s_dr_out )

    ELSE

      IF ( phenol_call == 0 )  &
!-------------------------------------------------------------------------------
!       Run includes phenology, but not dynamic vegetation
!       therefore call veg1 rather than veg2
!-------------------------------------------------------------------------------
        CALL veg1( land_pts,ntiles,can_model,a_step,phenol_period,l_phenol  &
                  ,REAL(timestep_len),satcon,g_leaf_acc,frac,lai,canht_ft   &
                  ,catch_snow,catch,infil_tile,z0_tile,z0h_tile_bare        &
                  ,g_leaf_day,g_leaf_phen,lai_phen )

    ENDIF  !  triffid_call

  ENDIF   !  routeonly

  END SUBROUTINE control
#endif
