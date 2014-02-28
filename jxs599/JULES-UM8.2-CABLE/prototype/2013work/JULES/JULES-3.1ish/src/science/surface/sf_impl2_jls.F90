! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SF_IMPL2-----------------------------------------------
!
!  Purpose: Calculate implicit correction to surface fluxes of heat,
!           moisture and momentum to be used by the unconditionally
!           stable and non-oscillatory BL numerical solver.  Also
!           calculates screen level temperature and humidity as well
!           as 10 m winds.
!
!--------------------------------------------------------------------

!    Arguments :-
SUBROUTINE sf_impl2 (                                             &

! IN values defining field dimensions and subset to be processed :
 off_x,off_y,row_length,rows,n_rows,land_pts,                     &

! IN soil/vegetation/land surface data :
 land_index,land_mask,nice,nice_use,                              &
 ntiles,tile_index,tile_pts,sm_levels,                            &
 canhc_tile,canopy,flake,smc,                                     &
 tile_frac,wt_ext_tile,                                           &
 fland,flandg,                                                    &

! IN sea/sea-ice data :
 di,ice_fract,di_ncat,ice_fract_ncat,k_sice,u_0,v_0,              &

! IN everything not covered so far :
 pstar,lw_down,rad_sice,sw_tile,timestep,                         &
 t_soil,qw_1,tl_1,u_1,v_1,rhokm_u_1,rhokm_v_1,GAMMA,              &
 gamma1,gamma2,alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile,   &

 dtrdz_charney_grid_1,du_1,dv_1,                                  &
 fqw_tile,epot_tile,fqw_ice,ftl_ice,                              &
 fraca,resfs,resft,rhokh,rhokh_tile,rhokh_sice,                   &
 rhokpm,rhokpm_pot,rhokpm_sice,                                   &
 dtstar_tile,dtstar,z1,                                           &
 z0hssi,z0mssi,z0h_tile,z0m_tile,cdr10m_u,cdr10m_v,               &
 chr1p5m,chr1p5m_sice,ct_ctq_1,ctctq1,dqw_1,dtl_1,                &
 dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,             &
 l_neg_tstar,l_correct,flandg_u,flandg_v,                         &
 l_sice_heatflux,                                                 &
 emis_tile,emis_soil,                                             &

! IN additional variables used to calculate atmospheric cooling
!    at the screen level
 l_co2_interactive, co2_mmr, co2_3d,                              &
 rho1, f3_at_p, ustargbm,                                         &

! IN STASH flags :-
 simlt,smlt,slh,sq1p5,st1p5,su10,sv10,                            &

! INOUT data :
 ti,ti_gb,tstar,                                                  &
 tstar_land,tstar_sea,tstar_sice,tstar_sice_cat,tstar_ssi,        &
 tstar_tile,snow_tile,                                            &
 le_tile,radnet_sice,radnet_tile,                                 &
 e_sea,fqw_1,ftl_1,ftl_tile,h_sea,olr,taux_1,tauy_1,              &
 taux_land,taux_land_star,taux_ssi,taux_ssi_star,tauy_land,       &
 tauy_land_star,tauy_ssi,tauy_ssi_star,                           &
 TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,                            &

! OUT Diagnostic not requiring STASH flags :
 ecan,ei_tile,esoil_tile,                                         &
 sea_ice_htf,surf_ht_flux,surf_ht_flux_land,surf_ht_flux_sice,    &
 surf_htf_tile,                                                   &

! OUT diagnostic requiring STASH flags :
 sice_mlt_htf,snomlt_surf_htf,latent_heat,                        &
 q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m,                     &

! OUT data required elsewhere in UM system :
 ecan_tile,ei,ei_sice,esoil,ext,snowmelt,melt_tile,rhokh_mix,     &
 error,                                                           &

! LOGICAL LTIMER
 lq_mix_bl,                                                       &
 l_flux_bc,                                                       &
 ltimer,                                                          &
!CABLE: needs these vars
      resp_s_tot      & 
        )

USE atm_fields_bounds_mod
USE theta_field_sizes, ONLY : t_i_length, t_j_length

USE csigma
USE c_lheat
USE c_r_cp
USE c_0_dg_c
USE surf_param, ONLY : ls, ip_scrndecpl2
USE switches  , ONLY : l_aggregate                                &
                      ,l_flake_model, l_tstar_sice_new            &
                      ,l_cable
USE nstypes   , ONLY : lake
USE lake_mod  , ONLY : surf_ht_flux_lake                          &
                      ,lake_h_ice

USE fluxes, ONLY : anthrop_heat, surf_ht_store

USE ancil_info, ONLY: nsmax,ssi_pts,sea_pts,sice_pts,             &
    sice_pts_ncat,ssi_index,sea_index,sice_index,sice_index_ncat, &
    fssi,sea_frac,sice_frac,sice_frac_ncat, &
    !CABLE
    dim_cs1, dim_cs2 
USE switches,    ONLY: iscrntdiag                                 &
                      ,can_model                                  &
                      ,l_snowdep_surf
USE prognostics, ONLY: nsnow                                      &
                      ,snowdepth,                                 &
! CABLE:
                      smcl, gc, gs, canopy_gb
use p_s_parms, ONLY : sthu, sthf, smvcst
USE trifctl, ONLY : g_leaf, npp, gpp, resp_s, resp_p_ft, resp_p, npp_ft, gpp_ft
USE cable_vars_mod, ONLY : conv_rain, conv_snow, sthf_tile, smcl_tile, &
      tsoil_tile, sthu_tile, isnow_flg3l, snow_rho1l, snow_rho3l, snage_tile,&
      snow_cond, snow_mass3l, snow_depth3l, snow_rho3l, snow_tmp3l
USE forcing, ONLY: ls_rain, ls_snow

!end CABLE   
USE jules_mod,   ONLY: snowdep_surf
USE snow_param,  ONLY: rho_snow_const                             &
                      ,cansnowtile
USE c_perma,     ONLY: rho_ice
USE solinc_data, ONLY: sky, l_skyview

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!  Inputs :-
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.

INTEGER                                                           &
  row_length                                                      &
             ! Local number of points on a row
, rows                                                            &
             ! Local number of rows in a theta field
, n_rows                                                          &
             ! Local number of rows in a v field
, off_x                                                           &
             ! Size of small halo in i
, off_y                                                           &
             ! Size of small halo in j.
,land_pts    ! IN No of land points

! (c) Soil/vegetation/land surface parameters (mostly constant).

LOGICAL                                                           &
 land_mask(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN T if land, F elsewhere.

INTEGER                                                           &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                  !    point in ROW_LENGTH,ROWS is the
!                                  !    Ith land point.

INTEGER                                                           &
 sm_levels                                                        &
                             ! IN No. of soil moisture levels
,ntiles                                                           &
                             ! IN No. of land tiles
,tile_index(land_pts,ntiles)                                      &
                             ! IN Index of tile points
,tile_pts(ntiles)                                                 &
                             ! IN Number of tile points
,nice                                                             &
                             ! IN Number of sea ice categories
,nice_use                    ! IN Number of sea ice categories used 
                             !    fully in surface exchange

REAL                                                              &
 canhc_tile(land_pts,ntiles)                                      &
                             ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
,canopy(land_pts,ntiles)                                          &
                             ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
,flake(land_pts,ntiles)                                           &
                             ! IN Lake fraction.
,smc(land_pts)                                                    &
                             ! IN Available soil moisture (kg/m2).
,tile_frac(land_pts,ntiles)                                       &
                             ! IN Tile fractions including
!                                  ! snow cover in the ice tile.
,wt_ext_tile(land_pts,sm_levels,ntiles)                           &
!                                  ! IN Fraction of evapotranspiration
!                                  !    extracted from each soil layer
!                                  !    by each tile.
,fland(land_pts)                                                  &
                             ! IN Land fraction on land pts.
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
!                                  ! IN Land fraction on all pts.
,flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                             ! IN Land fraction on U grid.
,flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)    &
                             ! IN Land fraction on V grid.
,emis_tile(land_pts,ntiles)                                       &
                             ! IN Emissivity for land tiles
,emis_soil(land_pts)         ! IN Emissivity of underlying soil

! (d) Sea/sea-ice data.

REAL                                                              &
 di(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                             ! IN "Equivalent thickness" of
!                            !     sea-ice GBM agregate (m).
,di_ncat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)&
                             ! IN "Equivalent thickness" of
!                            !     sea-ice categories (m).
,ice_fract(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
                             ! IN Fraction of gridbox covered by
!                            !     sea-ice (decimal fraction).
,ice_fract_ncat(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end,nice)                   &
                             ! IN Fraction of gridbox
!                            !  covered by sea-ice on categories.
,k_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice) &
                             ! IN sea ice effective conductivity in
!                             !     sfc layer on categories (W/m2/k)      
,u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end)         &
                             ! IN W'ly component of surface
!                                  !    current (m/s).
,v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                             ! IN S'ly component of surface
!                                  !    current (m/s).

! (f) Atmospheric + any other data not covered so far, incl control.

REAL                                                              &
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       &
                             ! IN Surface pressure (Pascals).
,lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                             ! IN Surface downward LW radiation
!                                  !    (W/m2).
,rad_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)&
                             ! IN Surface net SW and downward LW
!                                  !    radiation for sea-ice (W/sq m).
,sw_tile(land_pts,ntiles)                                         &
                             ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
,timestep                                                         &
                             ! IN Timestep (seconds).
,t_soil(land_pts,sm_levels)                                       &
                             ! IN Soil temperatures (K).
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                             ! IN Total water content
,tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                             ! IN Ice/liquid water temperature
,u_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end) &
                             ! IN W'ly wind component (m/s)
,v_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end) &
                             ! IN S'ly wind component (m/s)
,rhokm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                             ! IN Exchange coefficients for
!                                  !    momentum (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
,rhokm_u_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end)&
!                                  ! IN Exchange coefficients for
!                                  !    land momentum
!                                  !    (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
,rhokm_u_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end) &
                             ! IN Exchange coefficients for
!                                  !    mean sea momentum
!                                  !    (on U-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
,rhokm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
                             ! IN Exchange coefficients for
!                                  !    momentum (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
,rhokm_v_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)&
!                                  ! IN Exchange coefficients for
!                                  !    land momentum
!                                  !    (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")
,rhokm_v_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
!                                  ! IN Exchange coefficients for
!                                  !    mean sea momentum
!                                  !    (on V-grid, with 1st
!                                  !    and last rows undefined or, at
!                                  !    present, set to "missing data")

REAL GAMMA                   ! IN implicit weight in level 1

REAL                                                              &
 gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                             ! weights for new BL solver
,gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL                                                              &
 alpha1(land_pts,ntiles)                                          &
                             ! IN Mean gradient of saturated
!                                  !    specific humidity with respect
!                                  !    to temperature between the
!                                  !    bottom model layer and tile
!                                  !    surfaces
,alpha1_sice(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice_use)                  &
                             ! IN ALPHA1 for sea-ice.
,ashtf_prime(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice_use)                  &
                             ! IN Adjusted SEB coefficient for
                             !    sea-ice
,ashtf_prime_tile(land_pts,ntiles)                                &
                             ! IN Adjusted SEB coefficient for
                             !    land tiles.
,dtrdz_charney_grid_1(pdims%i_start:pdims%i_end,                  &
                      pdims%j_start:pdims%j_end)                  &
!                                  ! IN -g.dt/dp for model layers.
,fqw_tile(land_pts,ntiles)                                        &
                             ! IN Surface FQW for land tiles
,epot_tile(land_pts,ntiles)                                       &
                             ! IN surface tile potential
!                                  !    evaporation
,fqw_ice(tdims%i_start:tdims%i_end,                               &
         tdims%j_start:tdims%j_end,nice_use)                      &
                             ! IN Surface FQW for sea-ice
,ftl_ice(tdims%i_start:tdims%i_end,                               &
         tdims%j_start:tdims%j_end,nice_use)                      &
                             ! IN Surface FTL for sea-ice
,fraca(land_pts,ntiles)                                           &
                             ! IN Fraction of surface moisture
!                                  !    flux with only aerodynamic
!                                  !    resistance for snow-free land
!                                  !    tiles.
,resfs(land_pts,ntiles)                                           &
                             ! IN Combined soil, stomatal
!                                  !    and aerodynamic resistance
!                                  !    factor for fraction (1-FRACA) of
!                                  !    snow-free land tiles.
,resft(land_pts,ntiles)                                           &
                             ! IN Total resistance factor.
!                                  !    FRACA+(1-FRACA)*RESFS for
!                                  !    snow-free land, 1 for snow.
,rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! IN Grid-box surface exchange
!                                  !     coefficients
!                                  !    (not used for JULES)
,rhokh_tile(land_pts,ntiles)                                      &
                             ! IN Surface exchange coefficients
!                                  !    for land tiles
,rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  &
                             ! IN Surface exchange coefficients
!                                  !    for sea and sea-ice
,rhokpm(land_pts,ntiles)                                          &
                             ! IN Land surface exchange coeff.
!                                  !    (not used for JULES)
,rhokpm_pot(land_pts,ntiles)                                      &
                             ! IN Land surface exchange coeff.
!                                  !    for potential evaporation.
!                                  !    (not used for JULES)
,rhokpm_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) &
                             ! IN Sea-ice surface exchange coeff.
!                                  !    (not used for JULES)
,dtstar_tile(land_pts,ntiles)                                     &
                             ! IN Change in TSTAR over timestep
!                                  !     for land tiles
,dtstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)
                             ! IN Change is TSTAR over timestep
!                                  !     for sea-ice

 REAL                                                             &
 z1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                             ! IN Height of lowest level (i.e.
!                                  !    height of middle of lowest
!                                  !    layer).
,z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
,z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! IN Roughness lengths over sea (m)
,z0h_tile(land_pts,ntiles)                                        &
                             ! IN Tile roughness lengths for heat
!                                  !    and moisture (m).
,z0m_tile(land_pts,ntiles)                                        &
                             ! IN Tile roughness lengths for
!                                  !    momentum.
,cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                             ! IN Ratio of CD's reqd for
!                                  !    calculation of 10 m wind. On
!                                  !    U-grid; comments as per RHOKM.
,cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)    &
                             ! IN Ratio of CD's reqd for
!                                  !    calculation of 10 m wind. On
!                                  !    V-grid; comments as per RHOKM.
,chr1p5m(land_pts,ntiles)                                         &
                             ! IN Ratio of coefffs for calculation
!                                  !    of 1.5m temp for land tiles.
,chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)&
!                                  ! IN CHR1P5M for sea and sea-ice
!                                  !    (leads ignored).
,ct_ctq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                             ! IN Coefficient in T and q
!                                  !    tri-diagonal implicit matrix
,cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                             ! IN Coefficient in U tri-diagonal
!                                  !    implicit matrix
,cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
                             ! IN Coefficient in V tri-diagonal
!                                  !    implicit matrix
,dqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! IN Level 1 increment to q field
,dtl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! IN Level 1 increment to T field
,du_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end)&
                             ! IN Level 1 increment to u wind
!                                  !    field
,dv_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end)&
                             ! IN Level 1 increment to v wind
!                                  !    field
,ctctq1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
,dqw1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
,dtl1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
,du_star1(udims_s%i_start:udims_s%i_end,                          &
          udims_s%j_start:udims_s%j_end)                          &
,dv_star1(vdims_s%i_start:vdims_s%i_end,                          &
          vdims_s%j_start:vdims_s%j_end)                          &
,taux_land_star(udims%i_start:udims%i_end,                        &
                udims%j_start:udims%j_end)                        &
,taux_ssi_star(udims%i_start:udims%i_end,                         &
               udims%j_start:udims%j_end)                         &
,tauy_land_star(vdims%i_start:vdims%i_end,                        &
                vdims%j_start:vdims%j_end)                        &
,tauy_ssi_star(vdims%i_start:vdims%i_end,                         &
               vdims%j_start:vdims%j_end)
!                                  ! IN Additional arrays needed by the
!                                  !    uncond stable BL numerical solver
! IN Additional variables for screen-level diagnostics
LOGICAL, INTENT(IN) :: l_co2_interactive
                                ! Flag for interactive 3-D CO2
REAL, INTENT(IN)    :: co2_mmr
                                ! Initial or fixed mass mixing ratio
                                ! of CO2
REAL, INTENT(IN)    :: co2_3d                                         &
                      (tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)
                                ! 3-D field of CO2
REAL, INTENT(IN)    :: rho1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! Density on lowest level
REAL, INTENT(IN)    :: f3_at_p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! Coriolis parameter
REAL, INTENT(IN)    :: uStarGBM(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! GBM surface friction velocity

LOGICAL                                                           &
 ltimer                                                           &
                             ! IN Logical switch for TIMER diags
,l_neg_tstar                                                      &
                             ! IN Switch for -ve TSTAR error check
,l_flux_bc                                                        &
                             ! IN SCM logical for prescribed
                             !    surface flux forcing
,l_correct                                                        &
                             ! flag used by the new BL solver
,l_sice_heatflux             ! IN T: semi-implicit update of TI

LOGICAL                                                           &
 lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code
!  STASH flags :-

LOGICAL                                                           &
 simlt                                                            &
         ! IN Flag for SICE_MLT_HTF (q.v.)
,smlt                                                             &
         ! IN Flag for SNOMLT_SURF_HTF (q.v.)
,slh                                                              &
         ! IN Flag for LATENT_HEAT (q.v.)
,sq1p5                                                            &
         ! IN Flag for Q1P5M (q.v.)
,st1p5                                                            &
         ! IN Flag for T1P5M (q.v.)
,su10                                                             &
         ! IN Flag for U10M (q.v.)
,sv10    ! IN Flag for V10M (q.v.)

!  In/outs :-

REAL                                                              &
 ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)     &
                             ! IN   Sea-ice surface layer
                             !       temperature (K).
,ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! IN GBM ice surface temperature (K)
                             !     (Not used for JULES)
,tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! OUT   GBM surface temperature (K).
,tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  &
                             ! OUT   Land mean sfc temperature (K)
,tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
                             ! IN    Open sea sfc temperature (K).
,tstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  &
                             ! OUT Ice area mean sea ice surface temperature
,tstar_sice_cat(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end,nice_use)               &
                             ! INOUT   Sea-ice sfc temperature (K).
,tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)   &
                             ! INOUT Sea mean sfc temperature (K).
,tstar_tile(land_pts,ntiles)                                      &
                             ! INOUT Surface tile temperatures
,snow_tile(land_pts,ntiles)                                       &
                             ! INOUT Lying snow on tiles (kg/m2)
,le_tile(land_pts,ntiles)                                         &
                             ! INOUT Surface latent heat flux for
!                                  !     land tiles
,radnet_sice(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice_use)                  &
                             ! INOUT Surface net radiation on
!                                  !       sea-ice (W/m2)
,radnet_tile(land_pts,ntiles)                                     &
                             ! INOUT Surface net radiation on
!                                  !       land tiles (W/m2)
,e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! INOUT Evaporation from sea times
!                                  !       leads fraction. Zero over
!                                  !       land. (kg per square metre
!                                  !       per sec).
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! INOUT Moisture flux between layers
!                                  !       (kg per square metre per sec)
!                                  !       FQW(,1) is total water flux
!                                  !       from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! INOUT FTL(,K) contains net
!                                  !       turbulent sensible heat flux
!                                  !       into layer K from below; so
!                                  !       FTL(,1) is the surface
!                                  !       sensible heat, H.(W/m2)
,ftl_tile(land_pts,ntiles)                                        &
                             ! INOUT Surface FTL for land tiles
,h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! INOUT Surface sensible heat flux
!                                  !       over sea times leads fraction
!                                  !       (W/m2)
,olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)         &
                             ! IN    TOA - surface upward LW on
!                                  !       last radiation timestep
!                                  ! OUT   Corrected TOA outward LW
,taux_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end)      &
                             ! OUT   W'ly component of surface
!                                  !       wind stress (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
,taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                             ! INOUT W'ly component of surface
!                                  !       wind stress over land
!                                  !       (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
,taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                             ! INOUT W'ly component of surface
!                                  !       wind stress over mean sea
!                                  !       (N/sq m). (On
!                                  !       UV-grid with first and last
!                                  !       rows undefined or, at
!                                  !       present, set to missing data
,tauy_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)      &
                             ! OUT   S'ly component of surface
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
,tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
                             ! INOUT S'ly component of land sfc
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX
,tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
!                                  !       wind stress (N/sq m).  On
!                                  !       UV-grid; comments as per TAUX

REAL, INTENT(INOUT) :: TScrnDcl_SSI(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                            !    Decoupled screen-level temperature
                            !    over sea or sea-ice
REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_pts,ntiles)
                            !    Decoupled screen-level temperature
                            !    over land tiles
REAL, INTENT(INOUT) :: tStbTrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                            !    Time since the transition to stable
                            !    conditions


!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-

REAL                                                              &
 ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                             ! OUT Gridbox mean evaporation from
!                                  !     canopy/surface store (kg/m2/s).
!                                  !     Zero over sea.
,esoil_tile(land_pts,ntiles)                                      &
                             ! OUT ESOIL for snow-free land tiles
,sea_ice_htf(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end,nice)                      &
                             ! OUT Heat flux through sea-ice
!                                  !     (W/m2, positive downwards).
!                                  !     (Not used for JULES)
,surf_ht_flux(tdims%i_start:tdims%i_end,                          &
              tdims%j_start:tdims%j_end)                          &
!                                  ! OUT Net downward heat flux at
!                                  !     surface over land and sea-ice
!                                  !     fraction of gridbox (W/m2).
,surf_ht_flux_land(tdims%i_start:tdims%i_end,                     &
                   tdims%j_start:tdims%j_end)                     &
!                                  ! OUT Net downward heat flux at
!                                  !     surface over land
!                                  !     fraction of gridbox (W/m2).
,surf_ht_flux_sice(tdims%i_start:tdims%i_end,                     &
                   tdims%j_start:tdims%j_end,nice)                &
!                                  ! OUT Net category downward heat flux at
!                                  !     surface over sea-ice
!                                  !     fraction of gridbox (W/m2).
,surf_htf_tile(land_pts,ntiles)
!                                  ! OUT Net downward surface heat flux
!                                  !     on tiles (W/m2)

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

REAL                                                              &
 sice_mlt_htf(tdims%i_start:tdims%i_end,                          &
              tdims%j_start:tdims%j_end,nice)                     &
!                                  ! OUT Heat flux due to melting of
!                                  !     sea-ice (Watts per sq metre).
,snomlt_surf_htf(tdims%i_start:tdims%i_end,                       &
                 tdims%j_start:tdims%j_end)                       &
!                                  ! OUT Heat flux required for surface
!                                  !     melting of snow (W/m2).
,latent_heat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end) &
                             ! OUT Surface latent heat flux, +ve
!                                  !     upwards (Watts per sq m).
,q1p5m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! OUT Q at 1.5 m (kg water / kg air).
,q1p5m_tile(land_pts,ntiles)                                      &
                             ! OUT Q1P5M over land tiles.
,t1p5m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! OUT T at 1.5 m (K).
,u10m(udims%i_start:udims%i_end,udims%j_start:udims%j_end)        &
                             ! OUT U at 10 m (m per s).
,t1p5m_tile(land_pts,ntiles)                                      &
                             ! OUT T1P5M over land tiles.
,v10m(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                             ! OUT V at 10 m (m per s).

!-2 Genuinely output, needed by other atmospheric routines :-

REAL                                                              &
 ei(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)          &
                             ! OUT Sublimation from lying snow or
!                                  !     sea-ice (kg/m2/s).
,ei_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)     &
                             ! OUT Sublimation from lying snow
!                                  !     (kg/m2/s).
,ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)&
                             ! OUT Sublimation from sea-ice
!                                  !     (kg/m2/s).
,ei_tile(land_pts,ntiles)                                         &
                             ! OUT EI for land tiles.
,ecan_tile(land_pts,ntiles)                                       &
                             ! OUT ECAN for snow-free land tiles
,esoil(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! OUT Surface evapotranspiration
!                                  !     from soil moisture store
!                                  !     (kg/m2/s).
,ext(land_pts,sm_levels)                                          &
                             ! OUT Extraction of water from each
!                                  !     soil layer (kg/m2/s).
,snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)    &
                             ! OUT Snowmelt (kg/m2/s).
,melt_tile(land_pts,ntiles)                                       &
                             ! OUT Snowmelt on land tiles (kg/m2/s
,rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! OUT Exchange coeffs for moisture.
                             !     (Not used for JULES)

INTEGER                                                           &
 error          ! OUT 0 - AOK;
!                     !     1 to 7  - bad grid definition detected;

!  Workspace :-

REAL                                                              &
 elake_tile(land_pts,ntiles)                                      &
                             ! Lake evaporation.
,melt_ice_tile(land_pts,ntiles)                                   &
                             ! Ice melt on FLake lake tile (kg/m2/s)
,lake_ice_mass(land_pts)                                          &
                             ! areal density equivalent to
                             ! lake ice of a given depth (kg/m2)
,qim_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! Implicit value of first model level
!                                  ! humidity
,tim_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)       &
                             ! Implicit value of first model level
!                                  ! temperature
,tstar_rad4(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)  &
                             ! Effective surface radiative
!                                  ! temperature for land and sea-ice
,tstar_tile_old(land_pts,ntiles)                                  &
!                                  ! Tile surface temperatures at
!                                  ! beginning of timestep.
,tstar_sice_cat_old(tdims%i_start:tdims%i_end,                    &
           tdims%j_start:tdims%j_end,nice_use)                    &
                             ! Sea ice surface T at beginning of timestep
,tsurf(land_pts)                                                  &
                             ! Soil or snow surface layer temp
,sice_melt(tdims%i_start:tdims%i_end,                             &
           tdims%j_start:tdims%j_end,nice)                        &
                             !Melt at surface sea-ice category
,non_lake_frac(tdims%i_start:tdims%i_end,                         &
               tdims%j_start:tdims%j_end)                         &
                             ! total tile fraction for surface types
                             ! other than inland water

,dftl_sice_ncat(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end)                        &
                             ! Increment for ftl_ice from sea-ice 
                             ! melt calculated for each category,
                             ! but un-weighted by category fractions
,dfqw_sice_ncat(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end)                        &
                             ! Increment for fqw_ice from sea-ice 
                             ! melt calculated for each category,
                             ! but un-weighted by category fractions
,dei_sice_ncat(tdims%i_start:tdims%i_end,                         &
               tdims%j_start:tdims%j_end)                         &
                             ! Increment for ei_sice from sea-ice 
                             ! melt calculated for each category,
                             ! but un-weighted by category fractions
,ice_fract_cat_use(tdims%i_start:tdims%i_end,                      &
                   tdims%j_start:tdims%j_end,nice_use)   
                             ! Sea ice category fractions
                             ! If nice_use=1, this is the total ice 
                             ! fraction
REAL, ALLOCATABLE :: tstar_ssi_old(:,:)
                                   ! Sea and sea-ice surface temperature
                                   ! at beginning of timestep --
                                   ! Required only for decoupled diagnosis,
                                   ! so allocatable, and local since it is
                                   ! used only on the predictor step

REAL, ALLOCATABLE :: tstar_sic(:,:,:)
                             !Ice category surface temperature
                             ! Only used if nice_use EQ 1

! dummy arrays required for sea and se-ice to create universal
! routines for all surfaces
REAL                                                              &
 array_one(t_i_length*t_j_length)                                 &
                             ! Array of ones
,array_one_e_six(t_i_length*t_j_length)
                             ! Array of 1.0E6

REAL ::                                                           &
 surf_ht_flux_sice_sm(row_length,rows)                            &
                             ! Sea area mean seaice surface heat flux
,ei_sice_sm(row_length,rows)
                             ! Sea area mean sea ice sublimation


!  Local scalars :-

REAL                                                              &
 lat_ht     ! Latent heat of evaporation for snow-free land
!                 ! or sublimation for snow-covered land and ice.

INTEGER                                                           &
 i,j                                                              &
            ! LOCAL Loop counter (horizontal field index).
,k                                                                &
            ! LOCAL Tile pointer
,l                                                                &
            ! LOCAL Land pointer
,n          ! LOCAL Loop counter (tile index).

LOGICAL                                                           &
 l_sice_new_code   ! Controls the sea ice temperature calculation
                   ! (See comments in sea ice section below)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!CABLE: added this var here for implicit call
   REAL ::                                                                     &
      RESP_S_TOT(DIM_CS2)    
     
IF (lhook) CALL dr_hook('SF_IMPL2',zhook_in,zhook_handle)

error = 0

array_one(:)=1.0
array_one_e_six(:)=1.0e6

! Set up sea ice field depending on nice_use
IF (nice_use > 1) THEN
  ! Use all categories fully in surface exchange
  ice_fract_cat_use(:,:,:) = ice_fract_ncat(:,:,:)
ELSE  ! nice_use=1
  ice_fract_cat_use(:,:,1) = ice_fract(:,:)
END IF


! DEPENDS ON: im_sf_pt2
CALL im_sf_pt2 (                                                  &
 land_pts,land_index,ntiles,tile_index,tile_pts                   &
,flandg,tile_frac,snow_tile,nice_use,ice_fract,ice_fract_cat_use  &
,GAMMA,gamma1,gamma2,alpha1,alpha1_sice                           &
,ashtf_prime,ashtf_prime_tile                                     &
,resft,dtstar_tile,dtstar                                         &
,rhokm_u_1,rhokm_v_1,rhokh_tile,rhokh_sice                        &
,ct_ctq_1,ctctq1,dqw_1,dtl_1,dqw1_1,dtl1_1                        &
,cq_cm_u_1,cq_cm_v_1,du_1,dv_1,du_star1,dv_star1                  &
,flandg_u,flandg_v                                                &
,fqw_1,ftl_1                                                      &
,taux_1,taux_land,taux_land_star,taux_ssi,taux_ssi_star,tauy_1    &
,tauy_land,tauy_land_star,tauy_ssi,tauy_ssi_star                  &
,fqw_tile,epot_tile,ftl_tile,fqw_ice,ftl_ice,e_sea,h_sea          &
,l_correct,l_flux_bc,ltimer                                       &
)


!-----------------------------------------------------------------------

! Calculate surface scalar fluxes, temperatures only at the 1st call
! of the subroutine (first stage of the new BL solver) using standard
! MOSES2 physics and equations. These are the final values for this
! timestep and there is no need to repeat the calculation.
!-----------------------------------------------------------------------

IF ( .NOT. l_correct ) THEN
!-----------------------------------------------------------------------
!! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!-----------------------------------------------------------------------

!Cfpp$ Select(CONCUR)
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      ftl_1(i,j) = ftl_1(i,j)*cp
    END DO
  END DO

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      h_sea(i,j) = cp*h_sea(i,j)
    END DO
  END DO

  DO n=1,nice_use
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        ftl_ice(i,j,n) = cp*ftl_ice(i,j,n)
      END DO
    END DO
  END DO


  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      ftl_tile(l,n) = cp*ftl_tile(l,n)
    END DO
  END DO
iF( l_cable ) &
CALL cable_implicit_driver( LS_RAIN, CONV_RAIN, LS_SNOW, CONV_SNOW,       &
                                  DTL_1,DQW_1, T_SOIL, TSOIL_TILE, SMCL,       &
                                  SMCL_TILE, timestep, SMVCST,STHF, STHF_TILE, &
                                  STHU,STHU_TILE, snow_tile, SNOW_RHO1L,       &
                                  ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L,      &
                                  SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,           &
                                  FTL_1, FTL_TILE, FQW_1, FQW_TILE,    &
                                  TSTAR_TILE,        &
                                  SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE,    &
                                  EI_TILE, RADNET_TILE, SNAGE_TILE,   &
                                  CANOPY, GS, T1P5M_TILE, Q1P5M_TILE,     &
                                  CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,        &
                                  DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,   &
                                  RESP_S_TOT, RESP_P, RESP_P_FT,  &
                                  G_LEAF )   



!-----------------------------------------------------------------------
! Land surface calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Optional error check : test for negative top soil layer temperature
!-----------------------------------------------------------------------
  IF (l_neg_tstar) THEN
    DO l=1,land_pts
      IF (t_soil(l,1) < 0) THEN
        error = 1
        WRITE(6,*) '*** ERROR DETECTED BY ROUTINE SF_IMPL2 ***'
        WRITE(6,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '
        WRITE(6,*) 'LAND POINT ',l
      END IF
    END DO
  END IF

!-----------------------------------------------------------------------
!!   Diagnose the land surface temperature
!-----------------------------------------------------------------------

  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      tstar_tile_old(l,n) = tstar_tile(l,n)
IF( .NOT. l_cable ) &
      tstar_tile(l,n) = tstar_tile_old(l,n) + dtstar_tile(l,n)
    END DO
  END DO


!-----------------------------------------------------------------------
!! 7.  Surface evaporation components and updating of surface
!!     temperature (P245, routine SF_EVAP).
!-----------------------------------------------------------------------
IF( .NOT. l_cable ) &
! DEPENDS ON: sf_evap
  CALL sf_evap (                                                  &
    land_pts,ntiles,                                              &
    land_index,tile_index,tile_pts,sm_levels,ltimer,fland,        &
    ashtf_prime_tile,canopy,dtrdz_charney_grid_1,flake,fraca,     &
    snow_tile,resfs,resft,rhokh_tile,tile_frac,smc,wt_ext_tile,   &
    timestep,GAMMA,fqw_1,fqw_tile,ftl_1,ftl_tile,tstar_tile,      &
    ecan,ecan_tile,elake_tile,esoil,esoil_tile,ei_tile,ext        &
    )

!-----------------------------------------------------------------------
!!     Surface melting of sea-ice and snow on land tiles.
!-----------------------------------------------------------------------

  ei_land(:,:)=0.0
  snowmelt(:,:)=0.0

! Lake initialisation
  melt_ice_tile(:,:)=0.0

  DO n=1,ntiles
! DEPENDS ON: sf_melt
    CALL sf_melt (                                                &
      land_pts,land_index,                                        &
      tile_index(:,n),tile_pts(n),ltimer,flandg,                  &
      alpha1(:,n),ashtf_prime_tile(:,n),dtrdz_charney_grid_1,     &
      resft(:,n),rhokh_tile(:,n),tile_frac(:,n),timestep,GAMMA,   &
      ei_tile(:,n),fqw_1,ftl_1,fqw_tile(:,n),ftl_tile(:,n),       &
      tstar_tile(:,n),snow_tile(:,n),snowdep_surf(:,n),           &
      melt_tile(:,n)                                              &
      )

!-----------------------------------------------------------------------
! thermodynamic, flux contribution of melting ice on the FLake lake tile
!-----------------------------------------------------------------------
    IF (     (l_flake_model   ) &
        .AND.(.NOT.l_aggregate) &
        .AND.(n == lake       ) ) THEN

! lake_h_ice is only initialised if FLake is on.
  lake_ice_mass=lake_h_ice * rho_ice

! DEPENDS ON: sf_melt
    CALL sf_melt (                                                &
      land_pts,land_index,                                        &
      tile_index(:,n),tile_pts(n),ltimer,flandg,                  &
      alpha1(:,n),ashtf_prime_tile(:,n),dtrdz_charney_grid_1,     &
      resft(:,n),rhokh_tile(:,n),tile_frac(:,n),timestep,GAMMA,   &
      ei_tile(:,n),fqw_1,ftl_1,fqw_tile(:,n),ftl_tile(:,n),       &
      tstar_tile(:,n),lake_ice_mass,lake_ice_mass/rho_snow_const, &
      melt_ice_tile(:,n)                                          &
        )
    END IF

!-----------------------------------------------------------------------
!  Increment snow by sublimation and melt
!-----------------------------------------------------------------------
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      ei_land(i,j) = ei_land(i,j) + tile_frac(l,n)*ei_tile(l,n)
      snowmelt(i,j) = snowmelt(i,j) +                             &
                      tile_frac(l,n)*melt_tile(l,n)
    END DO

  END DO

  IF (smlt) THEN
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        snomlt_surf_htf(i,j) = lf*snowmelt(i,j)
      END DO
    END DO
  END IF

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
    surf_ht_flux_land(i,j) = 0.
   END DO
  END DO

  IF (     (l_flake_model   ) &
      .AND.(.NOT.l_aggregate) ) THEN
    DO j=tdims%j_start,tdims%j_end
     DO i=tdims%i_start,tdims%i_end
      surf_ht_flux_lake(i,j) = 0.0
! initialise the non-lake fraction to one, not zero,
! in case there should ever be more than one lake tile, see below
      non_lake_frac(    i,j) = 1.0
     END DO
    END DO
  ENDIF

  DO l=1,land_pts
    j=(land_index(l)-1)/t_i_length + 1
    i = land_index(l) - (j-1)*t_i_length
    tstar_land(i,j) = 0.
  END DO

  IF (l_skyview) THEN
    DO n=1,ntiles
      DO k=1,tile_pts(n)
        l = tile_index(k,n)
        j=(land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length
        radnet_tile(l,n) = sw_tile(l,n) +   emis_tile(l,n)*       &
          sky(i,j)*( lw_down(i,j) - sbcon*tstar_tile(l,n)**4 )
      END DO
    END DO
  ELSE
    DO n=1,ntiles
      DO k=1,tile_pts(n)
        l = tile_index(k,n)
        j=(land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length
        radnet_tile(l,n) = sw_tile(l,n) +   emis_tile(l,n)*       &
                   ( lw_down(i,j) - sbcon*tstar_tile(l,n)**4 )
      END DO
    END DO
  END IF

  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      le_tile(l,n) = lc*ecan_tile(l,n) + lc*esoil_tile(l,n) +     &
                     lc*elake_tile(l,n) + ls*ei_tile(l,n)
      surf_ht_store(l,n) = (canhc_tile(l,n)/timestep) *           &
                           (tstar_tile(l,n) - tstar_tile_old(l,n))
      surf_htf_tile(l,n) = radnet_tile(l,n) + anthrop_heat(l,n) - &
                          ftl_tile(l,n) -                         &
                          le_tile(l,n) -                          &
                          lf*(melt_tile(l,n)+melt_ice_tile(l,n))- &
                          surf_ht_store(l,n)
! separate out the lake heat flux for FLake
      IF (     (l_flake_model   ) &
          .AND.(.NOT.l_aggregate) &
          .AND.(n == lake       ) ) THEN
        surf_ht_flux_lake(i,j) = surf_htf_tile(l,n)
        non_lake_frac(    i,j) = non_lake_frac(i,j) - tile_frac(l,n)
      ELSE
        surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j)           &
                          + tile_frac(l,n) * surf_htf_tile(l,n)
      END IF
      tstar_land(i,j) = tstar_land(i,j)                           &
                 + tile_frac(l,n)*tstar_tile(l,n)
    END DO
  END DO

! normalise the non-lake surface heat flux
  IF (     (l_flake_model   ) &
      .AND.(.NOT.l_aggregate) ) THEN
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
! be careful about gridboxes that are all lake
        IF (non_lake_frac(i,j) > EPSILON(0.0)) THEN
          surf_ht_flux_land(i,j) =   surf_ht_flux_land(i,j)       &
                                   / non_lake_frac(i,j)
        ENDIF
      END DO
    END DO
  END IF


!-----------------------------------------------------------------------
! Optional error check : test for negative surface temperature
!-----------------------------------------------------------------------
  IF (l_neg_tstar) THEN
    DO l=1,land_pts
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      IF (tstar_land(i,j) < 0) THEN
        error = 1
        WRITE(6,*) '*** ERROR DETECTED BY ROUTINE SF_IMPL2 ***'
        WRITE(6,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',l
      END IF
    END DO
  END IF


!-----------------------------------------------------------------------
! Sea and sea-ice surface calculations
!-----------------------------------------------------------------------

! Set control logical (l_sice_new_code) to determine how to calculate the 
! sea ice surface temperature.  l_sice_new_code = T is the method that is 
! compatible with using the sea ice categories fully in the radiation and 
! surface exchange code (so nice_use=nice).  l_sice_new_code = F is the 
! old method that only uses the categories in the implicit surface exchange 
! code (so nice_use = 1).  

! NOTE that nice_use=nice=1 can use either scheme and the choice is
! determined by the logical l_tstar_sice_new set in switches.F90.  The 
! difference between the 2 methods is small but causes differences in the 
! results that are larger than bit level, hence the need to control which 
! method is used.

  IF (nice_use == 1 .AND. nice == 1) THEN
    ! If nice_use=nice=1 : Choice method determined by logical 
    ! l_tstar_sice_new
    l_sice_new_code = l_tstar_sice_new

  ELSE IF (nice_use /= nice) THEN
    ! Old calculation must be used
    l_sice_new_code = .FALSE.

  ELSE IF (nice > 1) THEN
    ! New calculation must be used
    l_sice_new_code = .TRUE.
  END IF

  IF (.not.l_sice_new_code) THEN
    ALLOCATE(tstar_sic(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice))
    tstar_sic(:,:,:)= 0.0
  END IF
  surf_ht_flux_sice_sm(:,:)=0.0
  sice_melt(:,:,:)=0.0
  ei_sice(:,:,:)=fqw_ice(:,:,:)
  sice_mlt_htf(:,:,:)=0.0
  sea_ice_htf(:,:,:)=0.0

!-----------------------------------------------------------------------
! Store old surface temperature for sea and sea-ice if using the
! decoupled diagnostic.
!-----------------------------------------------------------------------
  IF (IScrnTDiag == IP_ScrnDecpl2) THEN
    ALLOCATE(tstar_ssi_old(row_length,rows))
    tstar_ssi_old(:,:) = tstar_ssi(:,:)
  ELSE
    ALLOCATE(tstar_ssi_old(1,1))
  END IF

!-----------------------------------------------------------------------
! Diagnose the surface temperature for points with sea-ice
! Note that k_sice = 2.0*thermal conductivity/surface layer thickness
!-----------------------------------------------------------------------

  IF (l_sice_new_code) THEN       

    ! Update tstar_sice_cat using dtstar
    DO n=1,nice_use
      DO j=tdims%j_start,tdims%j_end
        DO i=tdims%i_start,tdims%i_end
          IF ( flandg(i,j) < 1.0 .AND. ice_fract_ncat(i,j,n) > 0 ) THEN
            tstar_sice_cat_old(i,j,n) = tstar_sice_cat(i,j,n)
            tstar_sice_cat(i,j,n) = tstar_sice_cat_old(i,j,n) + dtstar(i,j,n)
          ENDIF
        END DO
      END DO
    END DO

  ELSE  ! Use old code

    ! Update tstar_sic using the surface heat equation
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        IF ( flandg(i,j) < 1.0 .AND. ice_fract(i,j) > 0. ) THEN
          surf_ht_flux_sice_sm(i,j) = radnet_sice(i,j,1) -          &
             4.0*sbcon*(tstar_sice_cat(i,j,1)**3.0)*dtstar(i,j,1) - &
                       ftl_ice(i,j,1) - ls*fqw_ice(i,j,1)
          DO n=1,nice
            IF (ice_fract_ncat(i,j,n) > 0.) THEN
              tstar_sic(i,j,n) = ti(i,j,n) +                        &
                       surf_ht_flux_sice_sm(i,j)/k_sice(i,j,n)
            END IF
          END DO
        END IF
      END DO
    END DO
  END IF

  IF (l_sice_new_code) THEN

    DO n=1, nice_use
! DEPENDS ON: sf_melt
      CALL sf_melt (                                                &
        ssi_pts,ssi_index,                                          &
        sice_index_ncat(:,n),sice_pts_ncat(n),ltimer,fssi,          &
        alpha1_sice,ashtf_prime,dtrdz_charney_grid_1,               &
        array_one,rhokh_sice,sice_frac_ncat(:,n),timestep,GAMMA,    &
        ei_sice(:,:,n),fqw_1,ftl_1,fqw_ice(:,:,n),ftl_ice(:,:,n),   &
        tstar_sice_cat(:,:,n),array_one_e_six,                      &
        array_one_e_six/rho_snow_const,                             &
        sice_melt(:,:,n)                                            &
        )

      DO k=1,sice_pts_ncat(n)
        l = sice_index_ncat(k,n)
        j=(ssi_index(l)-1)/t_i_length + 1
        i = ssi_index(l) - (j-1)*t_i_length
        IF (simlt) sice_mlt_htf(i,j,n) = lf * sice_melt(i,j,n)
      END DO

    END DO

  ELSE ! Use old code

    DO n=1,nice

! Since sea-ice categories are not actually tiled for their surface
! fluxes here, then the increment to ftl_ice, fqw_ice and
! ei_sice are not correctly weighted in sf_melt. Hence need to keep
! the increments and update ftl_ice, fqw_ice and ei_sice with
! weighted contributions below
      dftl_sice_ncat(:,:)=0.0
      dfqw_sice_ncat(:,:)=0.0
      dei_sice_ncat(:,:)=0.0

! DEPENDS ON: sf_melt
      CALL sf_melt (                                                &
        ssi_pts,ssi_index,                                          &
        sice_index_ncat(:,n),sice_pts_ncat(n),ltimer,fssi,          &
        alpha1_sice,ashtf_prime,dtrdz_charney_grid_1,               &
        array_one,rhokh_sice,sice_frac_ncat(:,n),timestep,GAMMA,    &
        dei_sice_ncat,fqw_1,ftl_1,dfqw_sice_ncat,dftl_sice_ncat,    &
        tstar_sic(:,:,n),array_one_e_six,                           &
        array_one_e_six/rho_snow_const,                             &
        sice_melt(:,:,n)                                            &
        )

      DO k=1,sice_pts_ncat(n)
        l = sice_index_ncat(k,n)
        j=(ssi_index(l)-1)/t_i_length + 1
        i = ssi_index(l) - (j-1)*t_i_length
        IF (simlt) sice_mlt_htf(i,j,n) = lf * sice_melt(i,j,n)
! Add weighted increments to ftl_ice, fqw_ice and ei_sice
        ftl_ice(i,j,1)=ftl_ice(i,j,1)                               &
           +( sice_frac_ncat(l,n)/sice_frac(l) )*dftl_sice_ncat(i,j)
        fqw_ice(i,j,1)=fqw_ice(i,j,1)                               &
           +( sice_frac_ncat(l,n)/sice_frac(l) )*dfqw_sice_ncat(i,j)
        ei_sice(i,j,1)=ei_sice(i,j,1)                               &
           +( sice_frac_ncat(l,n)/sice_frac(l) )*dei_sice_ncat(i,j)
      END DO

    END DO

  END IF

!-----------------------------------------------------------------------
!!     Gridbox-mean surface temperature and net surface heat fluxes
!-----------------------------------------------------------------------

  IF (l_sice_new_code) THEN

    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        surf_ht_flux_sice_sm(i,j) = 0.
        DO n=1,nice
          surf_ht_flux_sice(i,j,n)= 0.
        END DO
        IF (flandg(i,j) < 1.0 .AND. ice_fract(i,j) > 0.) THEN
          tstar_sice(i,j)= 0.0
          DO n=1,nice
            IF (ice_fract_ncat(i,j,n) > 0.) THEN
              tstar_sice(i,j) = tstar_sice(i,j) +                   &
                                (ice_fract_ncat(i,j,n)/             &
                               ice_fract(i,j)) * tstar_sice_cat(i,j,n)
            ELSE
              tstar_sice_cat(i,j,n) = tfs   ! copy setting in sice_htf
            END IF
          END DO
          tstar_ssi(i,j)=(1.-ice_fract(i,j))*tstar_sea(i,j) +       &
                          ice_fract(i,j)*tstar_sice(i,j)


          DO n=1,nice
            IF (ice_fract_ncat(i,j,n) > 0.) THEN
              radnet_sice(i,j,n) = rad_sice(i,j,n) -                &
                                    sbcon*tstar_sice_cat(i,j,n)**4

              surf_ht_flux_sice(i,j,n) = radnet_sice(i,j,n) -       &  
                               ftl_ice(i,j,n) - ls*fqw_ice(i,j,n) - &
                               lf*sice_melt(i,j,n)
              surf_ht_flux_sice_sm(i,j) = surf_ht_flux_sice_sm(i,j)+&
                           (ice_fract_ncat(i,j,n)/ice_fract(i,j))*  &
                                surf_ht_flux_sice(i,j,n)
            END IF
          END DO

        ELSE IF (flandg(i,j) < 1.0) THEN  ! non-icy ocean point
          tstar_sice_cat(i,j,:) = tfs  
        END IF
        IF (NICE .EQ. 1) THEN 
          tstar_sice(i,j) = tstar_sice_cat(i,j,1)  ! Ensure these are the
                                                   ! same at all points
        END IF
      END DO
    END DO

  ELSE   ! Use old code

    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        surf_ht_flux_sice_sm(i,j) = 0.
        DO n=1,nice
          surf_ht_flux_sice(i,j,n)=0.
        END DO
        IF (flandg(i,j) < 1.0 .AND. ice_fract(i,j) > 0.) THEN
          tstar_sice(i,j)= 0.0
          DO n=1,nice
            IF (ice_fract_ncat(i,j,n) > 0.) THEN
              tstar_sice(i,j) = tstar_sice(i,j) +                   &
                                (ice_fract_ncat(i,j,n)/             &
                                 ice_fract(i,j)) * tstar_sic(i,j,n)
            END IF
          END DO
          tstar_ssi(i,j)=(1.-ice_fract(i,j))*tstar_sea(i,j) +       &
                          ice_fract(i,j)*tstar_sice(i,j)
          radnet_sice(i,j,1) = rad_sice(i,j,1) -                    &
                                    sbcon*tstar_sice(i,j)**4
          DO n=1,nice
            IF (ice_fract_ncat(i,j,n) > 0.) THEN
              surf_ht_flux_sice(i,j,n) = radnet_sice(i,j,1) -       &
                               4.0*sbcon*tstar_sice(i,j)**3 *       &
                               (tstar_sic(i,j,n)-tstar_sice(i,j)) - &
                               ftl_ice(i,j,1) - ls*fqw_ice(i,j,1) - &
                               lf*sice_melt(i,j,n)
              surf_ht_flux_sice_sm(i,j) = surf_ht_flux_sice_sm(i,j)+&
                           (ice_fract_ncat(i,j,n)/ice_fract(i,j))*  &
                            surf_ht_flux_sice(i,j,n)

            END IF
          END DO

        ELSE IF (flandg(i,j) < 1.0) THEN  ! non-icy ocean point
          tstar_sice(i,j) = tfs  
        END IF

        tstar_sice_cat(i,j,1) = tstar_sice(i,j)  ! Ensure these are the
                                                 ! same at all points
      END DO
    END DO

  END IF

! Convert sea and sea-ice fluxes to be fraction of grid-box
! (as required by sea and sea-ice modellers)
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      h_sea(i,j)=(1.0-ice_fract(i,j))*h_sea(i,j)
      e_sea(i,j)=(1.0-ice_fract(i,j))*e_sea(i,j)
      surf_ht_flux_sice_sm(i,j)=ice_fract(i,j)*surf_ht_flux_sice_sm(i,j)
      ei_sice_sm(i,j) = 0.0
      DO n=1,nice_use
        ei_sice(i,j,n)=ice_fract_cat_use(i,j,n)*ei_sice(i,j,n)
        ei_sice_sm(i,j)= ei_sice_sm(i,j) + ei_sice(i,j,n)
        ftl_ice(i,j,n)=ice_fract_cat_use(i,j,n)*ftl_ice(i,j,n)
        fqw_ice(i,j,n)=ice_fract_cat_use(i,j,n)*fqw_ice(i,j,n)
        rad_sice(i,j,n)=ice_fract_cat_use(i,j,n)*rad_sice(i,j,n)
        radnet_sice(i,j,n)=ice_fract_cat_use(i,j,n)*radnet_sice(i,j,n)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
! GBM diagnostic calculations
!-----------------------------------------------------------------------

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      qim_1(i,j)=qw_1(i,j) + dqw1_1(i,j)-ctctq1(i,j)*fqw_1(i,j)
      tim_1(i,j)=tl_1(i,j) + dtl1_1(i,j)-ctctq1(i,j)*ftl_1(i,j)/cp
   END DO
  END DO

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      tstar(i,j)=flandg(i,j)*tstar_land(i,j)                      &
        +(1.-flandg(i,j))*tstar_ssi(i,j)
      ei(i,j)=flandg(i,j)*ei_land(i,j)                            &
        +(1.-flandg(i,j))*ei_sice_sm(i,j)
      surf_ht_flux(i,j)=flandg(i,j)*surf_ht_flux_land(i,j)        &
        +(1.-flandg(i,j))*surf_ht_flux_sice_sm(i,j)
      rhokh_mix(i,j)=rhokh(i,j)
    END DO
  END DO

! TOA outward LW radiation after boundary layer

  tstar_rad4(:,:)=0.0

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      tstar_rad4(i,j) = tstar_rad4(i,j) + (1.0-flandg(i,j))*      &
                       (1.0-ice_fract(i,j))*tstar_sea(i,j)**4
      DO n=1,nice_use
        tstar_rad4(i,j) = tstar_rad4(i,j) + (1.0-flandg(i,j))*    &
                  ice_fract_cat_use(i,j,n)*tstar_sice_cat(i,j,n)**4 

      END DO
    END DO
  END DO

  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/t_i_length + 1
      i = land_index(l) - (j-1)*t_i_length
      tstar_rad4(i,j) = tstar_rad4(i,j) + flandg(i,j)*            &
                      tile_frac(l,n)*tstar_tile(l,n)**4
    END DO
  END DO

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      olr(i,j) = olr(i,j) + sbcon*tstar_rad4(i,j)
    END DO
  END DO


!-----------------------------------------------------------------------
!!     Specific humidity and temperature at 1.5 metres.
!-----------------------------------------------------------------------
! DEPENDS ON: screen_tq
  CALL screen_tq (                                                &
    land_pts,ntiles,                                              &
    land_index,tile_index,tile_pts,flandg,                        &
    sq1p5,st1p5,chr1p5m,chr1p5m_sice,pstar,qim_1,resft,           &
    tile_frac,tim_1,tstar_ssi,tstar_tile,                         &
    z0hssi,z0h_tile,z0mssi,z0m_tile,z1,                           &
    timestep,tstar_ssi_old,tstar_tile_old,                        &
    l_co2_interactive, co2_mmr, co2_3d,                           &
    f3_at_p, uStarGBM, rho1,                                      &
    TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,                         &
    q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,                            &
    lq_mix_bl                                                     &
    )

! Release space allocated for the transitional diagnostic.
  DEALLOCATE(tstar_ssi_old)
  if (.not.l_sice_new_code) DEALLOCATE(tstar_sic)

!-----------------------------------------------------------------------
!! 9.  Calculate surface latent heat flux.
!-----------------------------------------------------------------------

  IF (slh) THEN
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        latent_heat(i,j) = lc*fqw_1(i,j)                           &
                          + lf*(flandg(i,j)*ei_land(i,j) +         &
                             (1.-flandg(i,j))*ei_sice_sm(i,j))
      END DO
    END DO
  END IF


!-----------------------------------------------------------------------
! Rescale FTL_1 as it should be used to update the botom row of the
! discrete equation handled by the new BL solver at the next (2nd)
! stage of the scheme.
!-----------------------------------------------------------------------
!fpp$ Select(CONCUR)
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      ftl_1(i,j) = ftl_1(i,j)/cp
    END DO
  END DO

ELSE ! L_correct = true: 2nd stage of the scheme

!-----------------------------------------------------------------------
! Rescale to Watts/m^2 as this is the final call to the imp BL solver
! and FTL_1 will be used by stash diagnostics
!-----------------------------------------------------------------------
!fpp$ Select(CONCUR)
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      ftl_1(i,j) = cp*ftl_1(i,j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  U_V will be updated at 2nd stage of the scheme as the equations
!  providing the implicit surface stresses have been modified
!  consistently with the new scheme.
!-----------------------------------------------------------------------
! U component of 10m wind
  IF (su10) THEN
    DO j=udims%j_start,udims%j_end
      DO i=udims%i_start,udims%i_end
        u10m(i,j) = (u_1(i,j) + du_star1(i,j) + (du_1(i,j) -      &
                     cq_cm_u_1(i,j)*taux_1(i,j)) -                &
                     u_0(i,j))*cdr10m_u(i,j) + u_0(i,j)
      END DO
    END DO
  END IF

! V component of 10m wind
  IF (sv10) THEN
    DO j=vdims%j_start,vdims%j_end
      DO i=vdims%i_start,vdims%i_end
        v10m(i,j) = (v_1(i,j) + dv_star1(i,j) + (dv_1(i,j) -      &
                     cq_cm_v_1(i,j)*tauy_1(i,j)) -                &
                     v_0(i,j))*cdr10m_v(i,j) + v_0(i,j)
      END DO
    END DO
  END IF

! Correct surface stress diagnostics

  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
      taux_land(i,j) = taux_land(i,j) + taux_land_star(i,j)
      taux_ssi(i,j)  = taux_ssi(i,j)  + taux_ssi_star(i,j)
    END DO
  END DO

  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
      tauy_land(i,j) = tauy_land(i,j) + tauy_land_star(i,j)
      tauy_ssi(i,j)  = tauy_ssi(i,j)  + tauy_ssi_star(i,j)
    END DO
  END DO

END IF ! IF .NOT. L_correct

IF (lhook) CALL dr_hook('SF_IMPL2',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_impl2
