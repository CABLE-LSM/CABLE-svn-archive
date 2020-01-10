! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE trif_vars_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module holding various shared variables for TRIFFID
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Diagnostics
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE :: wp_fast_in_gb(:)
                      ! C input to fast-turnover wood product pool
                      ! (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wp_med_in_gb(:)
                      ! C input to medium-turnover wood product pool
                      ! (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wp_slow_in_gb(:)
                      ! C input to slow-turnover wood product pool
                      ! (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wp_fast_out_gb(:)
                      ! C output from fast-turnover wood product pool
                      ! (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wp_med_out_gb(:)
                      ! C output from medium-turnover wood product pool
                      ! (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wp_slow_out_gb(:)
                      ! C output from slow-turnover wood product pool
                      ! (kg/m2/360days).

REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_c_orig_pft(:,:)
                      ! Loss of vegetation carbon due to litter,
                      ! landuse change and fire (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_c_ag_pft(:,:)
                      ! Carbon flux from vegetation to wood product pools
                      ! due to land use change (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_orig_pft(:,:)
                      ! Loss of vegetation nitrogen due to litter,
                      ! landuse change and fire (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_ag_pft(:,:)
                      ! Nitrogen removed from system due to landuse change
                      ! this flux is removed from vegetation and not added
                      ! to any other store (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_c_fire_pft(:,:)
                      ! Loss of vegetation carbon due to fire (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_c_nofire_pft(:,:)
                      ! Loss of vegetation carbon due to litter and
                      ! landuse change (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: burnt_carbon_dpm(:)
                      ! Loss of DPM carbon due fire (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_fire_pft(:,:)
                      ! Loss of vegetation nitrogen due to fire
                      ! (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_nofire_pft(:,:)
                      ! Loss of vegetation nitrogen due to litter and
                      ! landuse change (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: burnt_carbon_rpm(:)
                      ! Loss of RPM carbon due fire (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: veg_c_fire_emission_gb(:)
                      ! Gridbox mean carbon flux to the atmosphere from fire
                      ! (kg/(m2 land)/yr).
REAL(KIND=real_jlslsm), ALLOCATABLE :: veg_c_fire_emission_pft(:,:)
                      ! Carbon flux to the atmosphere from fire per PFT
                      ! (kg/(m2 PFT)/yr).
REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_s_to_atmos_gb(:,:)
                      ! Soil-to-atmosphere respiration flux
                      ! [kg m-2 (360 day)-1].
REAL(KIND=real_jlslsm), ALLOCATABLE :: root_abandon_pft(:,:)
                      ! Root carbon moved to soil carbon during
                      ! landuse change (kg/(m2 PFT)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: root_abandon_gb(:)
                      ! Root carbon moved to soil carbon during
                      ! landuse change (kg/(m2 land)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: harvest_pft(:,:)
                      ! Carbon harvested from crops (kg/(m2 PFT)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: harvest_gb(:)
                      ! Gridbox mean carbon harvested from crops
                      ! (kg/(m2 land)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: root_abandon_n_pft(:,:)
                      ! Root nitrogen moved to soil nitrogen during
                      ! landuse change (kg/(m2 PFT)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: root_abandon_n_gb(:)
                      ! Root nitrogen moved to soil nitrogen during
                      ! landuse change (kg/(m2 land)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: harvest_n_pft(:,:)
                      ! Nitrogen harvested from crops (kg/(m2 PFT)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: harvest_n_gb(:)
                      ! Nitrogen harvested from crops: gridbox mean
                      ! (kg/(m2 land)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_fertiliser_pft(:,:)
                      ! Nitrogen available to crop PFTs in addition
                      ! to soil nitrogen (kg/(m2 PFT)/360days)
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_fertiliser_gb(:)
                      ! N available to crop PFTs in addition to
                      ! soil N: gridbox mean (kg/(m2 land)/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_leaf_pft(:,:)
                      ! Leaf N content scaled by LAI, in sf_stom (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_root_pft(:,:)
                      ! Root N content scaled by LAI_BAL, in sf_stom (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_stem_pft(:,:)
                      ! Stem N content scaled by LAI_BAL, in sf_stom (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_leaf_trif_pft(:,:)
                      ! Total Leaf N content (labile + allocated components)
                      ! scaled by lai_bal, in triffid (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_root_trif_pft(:,:)
                      ! Root N content scaled by LAI_BAL, in triffid (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_stem_trif_pft(:,:)
                      ! Stem N content scaled by LAI_BAL, in triffid (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_leaf_alloc_trif_pft(:,:)
                      ! Leaf N content allocated to leaf via phenology
                      ! (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_leaf_labile_trif_pft(:,:)
                      ! Leaf N content in labile leaf pool
                      ! (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_r_pft(:,:)
                      ! Root maintenance respiration (kg C/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_l_pft(:,:)
                      ! Leaf maintenance respiration (kg C/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lai_bal_pft(:,:)
                      ! Balanced lai from sf_stom.
REAL(KIND=real_jlslsm), ALLOCATABLE :: frac_past_gb(:)
                      ! Fraction of pasture.
REAL(KIND=real_jlslsm), ALLOCATABLE :: pc_s_pft(:,:)
                      ! Carbon available for spreading a PFT (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_ag_pft_diag(:,:)
                      ! Nitrogen on tiles lost through landuse
                      ! including harvest (kg/m2/360d).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_luc(:)
                      ! Nitrogen lost through landuse
                      ! including harvest (kg/m2/360d).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_pft_diag(:,:)
                      ! Nitrogen on tiles flux to soil N
                      ! including harvest (kg/m2/360d).
REAL(KIND=real_jlslsm), ALLOCATABLE :: fapar_diag_pft(:,:)
                      ! Fraction of Absorbed Photosynthetically Active
                      ! Radiation diagnostic.
REAL(KIND=real_jlslsm), ALLOCATABLE :: apar_diag_pft(:,:)
                      ! Absorbed Photosynthetically Active
                      ! Radiation diagnostic in W m-2.
REAL(KIND=real_jlslsm), ALLOCATABLE :: apar_diag_gb(:)
                      ! Absorbed Photosynthetically Active
                      ! Radiation diagnostic in W m-2.
REAL(KIND=real_jlslsm), ALLOCATABLE :: fao_et0(:)
                      ! FAO Penman-Monteith evapotranspiration for
                      ! reference crop (kg m-2 s-1).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_carbon_veg2_gb(:)
                      ! Diagnostic of error in land carbon
                      ! conservation in the veg2 routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_carbon_triffid_gb(:)
                      ! Diagnostic of error in land carbon
                      ! conservation in the triffid routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_veg_triffid_gb(:)
                      ! Diagnostic of error in vegetation carbon
                      ! conservation in the triffid routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_soil_triffid_gb(:)
                      ! Diagnostic of error in soil carbon
                      ! conservation in the triffid routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_prod_triffid_gb(:)
                      ! Diagnostic of error in wood product carbon
                      ! conservation in the triffid routine (kg m-2).

!-----------------------------------------------------------------------------
! Variables added for nitrogen conservation checks
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_nitrogen_triffid_gb(:)
                      ! Diagnostic of error in land nitrogen
                      ! conservation in the triffid routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_vegN_triffid_gb(:)
                      ! Diagnostic of error in vegetation nitrogen
                      ! conservation in the triffid routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_soilN_triffid_gb(:)
                      ! Diagnostic of error in soil nitrogen
                      ! conservation in the triffid routine (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: cnsrv_N_inorg_triffid_gb(:)
                      ! Diagnostic of error in inorganic nitrogen
                      ! conservation in the triffid routine (kg m-2).

!-----------------------------------------------------------------------------
! Variables added for ticket #7,#127 (nitrogen scheme)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE :: deposition_n_gb(:)
                      ! Nitrogen deposition (kg/m2/s).

REAL(KIND=real_jlslsm), ALLOCATABLE :: leafC_pft(:,:)
                      ! Leaf carbon on PFTs (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: rootC_pft(:,:)
                      ! Root carbon on PFTs (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: stemC_pft(:,:)
                      ! Stem carbon on PFTs (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: woodC_pft(:,:)
                      ! Wood carbon on PFTs (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: leafC_gbm(:)
                      ! Leaf carbon (GBM) (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: rootC_gbm(:)
                      ! Root carbon (GBM) (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: woodC_gbm(:)
                      ! Wood carbon (GBM) (kg/m2).

REAL(KIND=real_jlslsm), ALLOCATABLE :: droot_pft(:,:)
                      ! Increment in leaf carbon on PFTs
                      ! (kg m-2 per TRIFFID timestep).
REAL(KIND=real_jlslsm), ALLOCATABLE :: dleaf_pft(:,:)
                      ! Increment in leaf carbon on PFTs
                      ! (kg m-2 per TRIFFID timestep).
REAL(KIND=real_jlslsm), ALLOCATABLE :: dwood_pft(:,:)
                      ! Increment in wood carbon on PFTs
                      ! (kg m-2 per TRIFFID timestep).

REAL(KIND=real_jlslsm), ALLOCATABLE :: root_litC_pft(:,:)
                      ! Root litter C turnover on PFTs (kg/m2/360 day).
REAL(KIND=real_jlslsm), ALLOCATABLE :: leaf_litC_pft(:,:)
                      ! Leaf litter C turnover on PFTs (kg/m2/360 day).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wood_litC_pft(:,:)
                      ! Wood litter C turnover on PFTs (kg/m2/360 day).

REAL(KIND=real_jlslsm), ALLOCATABLE :: root_litN_pft(:,:)
                      ! Root litter N turnover on PFTs (kg/m2/360 day).
REAL(KIND=real_jlslsm), ALLOCATABLE :: leaf_litN_pft(:,:)
                      ! Leaf litter N turnover on PFTs (kg/m2/360 day).
REAL(KIND=real_jlslsm), ALLOCATABLE :: wood_litN_pft(:,:)
                      ! Wood litter N turnover on PFTs (kg/m2/360 day).

REAL(KIND=real_jlslsm), ALLOCATABLE :: litterC_pft(:,:)
                      ! Carbon in local litter production (kg/m2/360 day).
REAL(KIND=real_jlslsm), ALLOCATABLE :: litterN_pft(:,:)
                      ! Nitrogen in local litter production (kgN/m2/360 day).

REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_pft(:,:)
                      ! Nitrogen in total litter production (kgN/m2/360 day).

REAL(KIND=real_jlslsm), ALLOCATABLE :: n_uptake_growth_pft(:,:)
                      ! Vegetation N uptake for growth on PFTs
                      ! (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_demand_growth_pft(:,:)
                      ! Vegetation N demand for growth on PFTs
                      ! (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_demand_lit_pft(:,:)
                      ! Vegetation N demand for balanced litter
                      ! production on PFTs (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_demand_spread_pft(:,:)
                      ! Vegetation N demand for spreading on
                      ! PFTs(kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_uptake_spread_pft(:,:)
                      ! Vegetation N uptake for spreading in PFTs
                      ! (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_uptake_pft(:,:)
                      ! Vegetation N uptake on PFTs (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_demand_pft(:,:)
                      ! Vegetation N demand on PFTs(kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_uptake_gb(:)
                      ! Vegetation N uptake (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_demand_gb(:)
                      ! Vegetation N demand (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_veg_pft(:,:)
                      ! Veg N on PFTs (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_veg_gb(:)
                      ! Veg N (kg/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_loss_gb(:)
                      ! Other N Gaseous Loss (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_fix_pft(:,:)
                      ! Fixed N on PFTs (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_fix_gb(:)
                      ! Fixed N (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_leach_soilt(:,:)
                      ! Leached N (kg/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_gas_gb(:,:)
                      ! Mineralised N Gas Emmissions (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: dpm_ratio_gb(:)
                      ! Ratio of Decomposatble Plant Material to Resistant.
REAL(KIND=real_jlslsm), ALLOCATABLE :: dnveg_pft(:,:)
                      ! Increment in veg N on PFTs (kg m-2 per TRIFFID
                      ! timestep).
REAL(KIND=real_jlslsm), ALLOCATABLE :: dcveg_pft(:,:)
                      ! Increment in veg C on PFTs (kg m-2 per TRIFFID
                      ! timestep).
REAL(KIND=real_jlslsm), ALLOCATABLE :: dnveg_gb(:)
                      ! Increment in veg N (kg m-2 per TRIFFID timestep).
REAL(KIND=real_jlslsm), ALLOCATABLE :: dcveg_gb(:)
                      ! Increment in veg C (kg m-2 per TRIFFID timestep).

REAL(KIND=real_jlslsm), ALLOCATABLE :: immob_n_gb(:,:,:)
                      ! Immobilised N on soil pools (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: immob_n_pot_gb(:,:,:)
                      ! Unlimited immobilised N on soil pools
                      ! (kg/m2/360 days).

REAL(KIND=real_jlslsm), ALLOCATABLE :: minl_n_gb(:,:,:)
                      ! Mineralised N on soil pools (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: minl_n_pot_gb(:,:,:)
                      ! Unlimited mineralised N on soil pools
                      ! (kg/m2/360 days).

REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_s_diag_gb(:,:,:)
                      ! Carbon in diagnosed soil respiration after TRIFFID
                      ! on soil pools (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_s_pot_diag_gb(:,:,:)
                      ! Carbon in diagnosed unlimited soil respiration after
                      ! TRIFFID (kg/m2/360 days).

REAL(KIND=real_jlslsm), ALLOCATABLE :: fn_gb(:,:)
                      ! Nitrogen decomposition factor.

REAL(KIND=real_jlslsm), ALLOCATABLE :: lit_n_t_gb(:)
                      ! Total N litter flux (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: exudates_pft(:,:)
                      ! Unallocated C due to lack of N, added to soil
                      ! respiration on PFTs ((kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: exudates_gb(:)
                      ! Unallocated C due to lack of N, added to soil
                      ! respiration ((kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: npp_n(:,:)
                      ! NPP post N limitation on PFTs (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: npp_n_gb(:)
                      ! NPP post N limitation (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: gpp_pft_out(:,:)
                      ! Gross Primary Productivity on PFTs (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: gpp_gb_out(:)
                      ! Gross Primary Productivity (GBM) (kg/m2/360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: gpp_pft_acc(:,:)
                      ! Accumulated GPP on PFTs for calculating plant
                      ! respiration after N limitation on TRIFFID timesteps,
                      ! within TRIFFID itself (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: gpp_gb_acc(:)
                      ! Accumulated GPP (GBM) for calculating plant
                      ! respiration after N limitation on TRIFFID timesteps,
                      ! within TRIFFID itself (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_p_actual_pft(:,:)
                      ! Carbon in plant respiration on PFTs after NPP
                      ! reduction due to nitrogen limitation (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: resp_p_actual_gb(:)
                      ! Carbon in plant respiration (GBM) after NPP
                      ! reduction due to nitrogen limitation (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_leach_gb_acc(:)
                      ! Accumulated leached nitrogen term for outputting
                      ! on leached N on TRIFFID timesteps via
                      ! diagnostics_veg.F90 (kg/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: g_burn_pft_acc(:,:)
                      ! Burnt area disturbance accumulation variable
                      ! If this was in the UM it would be considered a
                      ! prognostic (m2/m2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: g_burn_pft(:,:)
                      ! Burnt area disturbance (m2/m2/360days).
REAL(KIND=real_jlslsm), ALLOCATABLE :: g_burn_gb(:)
                      ! Burnt area disturbance per gb (m2/m2/360days).

!------------------------------------------------------------------------------
! Variables to pass fluxes to a soil model.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_fertiliser_add(:,:,:)
                      ! Nitrogen added to each soil layer by fertiliser
                      ! (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_fix_add(:,:,:)
                      ! Nitrogen added to each soil layer by fixation
                      ! (kg m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE :: n_uptake_extract(:,:,:)
                      ! Nitrogen removed from each soil layer by plant
                      ! uptake (kg m-2).


END MODULE trif_vars_mod
