! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holds surface parameters for each Plant Functional Type (but
! not parameters that are only used by TRIFFID).



! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE pftparm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Radiation and albedo parameters.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                       &
 orient(:)                                                                    &
                 ! Flag for leaf orientation: 1 for horizontal,
!                    0 for spherical.
  ,fsmc_mod(:)
                   ! Flag for whether water stress is calculated from 
                   ! available water in layers weighted by root fraction (0) 
                   ! or
                   ! whether water stress calculated from available 
                   ! water in root zone (1)

REAL, ALLOCATABLE ::                                                          &
 albsnc_max(:)                                                                &
                 ! Snow-covered albedo for large LAI.
,albsnc_min(:)                                                                &
                 ! Snow-covered albedo for zero LAI.
,albsnf_maxu(:)                                                               &
                 ! Max Snow-free albedo (max LAI) when scaled to obs
,albsnf_max(:)                                                                &
                 ! Snow-free albedo for large LAI.
,albsnf_maxl(:)                                                               &
                 ! Min Snow-free albedo (max LAI) when scaled to obs
,alniru(:)                                                                    &
                 ! upper limit on alnir, when scaled to albedo obs
,alnir(:)                                                                     &
                 ! Leaf reflection coefficient for near infra-red.
,alnirl(:)                                                                    &
                 ! lower limit on alnir, when scaled to albedo obs
,alparu(:)                                                                    &
                 ! upper limit on alpar, when scaled to albedo obs
,alpar(:)                                                                     &
                 ! Leaf reflection coefficient for PAR.
,alparl(:)                                                                    &
                 ! lower limit on alpar, when scaled to albedo obs
,kext(:)                                                                      &
                 ! Light extinction coefficient - used to
!                    calculate weightings for soil and veg.
  ,kpar(:)                                                                    &
                   ! PAR Extinction coefficient
!                    (m2 leaf/m2 ground)
  ,lai_alb_lim(:)                                                             &
!                  ! Lower limit on permitted LAI in albedo
  ,omegau(:)                                                                  &
                   ! upper limit on omega, when scaled to albedo obs
  ,omega(:)                                                                   &
                   ! Leaf scattering coefficient for PAR.
  ,omegal(:)                                                                  &
                   ! lower limit on omega, when scaled to albedo obs
  ,omniru(:)                                                                  &
                   ! upper limit on omnir, when scaled to albedo obs
  ,omnir(:)                                                                   &
                   ! Leaf scattering coefficient for near infra-red.
  ,omnirl(:)
                   ! lower limit on omnir, when scaled to albedo obs

!-----------------------------------------------------------------------------
! Parameters for phoyosynthesis and respiration.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                       &
 c3(:)           ! Flag for C3 types: 1 for C3 Plants,
!                    0 for C4 Plants.

REAL, ALLOCATABLE ::                                                          &
 alpha(:)                                                                     &
                 ! Quantum efficiency of photosynthesis
!                    (mol CO2/mol PAR photons).
  ,dqcrit(:)                                                                  &
                   ! Critical humidity deficit (kg H2O/kg air), used with the
                   ! Jacobs closure.
  ,fd(:)                                                                      &
                   ! Dark respiration coefficient.
  ,f0(:)                                                                      &
                   ! CI/CA for DQ = 0, used with the Jacobs closure.
  ,g1_stomata(:)                                                              &
                   ! Parameter g1 for the Medlyn et al. (2011) model of
                   ! stomatal conductance (kPa**0.5) - see Eqn.11 of
                   ! doi: 10.1111/j.1365-2486.2012.02790.x.
  ,neff(:)                                                                    &
                  ! Constant relating VCMAX and leaf N (mol/m2/s)
!                   from Schulze et al. 1994
!                   (AMAX = 0.4e-3 * NL  - assuming dry matter is
!                   40% carbon by mass)
!                   and Jacobs 1994:
!                   C3 : VCMAX = 2 * AMAX ;
!                   C4 : VCMAX = AMAX  ..
  ,nl0(:)                                                                     &
                   ! Top leaf nitrogen concentration
!                    (kg N/kg C).
  ,nr_nl(:)                                                                   &
                   ! Ratio of root nitrogen concentration to
!                    leaf nitrogen concentration.
  ,ns_nl(:)                                                                   &
                   ! Ratio of stem nitrogen concentration to
!                    leaf nitrogen concentration.
  ,nr(:)                                                                      &
                   ! Root nitrogen concentration (kg N/kg C)
  ,nsw(:)                                                                     &
                   ! Stem nitrogen concentration (kg N/kg C)
  ,hw_sw(:)                                                                   &
                   ! Heart:Stemwood Ratio (kg N/kg N)
  ,can_struct_a(:)                                                            &
                   ! Pinty canopy structure factor corresponding to overhead
                   ! sun (dimensionless)
  ,r_grow(:)                                                                  &
                   ! Growth respiration fraction.
  ,tlow(:)                                                                    &
                   ! Lower temperature for photosynthesis (deg C).
  ,tupp(:)                                                                    &
                   ! Upper temperature for photosynthesis (deg C).
  ,lma(:)                                                                     &
                   ! Leaf mass by area (1/SLA), (kg leaf/ m2)
  ,nmass(:)                                                                   &
                   ! Leaf nitrogen dry weight (g N/g leaf)
  ,vsl(:)                                                                     &
                   ! Slope of the Narea to Vcmax relationship
                   ! from Kattge et al. (2009)
  ,vint(:)                                                                    &
                   ! Y intercept of the Narea to Vcmax relationship
                   ! from Kattge et al. (2009)
  ,kn(:)                                                                      &
                   ! Exponential for N profile in canopy, used with
                   ! can_rad_mod=4, 5 (decay is a function of layers).
  ,knl(:)                                                                     &
                   ! Decay coefficient for N profile in canopy, used with
                   ! can_rad_mod=6 (decay is a function of LAI).
  ,q10_leaf(:)                                                                &
                   ! Factor for leaf respiration.
  !---------------------------------------------------------------------------
  ! Parameters for the Farquhar photosynthesis model.
  ! Jmax is the potential rate of electron transport.
  ! Vcmax is the maximum rate of carboxylation of Rubisco.
  !---------------------------------------------------------------------------
  ,act_jmax(:)                                                                &
                   ! Activation energy for temperature response of Jmax
                   ! (J mol-1).
  ,act_vcmax(:)                                                               &
                   ! Activation energy for temperature response of Vcmax
                   ! (J mol-1).
  ,alpha_elec(:)                                                              &
                   ! Quantum yield of electron transport
                   ! (mol electrons/mol PAR photons).
  ,deact_jmax(:)                                                              &
                   ! Deactivation energy for temperature response of Jmax
                   ! (J mol-1). This describes the rate of decrease
                   ! above the optimum temperature.
  ,deact_vcmax(:)                                                             &
                   ! Deactivation energy for temperature response of Vcmax
                   ! and Jmax (J mol-1). This describes the rate of decrease
                   ! above the optimum temperature.
  ,ds_jmax(:)                                                                 &
                   ! Entropy factor for temperature reponse of Jmax 
                   ! (J mol-1 K-1).
  ,ds_vcmax(:)                                                                &
                   ! Entropy factor for temperature reponse of Vcmax 
                   ! (J mol-1 K-1).
  ,jv25_ratio(:)
                   ! Ratio of Jmax to Vcmax at 25 deg C
                   ! (mol electrons mol-1 CO2).

!-----------------------------------------------------------------------------
! Allometric and other parameters.
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 a_wl(:)                                                                      &
                 ! Allometric coefficient relating the target
!                    woody biomass to the leaf area index
!                    (kg C/m2)
  ,a_ws(:)                                                                    &
                   ! Woody biomass as a multiple of live
!                    stem biomass.
  ,b_wl(:)                                                                    &
                   ! Allometric exponent relating the target
!                    woody biomass to the leaf area index.
  ,eta_sl(:)                                                                  &
                   ! Live stemwood coefficient (kg C m-1 (m2 leaf)-1)
  ,sigl(:)
                   ! Specific density of leaf carbon
!                    (kg C/m2 leaf).

!-----------------------------------------------------------------------------
! Phenology parameters.
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 g_leaf_0(:)                                                                  &
                 ! Minimum turnover rate for leaves (/360days).

,dgl_dm(:)                                                                    &
                 ! Rate of change of leaf turnover rate with
!                    moisture availability.
  ,fsmc_of(:)                                                                 &
                   ! Moisture availability below which leaves
!                    are dropped.
  ,dgl_dt(:)                                                                  &
                   ! Rate of change of leaf turnover rate with
!                    temperature (/K)
  ,tleaf_of(:)
                   ! Temperature below which leaves are
!                    dropped (K)

!-----------------------------------------------------------------------------
! Parameters for hydrological, thermal and other "physical" characteristics.
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 catch0(:)                                                                    &
                 ! Minimum canopy capacity (kg/m2).
,dcatch_dlai(:)                                                               &
                 ! Rate of change of canopy capacity with LAI.
,infil_f(:)                                                                   &
                 ! Infiltration enhancement factor.
,glmin(:)                                                                     &
                 ! Minimum leaf conductance for H2O (m/s).
,dz0v_dh(:)                                                                   &
                 ! Rate of change of vegetation roughness
!                    length with height.
,z0v(:)                                                                       &
                 ! Specified vegetation roughness length.
!                    used with l_spec_veg_z0 = .true.
  ,rootd_ft(:)                                                                &
                   ! e-folding depth (m) of the root density.
  ,psi_close(:)                                                               &
                   ! soil matric potential (Pa) below which soil moisture 
                   ! stress factor fsmc is zero. Should be negative.
  ,psi_open(:)                                                                &
                   ! soil matric potential (Pa) above which soil moisture 
                   ! stress factor fsmc is one. Should be negative.
  ,fsmc_p0(:)                                                                 &
                   ! parameter in calculation of the
                   ! soil moisture at which the plant begins to experience
                   ! water stress
  ,emis_pft(:)                                                                &
                   !  Surface emissivity
  ,dust_veg_scj(:) ! Dust emission scaling factor for  each PFT



!-----------------------------------------------------------------------------
! Parameters for ozone damage
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 fl_o3_ct(:)                                                                  &
                 ! Critical flux of O3 to vegetation (nmol/m2/s).
,dfp_dcuo(:)     ! Plant type specific O3 sensitivity parameter 
                 ! (nmol-1 m2 s).

!-----------------------------------------------------------------------------
! Parameters for BVOC emissions
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 ci_st(:),                                                                    &
                 ! Internal CO2 partial pressure (Pa)
                 !   at standard conditions
 gpp_st(:),                                                                   &
                 ! Gross primary productivity (KgC/m2/s)
                 !   at standard conditions
 ief(:),                                                                      &
                 ! Isoprene Emission Factor (ugC/g/h)
                 ! See Pacifico et al., (2011) Atm. Chem. Phys.
 tef(:),                                                                      &
                 ! (Mono-)Terpene Emission Factor (ugC/g/h)
 mef(:),                                                                      &
                 ! Methanol Emission Factor (ugC/g/h)
 aef(:)
                 ! Acetone Emission Factor (ugC/g/h)

!-----------------------------------------------------------------------------
! Parameters for INFERNO combustion
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 ccleaf_min(:),                                                               &
                 ! Leaf minimum combustion completeness (kg/kg)
 ccleaf_max(:),                                                               &
                 ! Leaf maximum combustion completeness (kg/kg)
 ccwood_min(:),                                                               &
                 ! Wood (or Stem) minimum combustion completeness (kg/kg)
 ccwood_max(:),                                                               &
                 ! Wood (or Stem) maximum combustion completeness (kg/kg)
 avg_ba(:),                                                                   &
                 ! Average PFT Burnt Area per fire
 fire_mort(:)
                 ! Fire mortality per PFT

!-----------------------------------------------------------------------------
! Parameters for INFERNO emissions
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
 fef_co2(:),                                                                  &
                 ! Fire CO2 Emission Factor (g/kg)
                 ! See Thonicke et al., (2005,2010)
 fef_co(:),                                                                   &
                 ! Fire CO Emission Factor (g/kg)
 fef_ch4(:),                                                                  &
                 ! Fire CH4 Emission Factor (g/kg)
 fef_nox(:),                                                                  &
                 ! Fire NOx Emission Factor (g/kg)
 fef_so2(:),                                                                  &
                 ! Fire SO2 Emission Factor (g/kg)
 fef_oc(:),                                                                   &
                 ! Fire OC Emission Factor (g/kg)
 fef_bc(:)
                 ! Fire BC Emission Factor (g/kg)

END MODULE pftparm
