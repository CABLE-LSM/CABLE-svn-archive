! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holds surface parameters for each Plant Functional Type (but not parameters that are
! only used by TRIFFID or phenology.).



! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE pftparm

  IMPLICIT NONE

! If we are not in the UM, define all the allocatable arrays
  INTEGER, ALLOCATABLE ::                                           &
   c3(:)                                                            &
                   ! Flag for C3 types: 1 for C3 Plants,
!                    0 for C4 Plants.
  ,orient(:)
                   ! Flag for leaf orientation: 1 for horizontal,
!                    0 for spherical.

  REAL, ALLOCATABLE ::                                              &
   a_wl(:)                                                          &
                   ! Allometric coefficient relating the target
!                    woody biomass to the leaf area index
!                    (kg C/m2)
  ,a_ws(:)                                                          &
                   ! Woody biomass as a multiple of live
!                    stem biomass.
  ,albsnc_max(:)                                                    &
                   ! Snow-covered albedo for large LAI.
  ,albsnc_min(:)                                                    &
                   ! Snow-covered albedo for zero LAI.
  ,albsnf_max(:)                                                    &
                   ! Snow-free albedo for large LAI.
  ,alpha(:)                                                         &
                   ! Quantum efficiency
!                    (mol CO2/mol PAR photons).
  ,alnir(:)                                                         &
                   ! Leaf reflection coefficient for near
!                    infra-red.
  ,alpar(:)                                                         &
                   ! Leaf reflection coefficient for PAR.
  ,b_wl(:)                                                          &
                   ! Allometric exponent relating the target
!                    woody biomass to the leaf area index.
  ,catch0(:)                                                        &
                   ! Minimum canopy capacity (kg/m2).
  ,dcatch_dlai(:)                                                   &
                   ! Rate of change of canopy capacity with LAI.
  ,dgl_dm(:)                                                        &
                   ! Rate of change of leaf turnover rate with
!                    moisture availability.
  ,dgl_dt(:)                                                        &
                   ! Rate of change of leaf turnover rate with
!                    temperature (/K)
  ,dqcrit(:)                                                        &
                   ! Critical humidity deficit (kg H2O/kg air)
  ,dz0v_dh(:)                                                       &
                   ! Rate of change of vegetation roughness
!                    length with height.
  ,eta_sl(:)                                                        &
                   ! Live stemwood coefficient (kg C/m/LAI)
  ,fd(:)                                                            &
                   ! Dark respiration coefficient.
  ,fsmc_of(:)                                                       &
                   ! Moisture availability below which leaves
!                    are dropped.
  ,f0(:)                                                            &
                   !  CI/CA for DQ = 0..
  ,g_leaf_0(:)                                                      &
                   ! Minimum turnover rate for leaves (/360days).
  ,glmin(:)                                                         &
                   ! Minimum leaf conductance for H2O.
  ,infil_f(:)                                                       &
                   ! Infiltration enhancement factor.
  ,kext(:)                                                          &
                   ! Light extinction coefficient - used to
!                    calculate weightings for soil and veg.
  ,kpar(:)                                                          &
                   ! PAR Extinction coefficient
!                    (m2 leaf/m2 ground)
  ,neff(:)                                                          &
                  ! Constant relating VCMAX and leaf N (mol/m2/s)
!                   from Schulze et al. 1994
!                   (AMAX = 0.4E-3 * NL  - assuming dry matter is
!                   40% carbon by mass)
!                   and Jacobs 1994:
!                   C3 : VCMAX = 2 * AMAX ;
!                   C4 : VCMAX = AMAX  ..
  ,nl0(:)                                                           &
                   ! Top leaf nitrogen concentration
!                    (kg N/kg C).
  ,nr_nl(:)                                                         &
                   ! Ratio of root nitrogen concentration to
!                    leaf nitrogen concentration.
  ,ns_nl(:)                                                         &
                   ! Ratio of stem nitrogen concentration to
!                    leaf nitrogen concentration.
  ,omega(:)                                                         &
                   ! Leaf scattering coefficient for PAR.
  ,omnir(:)                                                         &
                   ! Leaf scattering coefficient for near
!                    infra-red.
  ,r_grow(:)                                                        &
                   ! Growth respiration fraction.
  ,rootd_ft(:)                                                      &
                   ! e-folding depth (m) of the root density.
  ,sigl(:)                                                          &
                   ! Specific density of leaf carbon
!                    (kg C/m2 leaf).
  ,tleaf_of(:)                                                      &
                   ! Temperature below which leaves are
!                    dropped (K)
  ,tlow(:)                                                          &
                   ! Lower temperature for photosynthesis (deg C)
  ,tupp(:)                                                          &
                   !  Upper temperature for photosynthesis (deg C)
  ,emis_pft(:)                                                      &
                   !  Surface emissivity
  ,z0hm_pft(:)                                                      &
                   ! z0h/z0m for PFTs
  ,dust_veg_scj(:) ! Dust emission scaling factor for  each PFT


! Parameters for ozone damage
  REAL, ALLOCATABLE ::                                              &
   fl_o3_ct(:)                                                      &
                   ! Critical flux of O3 to vegetation (nmol/m2/s).
  ,dfp_dcuo(:)     ! Fractional reduction of photosynthesis
                   ! with the cumulative uptake of O3 by
                   ! leaves (/mmol/m2).

! Parameters for BVOC emissions
  REAL, ALLOCATABLE ::                                              &
   ief(:),                                                          &
                   ! Isoprene Emission Factor (ugC/g/h)
                   ! See Pacifico et al., (2011) Atm. Chem. Phys.
   tef(:),                                                          &
                   ! (Mono-)Terpene Emission Factor (ugC/g/h)
   mef(:),                                                          &
                   ! Methanol Emission Factor (ugC/g/h)
   aef(:)
                   ! Acetone Emission Factor (ugC/g/h)


  CHARACTER(LEN=30), ALLOCATABLE ::                                 &
   pftname(:)        !  Name of each PFT

END MODULE pftparm
