! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE crop_vars_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module holding various variables for the crop model
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
!-----------------------------------------------------------------------------
! Variables
!-----------------------------------------------------------------------------
    sow_date_cpft(:,:),                                                       &
        ! Sowing date of each crop functional type
    latestharv_date_cpft(:,:),                                                &
        ! Sowing date of each crop functional type
    tt_veg_cpft(:,:),                                                         &
        ! Thermal requirement of stage 1 for crop pfts (degree days).
    tt_rep_cpft(:,:),                                                         &
        ! Thermal requirement of stage 2 for crop pfts (degree days).
    phot(:),                                                                  &
      !  Photoperiod (hours) for crop model
    dphotdt(:),                                                               &
      !  Rate of change of phot for crops (hours per day).
!-----------------------------------------------------------------------------
! Prognostics
!-----------------------------------------------------------------------------
    dvi_cpft(:,:),                                                            &
                !  Development index for crop tiles
    rootc_cpft(:,:),                                                          &
                !  Root carbon pool for crop tiles (kg m-2).
    harvc_cpft(:,:),                                                          &
                !  Carbon in 'harvest parts' pool for crop tiles (kg m-2).
    reservec_cpft(:,:),                                                       &
                !  Carbon in stem reserves pool for crop tiles (kg m-2).
    croplai_cpft(:,:),                                                        &
                !  Leaf area index for crop tiles
    cropcanht_cpft(:,:),                                                      &
                !  Canopy height for crop tiles (m).
    sthu_irr_soilt(:,:,:),                                                    &
                !  Unfrozen soil wetness over irrigation
!-----------------------------------------------------------------------------
! Diagnostics
!-----------------------------------------------------------------------------
    yield_diag_cpft(:,:),                                                     &
        ! Harvested carbon (kg m-2).
    stemc_diag_cpft(:,:),                                                     &
        ! Stem carbon (kg m-2).
    leafc_diag_cpft(:,:),                                                     &
        ! Leaf carbon (kg m-2).
    irrig_water_gb(:),                                                        &
        ! Addition of irrigation water to soil (kg/m2/s).
    nonyield_diag_cpft(:,:)
        ! Carbon leaving the crop model which is not yield (kg m-2).

INTEGER, ALLOCATABLE ::                                                       &
  harvest_trigger_cpft(:,:),                                                  &
      ! Trigger condition for the harvest in this timestep
  harvest_counter_cpft(:,:)
      ! 1 if timestep contains a harvest, 0 otherwise

!-----------------------------------------------------------------------------
! Irrigation variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  frac_irr_all(:,:),                                                          &
  frac_irr_soilt(:,:),                                                        &
      ! Irrigation fraction for each soil tile
  frac_irr_old_soilt(:,:)
      ! Previous irrigated fraction in grid box

LOGICAL :: frac_irr_all_tiles
      ! Switch to assign irrigation fraction to ALL tiles, or to pre-defined
      ! tiles ONLY
      ! set TRUE to reproduce original results

INTEGER :: nirrtile
      ! Nr of tiles that can have irrigated fraction
INTEGER, ALLOCATABLE :: irrtiles(:)
      ! Tiles that can have irrigated fraction
      ! Only used when frac_irr_all_tiles = .FALSE.

! Constants
INTEGER, PARAMETER ::                                                         &
  ndpy = 365,                                                                 &
      ! No. of days per year
  nyav = 3,                                                                   &
      ! No. of years averaged for crop plant estimates
  nday_crop = 150
      ! Cropping date

INTEGER :: iyear_old, startyr, startmon, startday, starttime
      ! Variables to store required datetime information

INTEGER ,ALLOCATABLE ::                                                       &
  plant_n_gb(:),                                                              &
      ! best plant date for non-rice
  icntmax_gb(:)
      ! counter for start date for non-rice

! Variables for calculating average temperature, precip and radiation
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  tl_1_day_av_gb(:),                                                          &
    ! Average air temperature for the current day (K).
  tl_1_day_av_use_gb(:,:,:),                                                  &
    ! Daily average air temperature (K).
  prec_1_day_av_gb(:),                                                        &
    ! Average precipitation rate for the current day (kg m-2 s-1).
  prec_1_day_av_use_gb(:,:,:),                                                &
    ! Daily average precipitation rate (kg m-2 s-1).
  rn_1_day_av_gb(:),                                                          &
    ! Average net radiation for the current day (W m-2).
  rn_1_day_av_use_gb(:,:,:)
    ! Daily average net radiation (W m-2).

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  frac_irr_surft(:,:),                                                        &
      ! Irrigation fraction in tile
  smc_irr_soilt(:,:),                                                         &
      ! Available moisture in the soil profile in irrig frac (mm)
  wt_ext_irr_surft(:,:,:),                                                    &
      ! Fraction of transpiration over irrigated fraction which is
      ! extracted from each soil layer by each tile
  gs_irr_surft(:,:),                                                          &
      ! Conductance for irrigated surface types (m s-1).
  dvimax_gb(:),                                                               &
      ! Maximum value of development index in the gridbox.
  gc_irr_surft(:,:),                                                          &
      ! Interactive canopy conductance (m s-1).
  resfs_irr_surft(:,:),                                                       &
      ! Combined soil, stomatal and aerodynamic resistance factor for
      ! fraction 1-FRACA
  ext_irr_soilt(:,:,:),                                                       &
      ! Extraction of water from each soil layer over irrigation
      ! (kg m-2 s-1).
  wt_ext_irr_gb(:,:),                                                         &
      ! Fraction of transpiration extracted from each soil layer
      ! over irrigated fraction.
  fsmc_irr_gb(:),                                                             &
      ! Soil moisture availability factor over irrigated fraction
  irrDaysDiag_gb(:)
      ! Number of days on which irrigation is applied

END MODULE crop_vars_mod
