#if !defined(UM_JULES)

SUBROUTINE init_vars_tmp()

USE Ancil_info, ONLY: ice_fract_ncat_sicat, ti_cat_sicat, lice_pts, lice_index
USE boundary_layer_var_mod, ONLY: dzl, zh
USE jules_deposition_mod, ONLY: dzl_const, l_deposition
USE jules_surface_mod, ONLY: cmass
USE top_pdm, ONLY: inlandout_atm_gb
USE Forcing, ONLY: u_0_ij,v_0_ij
USE Prognostics, ONLY: canht_pft, di_ncat_sicat, k_sice_sicat, lai_pft,       &
                        snow_mass_ij, snow_mass_sea_sicat, soot_ij, z0msea_ij
USE Aero
USE Orog
USE jules_vegetation_mod, ONLY: l_triffid, l_use_pft_psi, fsmc_shape,         &
                                l_phenol

USE Trifctl, ONLY: cv_gb, g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,   &
                   resp_s_acc_soilt, resp_w_acc_pft, c_veg_pft, cv_gb

USE trif_vars_mod, ONLY: leafc_pft, rootc_pft, woodc_pft, n_veg_pft,          &
                          n_leaf_pft, n_root_pft, n_stem_pft, lai_bal_pft,    &
                          n_veg_gb, wp_fast_out_gb, wp_med_out_gb,            &
                          wp_slow_out_gb,wp_fast_in_gb, wp_med_in_gb,         &
                          wp_slow_in_gb

USE trif, ONLY: lai_min

USE pftparm, ONLY: fsmc_mod, psi_close, psi_open

USE ancil_info, ONLY: land_pts, frac_surft, nsoilt
USE Coastal
USE jules_sea_seaice_mod
USE C_kappai
USE p_s_parms, ONLY: satcon_soilt, smvcst_soilt, smvcwt_soilt,                &
                      v_close_pft, v_open_pft, sathh_soilt,                   &
                      smvccl_soilt, bexp_soilt

USE hyd_psi_mod, ONLY: sthu_from_psi

USE jules_sea_seaice_mod, ONLY: l_ssice_albedo

USE calc_c_comps_triffid_mod, ONLY: calc_c_comps_triffid

USE CN_utils_mod, ONLY: calc_n_comps_triffid

USE update_mod, ONLY: bl_height, l_imogen

USE logging_mod, ONLY: log_info, log_debug, log_warn, log_error, log_fatal

USE jules_surface_types_mod, ONLY: npft, ncpft, nnpft

USE jules_soil_mod, ONLY: sm_levels, soil_props_const_z

USE dump_mod, ONLY: ancil_dump_read

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises various variables that may change their initialisation in
!   future versions
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
REAL(KIND=real_jlslsm) :: phen ! Phenological state (=LAI/lai_bal)

INTEGER :: i,l,n,ft,m  ! Loop counters


!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Initialise accumulated fluxes for TRIFFID and phenology.
! This is not necessary if these are read from a restart file - but at
! present they're not.
!-----------------------------------------------------------------------
g_leaf_acc_pft(:,:)       = 0.0
g_leaf_phen_acc_pft(:,:)  = 0.0
npp_acc_pft(:,:)          = 0.0
resp_s_acc_soilt(:,:,:,:) = 0.0
resp_w_acc_pft(:,:)       = 0.0

!-----------------------------------------------------------------------
! Set saturated hydraulic conductivity to zero at land ice points.
!-----------------------------------------------------------------------
IF ( lice_pts > 0 ) THEN
  CALL log_info("init_vars_tmp",                                              &
                "Setting satcon to zero at land ice points")
  DO i = 1,lice_pts
    l = lice_index(i)
    satcon_soilt(l,:,:) = 0.0
  END DO
END IF

!-----------------------------------------------------------------------
! Set surface velocity to be zero
!-----------------------------------------------------------------------
u_0_ij(:,:) = 0.0
v_0_ij(:,:) = 0.0

!-----------------------------------------------------------------------
! Set boundary layer height and layer separation (if needed).
! These might later be updated via prescribed data.
!-----------------------------------------------------------------------
zh(:,:) = bl_height
IF ( l_deposition ) dzl(:,:,:) = dzl_const

!-----------------------------------------------------------------------
! Set CO2 variables
!-----------------------------------------------------------------------
co2_3d_ij(:,:) = 0.0

!-----------------------------------------------------------------------
! Set coastal tiling variables
!-----------------------------------------------------------------------
tstar_sea_ij(:,:)       = 280.0
tstar_sice_sicat(:,:,:) = 270.0

!-----------------------------------------------------------------------
! Set orographic roughness variables
!-----------------------------------------------------------------------
h_blend_orog_ij(:,:) = 0.0
sil_orog_land_gb(:)  = 0.0
ho2r2_orog_gb(:)     = 0.0

!-----------------------------------------------------------------------
! Set up prognostics which are not currently in dump
!-----------------------------------------------------------------------
ti_cat_sicat(:,:,:)         = 270.0
z0msea_ij(:,:)              = z0hsea
snow_mass_ij(:,:)           = 0.0
ice_fract_ncat_sicat(:,:,:) = 0.0
di_ncat_sicat(:,:,:)        = 0.0
snow_mass_sea_sicat(:,:,:)  = 0.0
k_sice_sicat(:,:,:)         = 2.0 * kappai / de

!-----------------------------------------------------------------------
! Set up sea-ice parameter variables
!-----------------------------------------------------------------------
IF (l_ssice_albedo) THEN
  alpham = 0.65
  alphac = 0.80
  alphab = 0.57
  dtice = 2.00
ELSE
  alpham = 0.50
  alphac = 0.80
  alphab=-1.00
  dtice = 10.00
END IF

soot_ij(:,:) = 0.0

!-----------------------------------------------------------------------------
! Unless IMOGEN is on, cv will not have been initialised
! Unless TRIFFID is on, it will also not be used...
!-----------------------------------------------------------------------------
IF ( .NOT. l_imogen ) cv_gb(:) = 0.0

! Initialise c_veg to 0 as not doing so can cause issues with IMOGEN outputs
c_veg_pft(:,:) = 0.0

!-----------------------------------------------------------------------------
! Initialising TRIFFID Diagnostic Veg Pools and fluxes prior to first TRIFFID call
! This is necessary to ensure time average diagnostics are produced correctly.
!-----------------------------------------------------------------------------

IF (l_triffid) THEN
  rootC_pft       = 0.0
  woodC_pft       = 0.0
  leafC_pft       = 0.0
  c_veg_pft       = 0.0
  n_leaf_pft      = 0.0
  n_root_pft      = 0.0
  n_stem_pft      = 0.0
  lai_bal_pft     = 0.0
  cv_gb           = 0.0
  n_veg_gb        = 0.0
  wp_fast_out_gb  = 0.0
  wp_med_out_gb   = 0.0
  wp_slow_out_gb  = 0.0
  wp_fast_in_gb   = 0.0
  wp_med_in_gb    = 0.0
  wp_slow_in_gb   = 0.0

  IF ( .NOT. l_phenol) THEN
    DO n = 1,nnpft
      lai_pft(:,n) = lai_min(n)
    END DO
  END IF

  DO l = 1,land_pts
    DO n = 1,nnpft
      IF (frac_surft(l,n) > 0.0) THEN

        CALL calc_c_comps_triffid(n, canht_pft(l,n), lai_bal_pft(l,n),        &
                                  leafC_pft(l,n), rootC_pft(l,n),             &
                                  woodC_pft(l,n), c_veg_pft(l,n))

        IF (l_phenol) THEN
          phen = lai_pft(l,n) / lai_bal_pft(l,n)
          IF ( phen > 1.0 + TINY(1.0e0) ) THEN
            CALL log_warn("init_vars_tmp",                                    &
                          "lai_pft should be <= lai_bal_pft")
          END IF
        ELSE
          phen = 1.0
          lai_pft(l,n) = lai_bal_pft(l,n)
        END IF

        CALL calc_n_comps_triffid(l ,n, phen, lai_bal_pft(l,n),               &
                                  woodC_pft(l,n), rootC_pft(l,n),             &
                                  n_leaf_pft(l,n), n_root_pft(l,n), n_stem_pft(l,n))

        n_veg_pft(l,n) = n_leaf_pft(l,n) + n_root_pft(l,n) + n_stem_pft(l,n)

        cv_gb(l) = cv_gb(l) + frac_surft(l,n) * c_veg_pft(l,n)

        n_veg_gb(l) = n_veg_gb(l) + frac_surft(l,n) * n_veg_pft(l,n)
      END IF
    END DO
  END DO
END IF

!------------------------------------------------------------------
! Temporary initialisation of variables added during UM integration
!------------------------------------------------------------------

inlandout_atm_gb(:) = 0.0

!------------------------------------------------------------------
! Initialisation of v_close_soilt_pft and v_open_soilt_pft
!------------------------------------------------------------------

IF ( l_use_pft_psi ) THEN
  ! now set v_close_soilt_pft and v_open_soilt_pft
  DO ft = 1,npft

    !Set the current soil tile (see *NOTICE REGARDING SOIL TILING* in physiol)
    IF (nsoilt == 1) THEN
      !There is only 1 soil tile
      m = 1
    ELSE ! nsoilt == nsurft
      !Soil tiles map directly on to surface tiles
      m = ft
    END IF !nsoilt

    DO n = 1,sm_levels
      DO i = 1,land_pts
        v_close_pft(i,n,ft) = sthu_from_psi(psi_close(ft),                    &
                              sathh_soilt(i,m,n),bexp_soilt(i,m,n),           &
                              0.01) * smvcst_soilt(i,m,n)
        v_open_pft(i,n,ft) = sthu_from_psi(psi_open(ft),                      &
                             sathh_soilt(i,m,n), bexp_soilt(i,m,n),           &
                             0.01) * smvcst_soilt(i,m,n)
      END DO
    END DO
  END DO
END IF

IF ((fsmc_shape == 1) .AND. ANY(fsmc_mod == 1)) THEN
  IF ( ancil_dump_read%soil_props ) THEN
    CALL log_fatal("init_vars_tmp",                                           &
                   "fsmc_shape=fsmc_mod=1 can not currently be used " //      &
                   "when JULES_SOIL_PROPS has read_from_dump=T.")
    ! FIXME: could allow this if the values in the soil ancils are checked
  END IF
  IF ( .NOT. soil_props_const_z ) THEN
    CALL log_fatal("init_vars_tmp",                                           &
                   "fsmc_shape=fsmc_mod=1 requires " //                       &
                   "const_z=T in JULES_SOIL_PROPS.")
  END IF
END IF

RETURN

END SUBROUTINE init_vars_tmp
#endif
