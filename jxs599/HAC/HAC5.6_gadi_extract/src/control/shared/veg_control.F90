! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Top-level control routine for vegetation section
!
! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt

MODULE veg_control_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG_CONTROL_MOD'

CONTAINS
SUBROUTINE veg_control(                                                       &
  land_pts, nsurft, dim_cs1,                                                  &
  a_step, asteps_since_triffid,                                               &
  land_pts_trif, npft_trif,                                                   &
  phenol_period, triffid_period,                                              &
  l_phenol, l_triffid, l_trif_eq,                                             &
  timestep_real, frac_agr_gb, frac_past_gb, satcon_soilt,                     &
  g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,                           &
  resp_s_acc_gb_um, resp_w_acc_pft,                                           &
  cs_pool_gb_um, frac_surft, lai_pft, clay_soilt, z0m_soil_gb, canht_pft,     &
  catch_snow_surft, catch_surft, infil_surft, z0_surft, z0h_bare_surft,       &
  c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb, g_leaf_day_pft, g_leaf_phen_pft,  &
  lai_phen_pft, g_leaf_dr_out_pft, npp_dr_out_pft, resp_w_dr_out_pft,         &
  resp_s_dr_out_gb_um, qbase_l_soilt, sthf_soilt, sthu_soilt, w_flux_soilt,   &
  t_soil_soilt,cs_pool_soilt                                                  &
  )

!Module imports

!Import relevant subroutines
USE veg1_mod,                     ONLY: veg1
USE veg2_mod,                     ONLY: veg2
USE veg_soil_index_mod,           ONLY: get_veg_soil_index

#if ! defined(UM_JULES)
USE soil_biogeochem_control_mod,  ONLY: soil_biogeochem_control
#endif

!Import variables. Modules in alphabetical order
USE ancil_info,                   ONLY: dim_cslayer, nsoilt

USE conversions_mod,              ONLY: rsec_per_day, isec_per_day

USE jules_soil_biogeochem_mod,    ONLY: soil_model_ecosse, soil_bgc_model

USE jules_soil_mod,               ONLY: sm_levels

USE jules_surface_types_mod,      ONLY: npft, ntype, nnpft, soil

USE jules_vegetation_mod,         ONLY:                                       &
  i_veg_vn, i_veg_vn_1b, i_veg_vn_2b, frac_min

#if ! defined(UM_JULES)
USE model_time_mod,               ONLY: timestep_len
#endif

USE timestep_mod,                 ONLY: timestep

USE trifctl,                      ONLY: resp_s_acc_soilt, resp_s_dr_out_gb

USE trif_vars_mod,                ONLY:                                       &
  cnsrv_carbon_veg2_gb, cnsrv_veg_triffid_gb, cnsrv_soil_triffid_gb,          &
  cnsrv_prod_triffid_gb, cnsrv_carbon_triffid_gb,                             &
  cnsrv_vegN_triffid_gb, cnsrv_soilN_triffid_gb,                              &
  cnsrv_N_inorg_triffid_gb, cnsrv_nitrogen_triffid_gb,                        &
  deposition_n_gb

!Dr Hook etc
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!Arguments

INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  nsurft,                                                                     &
  dim_cs1,                                                                    &
  a_step,                                                                     &
  phenol_period,                                                              &
  triffid_period,                                                             &
  land_pts_trif,                                                              &
  npft_trif

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  timestep_real,                                                              &
  frac_agr_gb(land_pts),                                                      &
  satcon_soilt(land_pts,nsoilt,0:sm_levels),                                  &
    ! Saturated hydraulic conductivity of the soil surface (kg/m2/s).
  z0m_soil_gb(land_pts),                                                      &
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                 &
        ! Base flow from each soil layer (kg m-2 s-1).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
  w_flux_soilt(land_pts,nsoilt,0:sm_levels),                                  &
        ! Fluxes of water between layers (kg m-2 s-1).
  t_soil_soilt(land_pts,nsoilt,sm_levels)

LOGICAL, INTENT(IN) ::                                                        &
  l_phenol, l_triffid, l_trif_eq

INTEGER, INTENT(INOUT) ::                                                     &
  asteps_since_triffid

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  g_leaf_acc_pft(land_pts,npft),                                              &
  g_leaf_phen_acc_pft(land_pts,npft),                                         &
  npp_acc_pft(land_pts_trif,npft_trif),                                       &
  resp_s_acc_gb_um(land_pts_trif,dim_cs1),                                    &
  resp_w_acc_pft(land_pts_trif,npft_trif),                                    &
  cs_pool_gb_um(land_pts,dim_cs1),                                            &
  frac_surft(land_pts,ntype),                                                 &
  lai_pft(land_pts,npft),                                                     &
  frac_past_gb(land_pts),                                                     &
  canht_pft(land_pts,npft),                                                   &
  clay_soilt(land_pts,nsoilt,dim_cslayer),                                    &
  cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  catch_snow_surft(land_pts,nsurft),                                          &
  catch_surft(land_pts,nsurft),                                               &
  infil_surft(land_pts,nsurft),                                               &
  z0_surft(land_pts,nsurft),                                                  &
  z0h_bare_surft(land_pts,nsurft),                                            &
  c_veg_pft(land_pts,npft),                                                   &
  cv_gb(land_pts),                                                            &
  lit_c_pft(land_pts,npft),                                                   &
  lit_c_mn_gb(land_pts),                                                      &
  g_leaf_day_pft(land_pts,npft),                                              &
  g_leaf_phen_pft(land_pts,npft),                                             &
  lai_phen_pft(land_pts,npft),                                                &
  g_leaf_dr_out_pft(land_pts,npft),                                           &
  npp_dr_out_pft(land_pts,npft),                                              &
  resp_w_dr_out_pft(land_pts,npft),                                           &
  resp_s_dr_out_gb_um(land_pts,dim_cs1+1)

!Local variables

INTEGER ::                                                                    &
  l, n, m, trif_pts, phenol_call, triffid_call, nstep_trif,                   &
  trif_index(land_pts)

REAL(KIND=real_jlslsm) ::                                                     &
  resp_s_dr_out(land_pts,1,dim_cs1+1),                                        &
  frac_surft_start(land_pts,ntype),                                           &
         ! Fractions of surface types at the start of the timestep.
  frac_vs(land_pts)
        ! Fraction of gridbox covered by veg or soil.

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER                           :: errorstatus

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VEG_CONTROL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! initialize carbon conservation diagnostics
! otherwise they can be non-zero on non-triffid timesteps
IF ( l_triffid ) THEN
  DO l = 1, land_pts
    cnsrv_carbon_veg2_gb(l)      = 0.0
    cnsrv_carbon_triffid_gb(l)   = 0.0
    cnsrv_veg_triffid_gb(l)      = 0.0
    cnsrv_soil_triffid_gb(l)     = 0.0
    cnsrv_prod_triffid_gb(l)     = 0.0
    cnsrv_nitrogen_triffid_gb(l) = 0.0
    cnsrv_vegN_triffid_gb(l)     = 0.0
    cnsrv_soilN_triffid_gb(l)    = 0.0
    cnsrv_N_inorg_triffid_gb(l)  = 0.0
  END DO
END IF

!Determine triggering for the two veg schemes.
!There is a subtle difference between standalone and UM behaviour
phenol_call  = 1
triffid_call = 1

IF (l_phenol) THEN
#if defined(UM_JULES)
  phenol_call = MOD(REAL(a_step),                                             &
                    (REAL(phenol_period) * (rsec_per_day / timestep)))
#else
  phenol_call = MOD (a_step,                                                  &
                          phenol_period * isec_per_day / timestep_len)
#endif
END IF

IF (l_triffid) THEN
#if defined(UM_JULES)
  nstep_trif = INT(rsec_per_day *       triffid_period / timestep)
#else
  nstep_trif = INT(rsec_per_day * REAL(triffid_period) / timestep_real )
#endif
  IF ( asteps_since_triffid == nstep_trif ) THEN
    triffid_call = 0
  END IF
END IF


#if defined(UM_JULES)

cs_pool_soilt(:,1,1,:)    = cs_pool_gb_um(:,:)
resp_s_acc_soilt(:,1,1,:) = resp_s_acc_gb_um(:,:)
resp_s_dr_out(:,1,:)      = 0.0
!Can't copy from resp_s_dr_out_gb_um as it's INTENT(OUT) so set to zero

!Call to veg 1 or 2 as appropriate
IF ((phenol_call == 0) .OR. (triffid_call == 0)) THEN
  SELECT CASE ( i_veg_vn )
  CASE ( i_veg_vn_2b )
    !     ! Find total fraction of gridbox covered by vegetation and soil, and
    !     !  use this to set indices of land points on which TRIFFID may operate.

    CALL get_veg_soil_index( land_pts, frac_surft, trif_pts,                  &
                            trif_index, frac_vs )

    ! Running with nsoilt > 1 is not compatible with dynamic vegetation, so we can
    ! hard-code the arguments appropriately.
    CALL veg2( land_pts, nsurft, a_step,                                      &
               phenol_period, triffid_period,                                 &
               trif_pts, trif_index, timestep_real,                           &
               frac_agr_gb, frac_past_gb, frac_vs,                            &
               satcon_soilt(:,1,0), clay_soilt(:,1,:), z0m_soil_gb,           &
               l_phenol, l_triffid, l_trif_eq,                                &
               asteps_since_triffid,                                          &
               g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,              &
               resp_s_acc_soilt(:,1,1,:), resp_w_acc_pft,                     &
               cs_pool_soilt(:,1,1,:), frac_surft, lai_pft, canht_pft,        &
               catch_snow_surft, catch_surft, infil_surft,                    &
               z0_surft, z0h_bare_surft, c_veg_pft, cv_gb,                    &
               g_leaf_day_pft, g_leaf_phen_pft, g_leaf_dr_out_pft,            &
               lai_phen_pft, lit_c_pft, lit_c_mn_gb, npp_dr_out_pft,          &
               resp_w_dr_out_pft, resp_s_dr_out                               &
               )

  CASE ( i_veg_vn_1b )
    CALL veg1( land_pts, nsurft, a_step, phenol_period, timestep_real,        &
               satcon_soilt(:,1,0), z0m_soil_gb, l_phenol,                    &
               g_leaf_acc_pft, g_leaf_phen_acc_pft, frac_surft, lai_pft,      &
               canht_pft,                                                     &
               catch_snow_surft, catch_surft, infil_surft,                    &
               g_leaf_day_pft, g_leaf_phen_pft,                               &
               lai_phen_pft, z0_surft, z0h_bare_surft                         &
               )

  CASE DEFAULT ! i_veg_vn
    errorstatus = 10
    WRITE (cmessage,'(A,A,I6)') 'Vegetation scheme version value',            &
           'i_veg_vn = ',i_veg_vn
    CALL Ereport ('VEG_CTL', errorstatus, cmessage)

  END SELECT ! i_veg_vn
END IF

cs_pool_gb_um(:,:)       = cs_pool_soilt(:,1,1,:)
resp_s_acc_gb_um(:,:)    = resp_s_acc_soilt(:,1,1,:)
resp_s_dr_out_gb_um(:,:) = resp_s_dr_out(:,1,:)

#else

IF ( l_triffid ) THEN

  IF ( triffid_call == 0 .OR. soil_bgc_model == soil_model_ecosse ) THEN

    ! Find total fraction of gridbox covered by vegetation and soil, and use
    ! this to set indices of land points on which TRIFFID may operate.
    ! We also do this if ECOSSE is used so the veg and soil models operate
    ! on the same set of points.
    ! Note: This code is essentially a repeat of code in veg_control. That
    ! subroutine is expected to eventually also be used for standalone
    ! JULES, but in the meanwhile we also calculate trif_pts here.
    CALL get_veg_soil_index( land_pts, frac_surft, trif_pts,                  &
                            trif_index, frac_vs )
  END IF  !  triffid_call OR ecosse
END IF  !  l_triffid

! Save frac at start of timestep (for ECOSSE).
frac_surft_start(:,:) = frac_surft(:,:)

IF ( triffid_call == 0 ) THEN
  !     Run includes dynamic vegetation
  !
  ! Running with nsoilt > 1 is not compatible with dynamic vegetation, so we
  ! can hard-code the arguments appropriately by setting m = 1

  m = 1
  CALL veg2( land_pts, nsurft, a_step                                         &
            ,phenol_period, triffid_period                                    &
            ,trif_pts, trif_index, timestep_real                              &
            ,frac_agr_gb, frac_past_gb, frac_vs                               &
            ,satcon_soilt(:,:,0), clay_soilt(:,m,:), z0m_soil_gb              &
            ,l_phenol, l_triffid, l_trif_eq                                   &
            ,asteps_since_triffid                                             &
            ,g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft                 &
            ,resp_s_acc_soilt(:,m,:,:), resp_w_acc_pft                        &
            ,cs_pool_soilt(:,m,:,:), frac_surft, lai_pft, canht_pft           &
            ,catch_snow_surft, catch_surft, infil_surft                       &
            ,z0_surft, z0h_bare_surft, c_veg_pft, cv_gb                       &
            ,g_leaf_day_pft, g_leaf_phen_pft, g_leaf_dr_out_pft               &
            ,lai_phen_pft, lit_c_pft, lit_c_mn_gb, npp_dr_out_pft             &
            ,resp_w_dr_out_pft, resp_s_dr_out_gb )

ELSE

  IF ( phenol_call == 0 ) THEN
    ! Run includes phenology,  but not dynamic vegetation
    ! therefore call veg1 rather than veg2
    CALL veg1( land_pts, nsurft, a_step, phenol_period, timestep_real         &
              ,satcon_soilt(:,:,0), z0m_soil_gb, l_phenol                     &
              ,g_leaf_acc_pft, g_leaf_phen_acc_pft, frac_surft, lai_pft       &
              ,canht_pft, catch_snow_surft, catch_surft, infil_surft          &
              ,g_leaf_day_pft, g_leaf_phen_pft                                &
              ,lai_phen_pft, z0_surft, z0h_bare_surft )
  END IF

END IF  !  triffid_call
#endif


#if !defined(UM_JULES)
!Soil biogeochemistry (Standalone only)
! At present this only deals with the ECOSSE model of soil C and N.
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  CALL soil_biogeochem_control( land_pts, triffid_call, trif_pts,             &
        trif_index, deposition_n_gb, frac_surft_start,                        &
        qbase_l_soilt, sthf_soilt, sthu_soilt, w_flux_soilt, t_soil_soilt )
END IF
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE veg_control
END MODULE veg_control_mod