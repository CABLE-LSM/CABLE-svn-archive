#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology, 2017.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE soil_biogeochem_control_mod

!-----------------------------------------------------------------------------
! Description:
!   Control routine for soil biogeochemistry.
!   At present this is only called for the ECOSSE soil model.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE ancil_info, ONLY:                                                         &
  ! imported scalars
  nsoilt

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC soil_biogeochem_control

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName = 'SOIL_BIOGEOCHEM_CONTROL_MOD'

CONTAINS

!#############################################################################

SUBROUTINE soil_biogeochem_control( land_pts, triffid_call, vs_pts,           &
        vs_index, deposition_n_gb, frac_surft_start,                          &
        qbase_l_soilt, sthf_soilt, sthu_soilt, w_flux_soilt, t_soil_soilt )

USE ecosse_control_mod, ONLY:                                                 &
  ! imported procedures
  ecosse_control

USE ereport_mod, ONLY:                                                        &
  ! imported procedures
  ereport

USE jules_soil_biogeochem_mod, ONLY:                                          &
  ! imported scalar parameters
  soil_model_ecosse,                                                          &
  ! imported scalars
  soil_bgc_model

USE jules_soil_mod, ONLY:                                                     &
  ! imported scalars
  sm_levels

USE jules_soil_ecosse_mod, ONLY:                                              &
  ! imported scalars
  dt_soilc, l_soil_N

USE jules_surface_types_mod, ONLY:                                            &
  ! imported scalars
  npft, ntype

USE model_time_mod, ONLY:                                                     &
  ! imported scalars
  timestep_number=>timestep

USE soil_ecosse_vars_mod, ONLY:                                               &
  ! imported arrays
  n_amm_soilt, n_nit_soilt, soil_c_add, soil_n_add

USE timestep_mod, ONLY:                                                       &
  ! imported scalars
  timestep_len_real=>timestep

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Control routine for soil biogeochemistry.
!   At present it only deals with ECOSSE but other models could be added.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! Number of land points.
  triffid_call,                                                               &
    ! Indicates if the dynamic vegetation model was called earlier in
    ! the current timestep.
    ! 0 = TRIFFID was called.
    ! 1 = TRIFFID was not called.
  vs_pts
    ! The number of points with vegetation and/or soil.

!-----------------------------------------------------------------------------
! Array arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  vs_index(land_pts)
    ! Indices of points with vegetation and/or soil.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  deposition_n_gb(land_pts),                                                  &
    ! Atmospheric deposition of N to land surface (kg m-2 s-1).
 frac_surft_start(land_pts,ntype),                                            &
    !  Fractional coverage of surface types at start of timestep.
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                 &
    ! Lateral flux of water from each soil layer (kg m-2 s-1).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content as a fraction of saturation.
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Unfrozen soil moisture content as a fraction of saturation.
  w_flux_soilt(land_pts,nsoilt,0:sm_levels),                                  &
    ! Downward water flux at bottom of each soil layer (kg m-2 s-1).
  t_soil_soilt(land_pts,nsoilt,sm_levels)
    ! Subsurface temperature in each layer (K).

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SOIL_BIOGEOCHEM_CONTROL'

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  errorstatus,                                                                &
    ! Value of error flag to indicate a hard or fatal error.
  istep,                                                                      &
    ! Work.
  m,                                                                          &
    ! Soil tile number.
  nstep
    ! Number of JULES timesteps per soil model timestep.

LOGICAL ::                                                                    &
  l_call   ! .TRUE. on timesteps when the soil model is run.

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Work out if the soil model is to be called this timestep.
!-----------------------------------------------------------------------------
nstep = NINT( dt_soilc / timestep_len_real )
istep = MOD( timestep_number, nstep )
IF ( istep == 0 ) THEN
  l_call = .TRUE.
ELSE
  l_call = .FALSE.
END IF

!-----------------------------------------------------------------------------
! Update driving variables.
!-----------------------------------------------------------------------------
CALL increment_soil_drivers( land_pts, nstep, l_call, .FALSE.,                &
                             deposition_n_gb, qbase_l_soilt, sthf_soilt,      &
                             sthu_soilt, t_soil_soilt, w_flux_soilt )

!-----------------------------------------------------------------------------
! If using ECOSSE, initialise increment for mass conservation.
! This value is retained on timesteps when ECOSSE is not called.
!-----------------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  soil_c_add(:) = 0.0
  IF ( l_soil_N ) soil_n_add(:) = 0.0
END IF

!-----------------------------------------------------------------------------
! If the soil model is not to be called this timtestep, nothing more to do.
!-----------------------------------------------------------------------------
IF ( .NOT. l_call ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                       &
                          zhook_out,zhook_handle)
  RETURN
END IF

!-----------------------------------------------------------------------------
! Call the science control routine for the chosen model.
! At present this can only be ECOSSE.
!-----------------------------------------------------------------------------
IF ( nsoilt > 1 ) THEN
  errorstatus = 100
  CALL ereport( TRIM(RoutineName), errorstatus,                               &
                "Only coded for a single soil tile." )
END IF
! Set soil tile number for inorganic fields.
m = 1
CALL ecosse_control( land_pts, triffid_call, vs_pts, vs_index,                &
                     n_amm_soilt(:,m,:), n_nit_soilt(:,m,:),                  &
                     frac_surft_start )

!-----------------------------------------------------------------------------
! Reset accumulations of driving variables.
!-----------------------------------------------------------------------------
CALL increment_soil_drivers( land_pts, nstep, l_call, .TRUE.,                 &
                             deposition_n_gb, qbase_l_soilt, sthf_soilt,      &
                             sthu_soilt, t_soil_soilt, w_flux_soilt )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE soil_biogeochem_control

!#############################################################################

SUBROUTINE increment_soil_drivers(                                            &
       land_pts, nstep, l_call, l_zero_vars,                                  &
       deposition_n_gb, qbase_l_soilt, sthf_soilt, sthu_soilt,                &
       t_soil_soilt, w_flux_soilt )

! Description:
!  Increment driver variables for a soil model.

! Module imports.
USE ancil_info, ONLY:                                                         &
  ! imported scalars
  soil_pts,                                                                   &
  ! imported arrays
  soil_index

USE jules_soil_ecosse_mod, ONLY:                                              &
  ! imported scalars
  l_driver_ave

USE soil_ecosse_vars_mod, ONLY:                                               &
  ! imported arrays
  deposition_n_driver, qbase_l_driver, sthf_driver, sthu_driver,              &
  tsoil_driver, wflux_driver

USE jules_hydrology_mod, ONLY:                                                &
  ! imported scalars
  l_top

USE jules_soil_mod, ONLY:                                                     &
  ! imported scalars
  sm_levels

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! Number of land points.
  nstep
    ! Number of JULES timesteps per soil model timestep.

LOGICAL, INTENT(IN) ::                                                        &
  l_call,                                                                     &
    ! .TRUE. on timesteps when soil_model is to be run.
  l_zero_vars
    ! TRUE when this routine is called to reset accumulations to zero.

!-----------------------------------------------------------------------------
! Array arguments with intent(in).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  deposition_n_gb(land_pts),                                                  &
    ! Atmospheric deposition of N to land surface (kg m-2 s-1).
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                 &
    ! Lateral flux of water from each soil layer (kg m-2 s-1).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content as a fraction of saturation.
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Unfrozen soil moisture content as a fraction of saturation.
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Subsurface temperature in each layer (K).
  w_flux_soilt(land_pts,nsoilt,0:sm_levels)
    ! Downward water flux at bottom of each soil layer (kg m-2 s-1).

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INCREMENT_SOIL_DRIVERS'

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER :: j,l  ! Indices.

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_zero_vars ) THEN
  !-------------------------------------------------------------------------
  ! Reset accumulations.
  ! This is done even if l_driver_ave=F so that we can always add the next
  ! value.
  !-------------------------------------------------------------------------
  deposition_n_driver(:) = 0.0
  qbase_l_driver(:,:,:)  = 0.0
  sthf_driver(:,:,:)     = 0.0
  sthu_driver(:,:,:)     = 0.0
  tsoil_driver(:,:,:)    = 0.0
  wflux_driver(:,:,:)    = 0.0

ELSE
  !------------------------------------------------------------------------
  ! Add to accumulations.
  ! For instantaneous values this is only required when the soil_model is
  ! to be run, except that deposition is always averaged (so as to conserve
  ! mass across the atmos-land interface).
  !------------------------------------------------------------------------
  DO j = 1,soil_pts
    l = soil_index(j)
    deposition_n_driver(l) = deposition_n_driver(l) + deposition_n_gb(l)
    IF ( l_driver_ave .OR. l_call ) THEN
      sthf_driver(l,:,:)  = sthf_driver(l,:,:)  + sthf_soilt(l,:,:)
      sthu_driver(l,:,:)  = sthu_driver(l,:,:)  + sthu_soilt(l,:,:)
      tsoil_driver(l,:,:) = tsoil_driver(l,:,:) + t_soil_soilt(l,:,:)
      wflux_driver(l,:,:) = wflux_driver(l,:,:) + w_flux_soilt(l,:,:)
      ! We don't use the deepest layer of qbase_l_soilt.
      IF ( l_top ) THEN
        qbase_l_driver(l,:,:) = qbase_l_driver(l,:,:) +                       &
                                qbase_l_soilt(l,:,1:sm_levels)
      END IF
    END IF
  END DO
  !-------------------------------------------------------------------------
  ! Calculate time averages if the soil_model is to be run.
  !-------------------------------------------------------------------------
  IF ( l_call ) THEN
    deposition_n_driver(:) = deposition_n_driver(:) / REAL( nstep )
    IF ( l_driver_ave ) THEN
      sthf_driver(:,:,:)  = sthf_driver(:,:,:)  / REAL( nstep )
      sthu_driver(:,:,:)  = sthu_driver(:,:,:)  / REAL( nstep )
      tsoil_driver(:,:,:) = tsoil_driver(:,:,:) / REAL( nstep )
      wflux_driver(:,:,:) = wflux_driver(:,:,:) / REAL( nstep )
      IF ( l_top ) THEN
        qbase_l_driver(:,:,:) = qbase_l_driver(:,:,:) / REAL( nstep )
      END IF
    END IF
  END IF

END IF  !  l_zero_vars

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE increment_soil_drivers

!#############################################################################

END MODULE soil_biogeochem_control_mod
#endif
