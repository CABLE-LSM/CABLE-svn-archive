#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE init_parms_mod
CONTAINS
SUBROUTINE init_parms()

  !Use in relevant subroutines
USE sparm_mod, ONLY: sparm
USE infiltration_rate_mod, ONLY: infiltration_rate

USE theta_field_sizes, ONLY: t_i_length, t_j_length

USE jules_surface_mod, ONLY: l_aggregate

USE jules_sea_seaice_mod, ONLY: nice

USE ancil_info, ONLY: ssi_index, fssi_ij, ice_fract_ij,                       &
                       ice_fract_ncat_sicat,                                  &
                       sea_frac, sice_frac, sice_frac_ncat, sea_index,        &
                       sice_index, sice_frac, sea_frac, sice_pts_ncat,        &
                       sice_index_ncat, land_pts, nsurft, sea_pts,            &
                       sice_pts, ssi_pts, surft_pts, surft_index,             &
                       frac_surft, land_index

USE coastal, ONLY: tstar_sice_ij, tstar_land_ij, tstar_ssi_ij, tstar_sea_ij,  &
                   flandg, tstar_sice_sicat

USE prognostics, ONLY: tstar_surft, canht_pft, lai_pft

USE p_s_parms, ONLY: catch_surft, catch_snow_surft, infil_surft,              &
                      satcon_soilt, z0_surft, z0h_bare_surft, z0m_soil_gb

USE fluxes, ONLY: tstar_ij

USE u_v_grid, ONLY: dtrdz_charney_grid_1_ij

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
INTEGER :: i,j,l,n  ! Loop counters


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Calculate surface parameters.
!-----------------------------------------------------------------------------
CALL sparm(land_pts, nsurft, surft_pts, surft_index,                          &
           frac_surft, canht_pft, lai_pft, z0m_soil_gb,                       &
           catch_snow_surft, catch_surft, z0_surft, z0h_bare_surft)

CALL infiltration_rate(land_pts, nsurft, surft_pts, surft_index,              &
                       satcon_soilt, frac_surft, infil_surft)

!-----------------------------------------------------------------------
! Set up index for sea and sea-ice
!-----------------------------------------------------------------------
ssi_pts = 0
ssi_index(:) = 0
DO j = 1,t_j_length
  DO i = 1,t_i_length
    ssi_pts = ssi_pts + 1
    IF ( flandg(i,j) < 1.0 ) THEN
      ssi_index(ssi_pts) = (j-1) * t_i_length + i
    END IF
    fssi_ij(i,j) = 1.0 - flandg(i,j)
  END DO
END DO

!-----------------------------------------------------------------------
! Set sea ice fraction.
!-----------------------------------------------------------------------
DO j = 1,t_j_length
  DO i = 1,t_i_length
    ice_fract_ij(i,j) = 0.0
    tstar_sice_ij(i,j) = 0.0
    DO n = 1,nice
      ice_fract_ij(i,j) = ice_fract_ij(i,j) + ice_fract_ncat_sicat(i,j,n)
    END DO
    IF (ice_fract_ij(i,j) > 0.0) THEN
      DO n = 1,nice  !assuming nice=nice_use here
        tstar_sice_ij(i,j) = tstar_sice_ij(i,j)                               &
                  + ice_fract_ncat_sicat(i,j,n) * tstar_sice_sicat(i,j,n) /   &
                                                        ice_fract_ij(i,j)
      END DO
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
! Initialise sea and sea-ice indices
!-----------------------------------------------------------------------
sea_pts  = 0
sice_pts = 0
sea_index(:)  = 0
sice_index(:) = 0
sice_frac(:) = 0.0
sea_frac(:)  = 0.0
DO l = 1,ssi_pts
  j = (ssi_index(l) - 1) / t_i_length + 1
  i = ssi_index(l) - (j-1) * t_i_length
  IF ( ssi_index(l) > 0 ) THEN
    IF ( ice_fract_ij(i,j) > 0.0 ) THEN
      sice_pts = sice_pts + 1
      sice_index(sice_pts) = l
      sice_frac(l) = ice_fract_ij(i,j)
    END IF
    IF ( ice_fract_ij(i,j) < 1.0 ) THEN
      sea_pts = sea_pts + 1
      sea_index(sea_pts) = l
      sea_frac(l) = 1.0 - sice_frac(l)
    END IF
  END IF
END DO

sice_pts_ncat(:) = 0
sice_index_ncat(:,:) = 0
sice_frac_ncat(:,:) = 0.0
DO n = 1,nice
  DO l = 1,ssi_pts
    j = (ssi_index(l) - 1) / t_i_length + 1
    i = ssi_index(l) - (j-1) * t_i_length
    IF ( ssi_index(l) > 0 ) THEN
      IF ( ice_fract_ncat_sicat(i,j,n) > 0.0 ) THEN
        sice_pts_ncat(n) = sice_pts_ncat(n) + 1
        sice_index_ncat(sice_pts_ncat(n),n) = l
        sice_frac_ncat(l,n) = ice_fract_ncat_sicat(i,j,n)
      END IF
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
! Set up gridbox "prognostics".
!-----------------------------------------------------------------------

tstar_ij(:,:)      = 0.0
tstar_land_ij(:,:) = 0.0
tstar_ssi_ij(:,:)  = 0.0

DO l = 1,land_pts
  j = (land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  IF ( l_aggregate ) THEN
    tstar_land_ij(i,j) = TSTAR_surft(l,1)
  ELSE
    DO n = 1,nsurft
      tstar_land_ij(i,j) = tstar_land_ij(i,j) + frac_surft(l,n) *             &
                                                             TSTAR_surft(l,n)
    END DO
  END IF
END DO

tstar_ssi_ij(:,:) = (1.0 - ice_fract_ij(:,:)) * tstar_sea_ij(:,:)             &
               + ice_fract_ij(:,:) * tstar_sice_ij(:,:)
tstar_ij(:,:) = flandg(:,:) * tstar_land_ij(:,:)                              &
           + (1.0 - flandg(:,:)) * tstar_ssi_ij(:,:)

!-----------------------------------------------------------------------------
! Set up information on U, V and T grids (assume that att grids are the same)
!-----------------------------------------------------------------------------
dtrdz_charney_grid_1_ij(:,:) = 0.0

RETURN

END SUBROUTINE init_parms
END MODULE init_parms_mod
#endif
