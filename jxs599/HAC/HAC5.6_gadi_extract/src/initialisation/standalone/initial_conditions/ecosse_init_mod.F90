#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE ecosse_init_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Description:
!   This module contains code that is used to initialise the ECOSSE model, in
!   particular to ensure consistency between variables.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in BIOGEOCHEMISTRY
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.

PRIVATE
PUBLIC ecosse_init

CONTAINS

!#############################################################################
!#############################################################################

SUBROUTINE ecosse_init

USE ancil_info, ONLY:                                                         &
  ! imported scalars
  land_pts, nsoilt, nz_soilc=>dim_cslayer,                                    &
  ! imported arrays
  frac_surft

USE ecosse_prepare_mod, ONLY:                                                 &
  ! imported procedures
  stable_n_c

USE ecosse_utils_mod, ONLY:                                                   &
  ! imported procedures
  adjust_soil, get_residual_n

USE jules_soil_ecosse_mod, ONLY:                                              &
  ! imported scalars
  l_soil_N

USE p_s_parms, ONLY:                                                          &
  ! imported arrays
  soil_pH_soilt

USE prognostics, ONLY:                                                        &
  ! imported arrays
  cs_pool_soilt, ns_pool_gb

USE soil_ecosse_vars_mod, ONLY:                                               &
  ! imported arrays
  n_amm_soilt, n_nit_soilt

USE veg_soil_index_mod, ONLY:                                                 &
  ! imported procedures
  get_veg_soil_index

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
LOGICAL, PARAMETER ::                                                         &
  l_init_true = .TRUE.
    ! For use as an argument to subroutine adjust_soil, to indicate to that
    ! routine that it is being called during initialisation.

!-----------------------------------------------------------------------------
! Local scalar variables
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  vs_pts,                                                                     &
    ! The number of points with veg and/or soil.
  i, l,                                                                       &
    ! Indices.
  s ! Soil tile number.

!-----------------------------------------------------------------------------
! Local array variables
!-----------------------------------------------------------------------------
INTEGER :: vs_index(land_pts)      !  indices of veg/soil points

REAL(KIND=real_jlslsm) ::                                                     &
  biohum_nc(land_pts,nz_soilc),                                               &
    ! Stable N:C ratio of biomass and humus pools.
  co2_dummy(land_pts),                                                        &
    ! Local or "dummy" variable so we can provide an argument. Notionally this
    ! is the C in the CO2 flux from soil to atmosphere (kg m-2 s-1).
  frac_vs(land_pts),                                                          &
    ! Fraction of gridbox covered by veg or soil.
  residual_n(land_pts,nz_soilc)
    ! Minimum-allowed (residual) inorganic N amount (kg m-2).

!-----------------------------------------------------------------------------
!end of header

!-----------------------------------------------------------------------------
! Get index for points with soil and/or vegetation.
!-----------------------------------------------------------------------------
CALL get_veg_soil_index( land_pts, frac_surft, vs_pts,                        &
                         vs_index, frac_vs )

!-----------------------------------------------------------------------------
! Loop over soil tiles.
! Note that many of these variables don't yet support soil tiling, so this
! loop is rather academic!
!-----------------------------------------------------------------------------
DO s = 1,nsoilt

  IF ( l_soil_N ) THEN

    !-------------------------------------------------------------------------
    ! Calculate stable N:C of biomass and humus pools.
    !-------------------------------------------------------------------------
    CALL stable_n_c( land_pts, vs_pts, vs_index, soil_pH_soilt(:,s,:),        &
                     biohum_nc )

    !------------------------------------------------------------------------
    ! Ensure that inorganic N pools are not smaller than the minimum allowed.
    ! This should only really be done during a "hard" (i.e. comprehensive)
    ! reconfiguration, but for now we will always do it.
    !------------------------------------------------------------------------
    ! Calculate the minimum-allowed amount of inorganic N.
    CALL get_residual_n( land_pts, vs_pts, vs_index, residual_n )

    DO i = 1,vs_pts
      l = vs_index(i)
      n_amm_soilt(l,s,:) = MAX( n_amm_soilt(l,s,:), residual_n(l,:) )
      n_nit_soilt(l,s,:) = MAX( n_nit_soilt(l,s,:), residual_n(l,:) )
    END DO

  END IF  !  l_soil_N

  ! Enforce minimum-allowed pool sizes (which also ensures C:N of biomass and
  ! humus pools). Again this should only really be done during a "hard"
  ! (i.e. comprehensive) reconfiguration but we will do it now to avoid
  ! issues such as zero-sized pools.
  co2_dummy(:) = 0.0
  CALL adjust_soil( land_pts, vs_pts, l_init_true, vs_index, biohum_nc,       &
                    cs_pool_soilt(:,s,:,1), cs_pool_soilt(:,s,:,2),           &
                    cs_pool_soilt(:,s,:,3), cs_pool_soilt(:,s,:,4),           &
                    ns_pool_gb(:,:,1), ns_pool_gb(:,:,2),                     &
                    ns_pool_gb(:,:,3), ns_pool_gb(:,:,4), co2_dummy )

END DO  !  soil tiles

END SUBROUTINE ecosse_init

!#############################################################################
!#############################################################################

END MODULE ecosse_init_mod
#endif
