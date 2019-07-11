#if !defined(UM_JULES)
! Module containing variables required for orographic roughness enhancement.

  MODULE orog

  REAL, DIMENSION(:),   ALLOCATABLE :: SIL_OROG_LAND
!                                    ! Silhouette area of unresolved orography
!                                    !   per unit horizontal area on land points only
  REAL, DIMENSION(:),   ALLOCATABLE :: HO2R2_OROG
!                                    ! Standard Deviation of orography
!                                    !   equivalent to peak to trough height of unresolved orography
  REAL, DIMENSION(:,:), ALLOCATABLE :: H_BLEND_OROG
!                                    ! Blending height used as part of effective roughness scheme
  REAL, DIMENSION(:,:), ALLOCATABLE :: Z0M_EFF
!                                    ! Effective grid-box roughness length for momentum

  END MODULE orog
#endif
