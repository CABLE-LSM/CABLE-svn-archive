#if !defined(UM_JULES)
! Module containing variables required for orographic roughness enhancement.

MODULE orog

IMPLICIT NONE

REAL, DIMENSION(:),   ALLOCATABLE :: sil_orog_land_gb
              ! Silhouette area of unresolved orography
              !   per unit horizontal area on land points only
REAL, DIMENSION(:),   ALLOCATABLE :: ho2r2_orog_gb
              ! Standard Deviation of orography
              ! equivalent to peak to trough height of unresolved orography
REAL, DIMENSION(:,:), ALLOCATABLE :: h_blend_orog_ij
              ! Blending height used as part of effective roughness scheme
REAL, DIMENSION(:,:), ALLOCATABLE :: z0m_eff_ij
              ! Effective grid-box roughness length for momentum

END MODULE orog
#endif
