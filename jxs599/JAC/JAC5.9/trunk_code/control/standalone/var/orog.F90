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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OROG'

CONTAINS

SUBROUTINE orog_alloc(land_pts,t_i_length,t_j_length)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts,t_i_length,t_j_length

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OROG_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( sil_orog_land_gb(land_pts))
ALLOCATE( ho2r2_orog_gb(land_pts))
ALLOCATE( h_blend_orog_ij(t_i_length,t_j_length))
ALLOCATE( z0m_eff_ij(t_i_length,t_j_length))

sil_orog_land_gb(:)  = 0.0
ho2r2_orog_gb(:)  = 0.0
h_blend_orog_ij(:,:)  = 0.0
z0m_eff_ij(:,:)  = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE orog_alloc

END MODULE orog
#endif
