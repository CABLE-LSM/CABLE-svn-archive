#if !defined(UM_JULES)
! Module containing variables required on alternative grids

MODULE u_v_grid

IMPLICIT NONE

REAL, DIMENSION(:,:), ALLOCATABLE :: u_0_p_ij
  ! W'ly component of surface current (m/s). P grid
REAL, DIMENSION(:,:), ALLOCATABLE :: v_0_p_ij
  ! S'ly component of surface current (m/s). P grid
REAL, DIMENSION(:,:), ALLOCATABLE :: u_1_p_ij
  ! U_1 on P-grid
REAL, DIMENSION(:,:), ALLOCATABLE :: v_1_p_ij
  ! V_1 on P-grid
REAL, DIMENSION(:,:), ALLOCATABLE :: dtrdz_charney_grid_1_ij
  ! -g.dt/dp for model layers

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='U_V_GRID'

CONTAINS
  
SUBROUTINE u_v_grid_alloc(t_i_length,t_j_length)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: t_i_length,t_j_length

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='U_V_GRID_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( u_0_p_ij(t_i_length,t_j_length))
ALLOCATE( v_0_p_ij(t_i_length,t_j_length))
ALLOCATE( u_1_p_ij(t_i_length,t_j_length))
ALLOCATE( v_1_p_ij(t_i_length,t_j_length))
ALLOCATE( dtrdz_charney_grid_1_ij(t_i_length,t_j_length))

u_0_p_ij(:,:)    = 0.0
v_0_p_ij(:,:)    = 0.0
u_1_p_ij(:,:)    = 0.0
v_1_p_ij(:,:)    = 0.0
dtrdz_charney_grid_1_ij(:,:) = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE u_v_grid_alloc

END MODULE u_v_grid
#endif
