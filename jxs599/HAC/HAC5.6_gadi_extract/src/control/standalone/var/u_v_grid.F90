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

END MODULE u_v_grid
#endif
