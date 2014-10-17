#if !defined(UM_JULES)
! Module containing variables required on alternative grids

  MODULE u_v_grid

  REAL, DIMENSION(:,:), ALLOCATABLE :: U_0_P
!                                    ! W'ly component of surface current (m/s). P grid
  REAL, DIMENSION(:,:), ALLOCATABLE :: V_0_P
!                                    ! S'ly component of surface current (m/s). P grid
  REAL, DIMENSION(:,:), ALLOCATABLE :: U_1_P
!                                    ! U_1 on P-grid
  REAL, DIMENSION(:,:), ALLOCATABLE :: V_1_P
!                                    ! V_1 on P-grid
  REAL, DIMENSION(:,:), ALLOCATABLE :: DTRDZ_CHARNEY_GRID_1
!                                    ! -g.dt/dp for model layers

  END MODULE u_v_grid
#endif
