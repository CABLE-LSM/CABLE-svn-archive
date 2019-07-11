#if !defined(UM_JULES)

MODULE screen

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module containing screen T&Q and 10m U&V variables
!-----------------------------------------------------------------------------

  REAL, DIMENSION(:,:), ALLOCATABLE :: Q1P5M
!                                    ! Q at 1.5 m (kg water / kg air)
  REAL, DIMENSION(:,:), ALLOCATABLE :: Q1P5M_TILE
!                                    ! Q1P5M over land tiles
  REAL, DIMENSION(:,:), ALLOCATABLE :: T1P5M
!                                    ! T at 1.5 m (K)
  REAL, DIMENSION(:,:), ALLOCATABLE :: T1P5M_TILE
!                                    ! T1P5M over land tiles
  REAL, DIMENSION(:,:), ALLOCATABLE :: U10M
!                                    ! U at 10 m (m per s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: V10M
!                                    ! V at 10 m (m per s)

END MODULE screen
#endif
