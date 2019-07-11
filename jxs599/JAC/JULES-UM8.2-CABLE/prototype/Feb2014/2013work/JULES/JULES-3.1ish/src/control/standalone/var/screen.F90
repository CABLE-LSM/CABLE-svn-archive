#if !defined(UM_JULES)
!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/var/screen.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

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
