#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


LOGICAL FUNCTION is_master_task()

  IMPLICIT NONE
  
  INCLUDE "mpif.h"
  
!-----------------------------------------------------------------------------
! Description:
!   Returns .TRUE. if the current task is the master task, .FALSE. otherwise
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: task_id  ! The id of this task
  
  INTEGER :: error  ! Error indicator for MPI calls
                    ! This is ignored as (most) MPI implementations fail
                    ! rather than returning actual error codes


!-----------------------------------------------------------------------------

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, task_id, error)
  
  is_master_task = (task_id == MASTER_TASK_ID)

  RETURN

END FUNCTION is_master_task
#endif
