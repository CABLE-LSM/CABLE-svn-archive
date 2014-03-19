#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE parallel_mod

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module provides utilities used for running in parallel using MPI
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Module parameters
!-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: MASTER_TASK_ID = 0

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
! These are set by calling decompose_domain
  INTEGER :: ntasks = 0   ! The number of MPI tasks available
  INTEGER :: task_id = -1  ! The ID of this task
  INTEGER :: task_npoints = 0  ! The number of points each task will handle

  INTEGER :: mpi_type_global_col  ! MPI type for columns in the full grid
  INTEGER :: mpi_type_local_col   ! MPI type for columns in the task local grid

! Used in scatter/gather calls
  INTEGER, ALLOCATABLE :: counts(:)  ! The number of blocks sent to each
                                     ! task. This is 1 for every task
  INTEGER, ALLOCATABLE :: offsets(:)  ! The offset, in numbers of blocks,
                                      ! for each task's block in the global
                                      ! array

!-----------------------------------------------------------------------------
! Visibilities
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC MASTER_TASK_ID, is_master_task, decompose_domain,                    &
                         scatter_land_field, gather_land_field


CONTAINS


#include "is_master_task.inc"
#include "decompose_domain.inc"
#include "scatter_land_field.inc"
#include "gather_land_field.inc"

END MODULE parallel_mod
#endif
