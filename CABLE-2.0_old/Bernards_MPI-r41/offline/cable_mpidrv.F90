! cable_mpidrv.f90
!
! MPI offline driver for CABLE land surface scheme, April 2011.
! Bare bones scaffolding...
! Maciej Golebiewski, CSIRO Advanced Scientific Computing.
!
! Please send bug reports to Bernard.Pak@csiro.au
!
PROGRAM mpi_driver

  USE mpi

  USE cable_mpicommon
  USE cable_mpimaster
  USE cable_mpiworker

  IMPLICIT NONE

  INTEGER :: comm, np, rank, ierr

  CALL MPI_Init (ierr)
  CALL MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  CALL MPI_Comm_size (comm, np, ierr)

  IF (np < 2) THEN
     WRITE (*,*) 'This program needs at least 2 processes to run!'
     CALL MPI_Abort (comm, 0, ierr)
  END IF

  CALL MPI_Comm_rank (comm, rank, ierr)

  IF (rank == 0) THEN
          CALL mpidrv_master (comm)
  ELSE
          CALL mpidrv_worker (comm)
  END IF

  CALL MPI_Finalize (ierr)

END PROGRAM mpi_driver

