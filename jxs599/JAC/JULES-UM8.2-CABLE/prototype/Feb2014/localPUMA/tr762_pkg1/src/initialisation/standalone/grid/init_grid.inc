#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_grid(nml_dir)

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model grids
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists

! Work variables
  INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------


! Open the grid namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'model_grid.nml'),        &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_grid",                                               &
                   "Error opening namelist file model_grid.nml " //           &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Defer to specialised routines to initialise the different aspect of the grid
  CALL init_input_grid()
  CALL init_latlon()
  CALL init_land_frac()
  CALL init_model_grid()
  CALL init_surf_hgt()

  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_grid",                                               &
                   "Error closing namelist file model_grid.nml " //           &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  RETURN

END SUBROUTINE init_grid
#endif
