#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/timeseries/file_ts_write_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_ts_write_var(file, var_id, values, write_subgrid,             &
                                                   points, x_coords, y_coords)

  USE file_gridded_mod, ONLY : file_gridded_write_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes data given as a 3d array to the given variable in the given file
!   Data has shape (x,y,levels) where dimensions can have size 1
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file
                                 ! The file to write data to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, INTENT(IN) :: values(:,:,:)
                                 ! The values to write to the file
                                 ! The dimensions are x,y,levels where
                                 ! the size of the y and/or levels
                                 ! dimensions can be 1, depending on the
                                 ! grid configuration, points requested
                                 ! or variable

  LOGICAL, INTENT(IN) :: write_subgrid
                                 ! T - write to a subgrid of the full
                                 !     grid of the file
                                 ! F - write the full grid

  INTEGER, INTENT(IN), TARGET, OPTIONAL :: points(:)
                                 ! The indices of points to select from a
                                 ! 1d grid

  INTEGER, INTENT(IN), TARGET, OPTIONAL :: x_coords(:), y_coords(:)
                                 ! The x and y coordinates of points to
                                 ! select from a 2d grid


! Work variables
  INTEGER :: actual_var_id  ! The id of the variable in the underlying file

  INTEGER, POINTER ::                                                         &
    points_local(:), x_coords_local(:), y_coords_local(:)
      ! Local version of points, x_coords and y_coords
      ! These will be pointed to points, x_coords and y_coords respectively
      ! when those arrays are provided, otherwise they will contain one
      ! out of range point - this will trigger an error if passed with
      ! extract_subgrid=T, but not if extract_subgrid=F


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Populate local arrays depending on what has been supplied
!-----------------------------------------------------------------------------
  IF ( PRESENT(points) ) THEN
    points_local => points
  ELSE
    ALLOCATE(points_local(1))
    points_local(1) = -1
  END IF

  IF ( PRESENT(x_coords) ) THEN
    x_coords_local => x_coords
  ELSE
    ALLOCATE(x_coords_local(1))
    x_coords_local(1) = -1
  END IF

  IF ( PRESENT(y_coords) ) THEN
    y_coords_local => y_coords
  ELSE
    ALLOCATE(y_coords_local(1))
    y_coords_local(1) = -1
  END IF

!-----------------------------------------------------------------------------
! Defer the actual writing to the file_gridded routine
!-----------------------------------------------------------------------------
  actual_var_id = file%vars(var_id)%current_id

  CALL file_gridded_write_var(                                                &
    file%open_file, actual_var_id, values,                                    &
    write_subgrid, points_local, x_coords_local, y_coords_local               &
  )

  IF ( .NOT. PRESENT(points) ) DEALLOCATE(points_local)
  NULLIFY(points_local)

  IF ( .NOT. PRESENT(x_coords) ) DEALLOCATE(x_coords_local)
  NULLIFY(x_coords_local)

  IF ( .NOT. PRESENT(y_coords) ) DEALLOCATE(y_coords_local)
  NULLIFY(y_coords_local)

  RETURN

END SUBROUTINE file_ts_write_var
#endif
