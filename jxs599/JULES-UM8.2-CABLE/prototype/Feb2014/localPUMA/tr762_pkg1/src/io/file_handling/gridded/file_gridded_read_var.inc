#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_gridded_read_var(file, var_id, extract_subgrid, subgrid) RESULT(cube)

  USE grid_utils_mod, ONLY : subgrid_info, OPERATOR(.ne.), subgrid_extract

  USE file_mod, ONLY : file_read_var

  USE data_cube_mod, ONLY : data_cube, cube_create

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads data from the given variable in the given file as a data cube
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_gridded), INTENT(INOUT) :: file
                                 ! The file to read data from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  LOGICAL, INTENT(IN) :: extract_subgrid
                                 ! T - extract a subgrid to return from the
                                 !     full grid of the file
                                 ! F - return the full grid
  TYPE(subgrid_info), OPTIONAL, INTENT(IN) :: subgrid  ! The subgrid to extract


! Return type
  TYPE(data_cube) :: cube  ! The data cube read from file


! Work variables
  TYPE(var_gridded) :: var  ! The variable we are reading from

  REAL, ALLOCATABLE :: values_3d(:,:,:)  ! 3D version of the cube values to
                                         ! simplify calculations using subgrids

  REAL, ALLOCATABLE :: data(:,:,:)  ! Used to hold the full grid of data if
                                    ! extracting a subgrid in memory
                                    ! Uses a combined levels dimension to
                                    ! simplify calculations

  INTEGER :: i, j, k, x, y  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check if a subgrid has been requested and that the parent grid is the same
! as the grid the file is on
!-----------------------------------------------------------------------------
  IF ( extract_subgrid ) THEN
    IF ( .NOT. PRESENT(subgrid) )                                             &
      CALL log_fatal("file_gridded_read_var",                                 &
                     "Subgrid extraction has been requested but no " //       &
                     "subgrid has been given")

! Check that the subgrid's parent grid matches the file's grid
    IF ( subgrid%parent /= file%grid )                                        &
      CALL log_fatal("file_gridded_read_var",                                 &
                     "Subgrid is inconsistent with grid in file")
  END IF  ! extract_subgrid

!-----------------------------------------------------------------------------
! Extract information on the variable being read into a local variable for
! convenience
!-----------------------------------------------------------------------------
  var = file%vars(var_id)

!-----------------------------------------------------------------------------
! Create a cube of the correct size
!-----------------------------------------------------------------------------
  IF ( extract_subgrid ) THEN
    cube = cube_create((/ subgrid%nx, subgrid%ny, var%lev_sizes /))
  ELSE
    cube = cube_create((/ file%grid%nx, file%grid%ny, var%lev_sizes /))
  END IF

!-----------------------------------------------------------------------------
! Actually read the requested data
!
! The decision process for what to actually read from file is as follows:
!   * If no subgrid has been specified, read everything (obviously!)
!
!   * If a subgrid has been specified using a region, then read only that
!     region from file (easy to specify with start and count)
!
!   * If a subgrid has been specified using points, then one of two things
!     happens:
!       1. If only a small proportion of the full grid is to be extracted, we
!          read the values from file one at a time
!       2. If a significant proportion of the full grid is to extracted, we
!          read all the values and perform the extraction in memory
!
! This process provides the best trade-off between memory usage and time
! spent in I/O
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! If we are not extracting a subgrid, just read all the data into the cube
! and return
!-----------------------------------------------------------------------------
  IF ( .NOT. extract_subgrid ) THEN
    IF ( file%grid%is_1d ) THEN
      CALL file_read_var(                                                     &
        file%fh, var%id, cube%values,                                         &
! We want to read the whole of the grid and levels dimensions, so set start
! and count accordingly
! For a 1D grid, we only need to consider the x axis of the grid
        (/ 1,            (1, i = 1,SIZE(var%lev_sizes)) /),                   &
        (/ file%grid%nx, var%lev_sizes /)                                     &
      )
    ELSE
      CALL file_read_var(                                                     &
        file%fh, var%id, cube%values,                                         &
! We want to read the whole of the grid and levels dimensions, so set start
! and count accordingly
! For a 2D grid, we only need to consider both axes of the grid
        (/ 1,            1,            (1, i = 1,SIZE(var%lev_sizes)) /),     &
        (/ file%grid%nx, file%grid%ny, var%lev_sizes /)                       &
      )
    END IF

    RETURN
  END IF

!-----------------------------------------------------------------------------
! If we get to here, we know we have to perform a subgrid extraction
!-----------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED(subgrid%points) ) THEN

! If the subgrid is specified using a region, just extract the region from file

    IF ( file%grid%is_1d ) THEN
      CALL file_read_var(                                                     &
        file%fh, var%id, cube%values,                                         &
! We want to read the whole of the levels dimensions, but only the correct
! chunk of the grid dimensions, so set start and count accordingly
! For a 1D grid, we only need to consider the 'x axis' of the grid
        (/ subgrid%x_start, (1, i = 1,SIZE(var%lev_sizes)) /),                &
        (/ subgrid%nx,      var%lev_sizes /)                                  &
      )
    ELSE
      CALL file_read_var(                                                     &
        file%fh, var%id, cube%values,                                         &
! We want to read the whole of the levels dimensions, but only the correct
! chunk of the grid dimensions, so set start and count accordingly
! For a 1D grid, we only need to consider the 'x axis' of the grid
        (/ subgrid%x_start, subgrid%y_start, (1, i = 1,SIZE(var%lev_sizes)) /),&
        (/ subgrid%nx,      subgrid%ny,      var%lev_sizes /)                 &
      )
    END IF

  ELSE

! If the subgrid is specified using points, then check whether we will do the
! extraction directly from file or in memory
! For now, we only extract directly from file if reading <= 1% of points
! Testing has shown this to be around the point where extracting in memory becomes
! faster for large files

! For this section, it is easier to deal with the data as if it has two grid
! dimensions and a combined z dimension that represents all the levels combined
    ALLOCATE(values_3d(subgrid%nx, subgrid%ny, PRODUCT(var%lev_sizes)))

    IF ( ( 100 * subgrid%nx * subgrid%ny ) > ( file%grid%nx * file%grid%ny ) ) THEN
! If the extraction is to take place in memory, then allocate and read the
! full grid of data
! Again, we use a combined z dimension here
! See above for explaination of start and count
      ALLOCATE(data(file%grid%nx, file%grid%ny, PRODUCT(var%lev_sizes)))
      IF ( file%grid%is_1d ) THEN
        CALL file_read_var(                                                   &
          file%fh, var%id, data,                                              &
          (/ 1,            (1, i = 1,SIZE(var%lev_sizes)) /),                 &
          (/ file%grid%nx, var%lev_sizes /)                                   &
        )
      ELSE
        CALL file_read_var(                                                   &
          file%fh, var%id, data,                                              &
          (/ 1,            1,            (1, i = 1,SIZE(var%lev_sizes)) /),   &
          (/ file%grid%nx, file%grid%ny, var%lev_sizes /)                     &
        )
      END IF

! Extract the points from the data - we use a utility routine to do this
! one level at a time
      DO k = 1,SIZE(data, 3)
        values_3d(:,:,k) = subgrid_extract(subgrid, data(:,:,k))
      END DO

      DEALLOCATE(data)

    ELSE

! Read the points one at a time directly from file
      DO j = 1,subgrid%ny
        DO i = 1,subgrid%nx
! Translate the index in subgrid%points into x and y coordinates in the file grid
          y = (subgrid%points(i,j) - 1) / file%grid%nx + 1
          x = subgrid%points(i,j) - (y - 1) * file%grid%nx

          IF ( file%grid%is_1d ) THEN
            CALL file_read_var(file%fh, var%id, values_3d(i,j,:),             &
! For a 1D grid, we only need to specify the x index
                               (/ x, (1, i = 1,SIZE(var%lev_sizes)) /),       &
                               (/ 1, var%lev_sizes /))
          ELSE
            CALL file_read_var(file%fh, var%id, values_3d(i,j,:),             &
! For a 2D grid, we specify the x and y indices
                               (/ x, y, (1, i = 1,SIZE(var%lev_sizes)) /),    &
                               (/ 1, 1, var%lev_sizes /))
          END IF
        END DO
      END DO
    END IF

! Convert the 3d values back into a 1d array
    cube%values(:) = RESHAPE(values_3d, (/ SIZE(values_3d) /))
    DEALLOCATE(values_3d)

  END IF  ! subgrid specified using region or points

  RETURN

END FUNCTION file_gridded_read_var
#endif
