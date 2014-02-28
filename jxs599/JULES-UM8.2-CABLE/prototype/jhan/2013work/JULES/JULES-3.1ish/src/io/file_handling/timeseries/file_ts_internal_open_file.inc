#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/timeseries/file_ts_internal_open_file.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_ts_internal_open_file(file, file_name)

  USE io_constants, ONLY : MODE_READ, MODE_WRITE, ATTR_GLOBAL

  USE datetime_mod, ONLY : datetime_to_string

  USE dictionary_mod, ONLY : DICT_KEY_LEN, DICT_CHAR_VAL_LEN, dict_get

  USE file_gridded_mod, ONLY : file_gridded_open, file_gridded_def_grid,      &
                               file_gridded_def_dim,                          &
                               file_gridded_def_record_dim,                   &
                               file_gridded_def_var, file_gridded_def_attr,   &
                               file_gridded_enddef, file_gridded_close

  USE file_mod, ONLY : file_def_var, file_def_attr

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   INTERNAL PROCEDURE TO file_ts_mod
!   Sets the currently open file to that specified by file_name, including
!   defining dimensions, variables etc.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file  ! The file to update
  CHARACTER(len=*), INTENT(IN) :: file_name  ! The file to open and set as
                                             ! current file


! Work variables
  INTEGER :: levels_dim_id  ! ID to use for levels dimension of when defining
                            ! variable
  INTEGER :: attr_var_id    ! Variable ID to use when defining attributes

  CHARACTER(len=DICT_KEY_LEN) :: key  ! Used when iterating over attribute
  INTEGER :: int_val                  ! dictionaries
  REAL :: real_val
  CHARACTER(len=DICT_CHAR_VAL_LEN) :: char_val

  INTEGER, ALLOCATABLE :: zero_element_array(:)  ! This array is allocated
                                                 ! to have size 0 - used to
                                                 ! indicate passing no
                                                 ! dimensions to file_def_var
                                                 ! for the time index


  INTEGER :: n_non_time_vars  ! The number of non-time-varying variables
  INTEGER :: nx, ny  ! Number of x and y points in the grid - used to
                     ! dimension array below
  REAL, ALLOCATABLE :: non_time_varying_data(:,:,:,:)
                           ! Array to hold the non-time-varying data from the
                           ! previous file until it can be copied into the
                           ! next file
                           ! Dimensions are: variable, x, y, levels
  INTEGER :: var_nlevs  ! The size of the vertical levels dimension for the
                        ! variable currently being processed
                        ! Only use when copying non-time-varying variables
                        ! from one file to another

  INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------


  IF ( file%has_open_file ) THEN
    IF ( file%mode == MODE_WRITE ) THEN
!-----------------------------------------------------------------------------
! If we are in write mode, we want to copy the values of all the variables
! not using the time dimensions into the new file
! To do that, we have to grab their values here
!-----------------------------------------------------------------------------
! Only bother doing this if we have non-time-varying variables
      n_non_time_vars = COUNT(.NOT. file%vars(1:file%nvars)%use_time)
      IF ( n_non_time_vars > 0 ) THEN
! Work out what the grid dimensions will be
        IF ( file%grid%is_1d ) THEN
          nx = file%grid%grid_npoints
          ny = 1
        ELSE
          nx = file%grid%grid_nx
          ny = file%grid%grid_ny
        END IF

! We allocate (potentially) too much data, since we only need it for this routine
        ALLOCATE(non_time_varying_data(                                       &
          n_non_time_vars, nx, ny,                                            &
          MAXVAL( (/ file%dims(1:file%ndims)%length, 1 /) )                   &
        ))

! Gather the data for the non-time varying variables
        j = 1
        DO i = 1,file%nvars
          IF ( .NOT. file%vars(i)%use_time ) THEN
! Get the number of vertical levels that the variable has so we can
! select the correct number of elements in the data array
            IF ( file%vars(i)%has_levels ) THEN
              var_nlevs = file%dims( file%vars(i)%levels_dim )%length
            ELSE
              var_nlevs = 1
            END IF
! We are not using a subgrid, so give extract_subgrid=F
            CALL file_ts_read_var(                                            &
              file, i, non_time_varying_data(j,:,:,1:var_nlevs), .FALSE.      &
            )
! Advance the non-time variable counter
            j = j + 1
          END IF
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------------
! Close the currently open file
!-----------------------------------------------------------------------------
    CALL file_gridded_close(file%open_file)
  END IF

!-----------------------------------------------------------------------------
! Open the new file
!-----------------------------------------------------------------------------
  file%open_file = file_gridded_open(file_name, file%mode)
  file%has_open_file = .TRUE.

!-----------------------------------------------------------------------------
! Define the grid on the newly opened file
!-----------------------------------------------------------------------------
! It doesn't matter that we are providing values for all the optional arguments
! Those that are not required will be ignored
  CALL file_gridded_def_grid(file%open_file, file%grid%is_1d,                 &
                             file%grid%grid_name, file%grid%grid_npoints,     &
                             file%grid%grid_x_name, file%grid%grid_nx,        &
                             file%grid%grid_y_name, file%grid%grid_ny)

!-----------------------------------------------------------------------------
! Define the vertical level dimensions and update the dimension ids to the
! ones used by the newly opened file
!-----------------------------------------------------------------------------
  DO i = 1,file%ndims
    file%dims(i)%current_id = file_gridded_def_dim(                           &
      file%open_file, file%dims(i)%name, file%dims(i)%length                  &
    )
  END DO

!-----------------------------------------------------------------------------
! Define the time dimension as a record dimension on the open file and update
! the stored id to the one used by the newly opened file
!-----------------------------------------------------------------------------
  file%time_dim%current_id = file_gridded_def_record_dim(                     &
    file%open_file, file%time_dim%name                                        &
  )

!-----------------------------------------------------------------------------
! Define the time index if we are in write mode
!-----------------------------------------------------------------------------
! Since the time index does not want grid semantics, we have to go down to the
! raw file routines
! We pass an array with 0 elements as the non-record dimensions, but indicate
! we want to use the record dim
  IF ( file%mode == MODE_WRITE ) THEN
    ALLOCATE(zero_element_array(0))
    file%time_index_var_id = file_def_var(                                    &
      file%open_file%fh, TIME_INDEX_VAR_NAME, zero_element_array, .TRUE.      &
    )
    DEALLOCATE(zero_element_array)

! Define it's long name and units
    CALL file_def_attr(                                                       &
      file%open_file%fh, file%time_index_var_id, "long_name", "Time of data"  &
    )
    CALL file_def_attr(                                                       &
      file%open_file%fh, file%time_index_var_id, "units",                     &
      "seconds since " // datetime_to_string(file%data_start)                 &
    )
  END IF

!-----------------------------------------------------------------------------
! Define the variables and update the variable ids to the ones used by the
! newly opened file
!-----------------------------------------------------------------------------
  DO i = 1,file%nvars
! Exploit the fact that passing -ve number as levels_dim_id is the same as
! not giving one
    levels_dim_id = -1
    IF ( file%vars(i)%has_levels )                                            &
      levels_dim_id = file%dims( file%vars(i)%levels_dim )%current_id

! Define the variable, indicating whether to use the record dimension
    file%vars(i)%current_id = file_gridded_def_var(                           &
      file%open_file, file%vars(i)%name, levels_dim_id,                       &
      file%vars(i)%use_time                                                   &
    )

! Define the attributes for the variable
! First the real valued attributes
    DO j = 1,file%vars(i)%attrs_real%length
      key = file%vars(i)%attrs_real%keys(j)
      CALL dict_get(file%vars(i)%attrs_real, key, real_val)

      CALL file_gridded_def_attr(                                             &
        file%open_file, file%vars(i)%current_id, key, real_val                &
      )
    END DO

! Next, integer valued attributes
    DO j = 1,file%vars(i)%attrs_int%length
      key = file%vars(i)%attrs_int%keys(j)
      CALL dict_get(file%vars(i)%attrs_int, key, int_val)

      CALL file_gridded_def_attr(                                             &
        file%open_file, file%vars(i)%current_id, key, int_val                 &
      )
    END DO

! Lastly, character valued attributes
    DO j = 1,file%vars(i)%attrs_char%length
      key = file%vars(i)%attrs_char%keys(j)
      CALL dict_get(file%vars(i)%attrs_char, key, char_val)

      CALL file_gridded_def_attr(                                             &
        file%open_file, file%vars(i)%current_id, key, char_val                &
      )
    END DO
  END DO

!-----------------------------------------------------------------------------
! Lastly, define global attributes
!-----------------------------------------------------------------------------
! First the real valued attributes
  DO i = 1,file%attrs_real%length
    key = file%attrs_real%keys(i)
    CALL dict_get(file%attrs_real, key, real_val)

    CALL file_gridded_def_attr(file%open_file, ATTR_GLOBAL, key, real_val)
  END DO

! Next, integer valued attributes
  DO i = 1,file%attrs_int%length
    key = file%attrs_int%keys(i)
    CALL dict_get(file%attrs_int, key, int_val)

    CALL file_gridded_def_attr(file%open_file, ATTR_GLOBAL, key, int_val)
  END DO

! Lastly, character valued attributes
  DO i = 1,file%attrs_char%length
    key = file%attrs_char%keys(i)
    CALL dict_get(file%attrs_char, key, char_val)

    CALL file_gridded_def_attr(file%open_file, ATTR_GLOBAL, key, char_val)
  END DO

!-----------------------------------------------------------------------------
! Take the file out of define mode, so that it is ready to read from/write to
!-----------------------------------------------------------------------------
  CALL file_gridded_enddef(file%open_file)

!-----------------------------------------------------------------------------
! Now everything is defined on the new file, we write the data that we
! gathered for non-time-varying fields before we closed the old file
!-----------------------------------------------------------------------------
  IF ( ALLOCATED(non_time_varying_data) ) THEN
! Write the data for the non-time varying variables
    j = 1
    DO i = 1,file%nvars
      IF ( .NOT. file%vars(i)%use_time ) THEN
! Get the number of vertical levels that the variable has so we can
! select the correct number of elements in the data array
        IF ( file%vars(i)%has_levels ) THEN
          var_nlevs = file%dims( file%vars(i)%levels_dim )%length
        ELSE
          var_nlevs = 1
        END IF
! We are not using a subgrid, so give write_subgrid=F
        CALL file_ts_write_var(                                               &
          file, i, non_time_varying_data(j,:,:,1:var_nlevs), .FALSE.          &
        )
! Advance the non-time variable counter
        j = j + 1
      END IF
    END DO

    DEALLOCATE(non_time_varying_data)
  END IF

  RETURN

END SUBROUTINE file_ts_internal_open_file
#endif
