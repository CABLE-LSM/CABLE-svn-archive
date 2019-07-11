#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE fill_variables_from_file(file_name, identifiers, sdf_names)

  USE io_constants, ONLY : MODE_READ, MAX_SDF_NAME_LEN, MAX_DIM_FILE, MAX_DIM_VAR

  USE dictionary_mod, ONLY : dict, dict_create, dict_get, dict_set,           &
                                   dict_has_key, dict_free

  USE data_cube_mod, ONLY : data_cube, cube_free

  USE model_interface_mod, ONLY : get_var_id, get_var_levs_dims, populate_var

  USE file_gridded_mod, ONLY : file_gridded, file_gridded_open,               &
                               file_gridded_def_grid, file_gridded_def_dim,   &
                               file_gridded_def_var, file_gridded_enddef,     &
                               file_gridded_read_var, file_gridded_close


  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Fills the model variables specified by the given identifiers using
!   variables from the given file specified by sdf_names
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Argument types
!-----------------------------------------------------------------------------
  CHARACTER(len=*) :: file_name
      ! The name of the file to open
  CHARACTER(len=*) :: identifiers(:)
      ! The model identifiers of the variables to fill
  CHARACTER(len=*) :: sdf_names(:)
      ! The name of the variable in the file for each identifier


!-----------------------------------------------------------------------------
! Work variables
!-----------------------------------------------------------------------------
! Variables that need to be carried over from the definition phase to the
! reading phase
  INTEGER :: nvars  ! The number of variables requested
  INTEGER :: var_ids(SIZE(identifiers))
                    ! The variable ids as used by model_interface_mod
  INTEGER :: var_file_ids(SIZE(identifiers))
                    ! The ids of the variables in the opened file
  TYPE(file_gridded) :: file  ! The file object that we are reading from


! Variables used during definition loop
  INTEGER :: ndims  ! The number of levels dimensions the current variable has
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_VAR)
                    ! The names of the dimensions for the current variable
  INTEGER :: dim_sizes(MAX_DIM_VAR)
                    ! The sizes of the dimensions for the current variable
  INTEGER :: dim_ids(MAX_DIM_VAR)
                    ! The ids of the dimensions for the current variable
  TYPE(dict) :: file_dim_ids  ! Dictionary containing the dimension ids in file
                              ! Maps dim_name => dim_id

! Workspace cube that can be deallocated to avoid memory leaks
  TYPE(data_cube) :: data


  INTEGER :: i, j  ! Loop counters


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  file_dim_ids = dict_create(MAX_DIM_FILE, INT(1))

!-----------------------------------------------------------------------------
! Check arguments
!-----------------------------------------------------------------------------
  nvars = SIZE(identifiers)
  IF ( nvars /= SIZE(sdf_names) )                                             &
    CALL log_fatal("fill_variables_from_file",                                &
                   "identifiers and sdf_names must have the same number " //  &
                   "of elements")

!-----------------------------------------------------------------------------
! Map the string identifiers to their integer ids
!-----------------------------------------------------------------------------
  DO i = 1,nvars
    var_ids(i) = get_var_id(identifiers(i))
  END DO

!-----------------------------------------------------------------------------
! Open the given file and define its grid
!-----------------------------------------------------------------------------
  file = file_gridded_open(file_name, MODE_READ)

! Define the grid
  CALL file_gridded_def_grid(file, grid)

!-----------------------------------------------------------------------------
! Define the required dimensions and variables
!-----------------------------------------------------------------------------
  DO i = 1,nvars
! Get the levels dims used by this variable - we only care about the names
! used in input files
    CALL get_var_levs_dims(var_ids(i), ndims=ndims,                           &
                           dim_names_in=dim_names, dim_sizes=dim_sizes)

    DO j = 1,ndims
! If it has not yet been defined, define the dimension, storing its id
      IF ( .NOT. dict_has_key(file_dim_ids, dim_names(j)) )                   &
        CALL dict_set(                                                        &
          file_dim_ids, dim_names(j),                                         &
          file_gridded_def_dim(file, dim_names(j), dim_sizes(j))              &
        )

! Get the dimension id from the dict and add it to the list for this variable
      CALL dict_get(file_dim_ids, dim_names(j), dim_ids(j))
    END DO

! Create the variable and store its id
    var_file_ids(i) = file_gridded_def_var(                                   &
      file, sdf_names(i), dim_ids(1:ndims), .FALSE.                           &
    )
  END DO

!-----------------------------------------------------------------------------
! We have finished defining things!
!-----------------------------------------------------------------------------
  CALL dict_free(file_dim_ids)
  CALL file_gridded_enddef(file)

!-----------------------------------------------------------------------------
! Read each variable and fill the appropriate model variable
!-----------------------------------------------------------------------------
  DO i = 1,nvars
! Read the data cube from file and pass it to populate_var
    data = file_gridded_read_var(file, var_file_ids(i), use_subgrid, subgrid)
    CALL populate_var(var_ids(i), data)
    CALL cube_free(data)
  END DO

  CALL file_gridded_close(file)

  RETURN

END SUBROUTINE fill_variables_from_file
#endif
