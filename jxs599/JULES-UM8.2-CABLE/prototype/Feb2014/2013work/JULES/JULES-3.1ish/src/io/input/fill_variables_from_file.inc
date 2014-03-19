#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/input/fill_variables_from_file.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE fill_variables_from_file(file_name, identifiers, sdf_names)

  USE io_constants, ONLY : MODE_READ, MAX_SDF_NAME_LEN, MAX_DIM_FILE

  USE dictionary_mod, ONLY : dict, dict_create, dict_get, dict_set, dict_free

  USE model_interface_mod, ONLY : get_var_id, get_var_levs_dim, populate_var

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
  INTEGER :: nvars  ! The number of variables requested
  INTEGER :: var_ids(SIZE(identifiers))
                                 ! The variable id from model_interface_mod
                                 ! for each requested variable
  INTEGER :: var_file_ids(SIZE(identifiers))
                                 ! The variable id in the file for each
                                 ! requested variable
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_levs_dims(SIZE(identifiers))
                                 ! The name of the levels dimension for each
                                 ! requested variable
  INTEGER :: var_nlevs(SIZE(identifiers))
                                 ! The size of the levels dimension for each
                                 ! requested variable

  TYPE(dict) :: file_dims  ! Dictionary containing the dimensions to define
                           ! in the file
                           ! Maps dim_name => dim_size

  TYPE(dict) :: file_dim_ids  ! Dictionary containing the dimension ids in file
                              ! Maps dim_name => dim_id

  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name  ! Temporary variables used to
  INTEGER :: dim_size                          ! populate dimension dictionaries
  INTEGER :: dim_id
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name_ignore  ! Used as a dummy argument

  INTEGER :: levels_dim_id  ! The id in the gridded file of the vertical levels
                            ! dimension for the variable being processed
                            ! - used when defining variables

  TYPE(file_gridded) :: file  ! The file object that we are reading from

  INTEGER :: data_nx, data_ny  ! The sizes of the x and y dimensions for data
  REAL, ALLOCATABLE :: data(:,:,:)
      ! Variable to temporarily store data on the model grid

  INTEGER :: i  ! Loop counters


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  file_dims    = dict_create(MAX_DIM_FILE, INT(1))
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
! Process
!-----------------------------------------------------------------------------
! Map the string identifiers to their integer ids
  DO i = 1,nvars
    var_ids(i) = get_var_id(identifiers(i))
  END DO

! Make a first pass over the requested variables to work out what dimensions
! we expect to have in the file
  DO i = 1,nvars
    CALL get_var_levs_dim(var_ids(i), dim_name, dim_name_ignore, dim_size)

! Store the levels dim name and size for this variable for use later
    var_levs_dims(i) = dim_name
    var_nlevs(i) = dim_size

! For a variable with no levels dimension, the returned dim_name will be the
! empty string - we don't want to add this to the dictionary
    IF ( LEN_TRIM(dim_name) < 1 ) CYCLE

! Otherwise store it in the dictionary
! Duplicate keys will overwrite, so we don't have to worry about duplication
    CALL dict_set(file_dims, dim_name, dim_size)
  END DO

!-----------------------------------------------------------------------------
! Open the given file and define its grid and vertical dimensions
!-----------------------------------------------------------------------------
  file = file_gridded_open(file_name, MODE_READ)

! Define the grid
  CALL file_gridded_def_grid(                                                 &
    file, grid_is_1d, grid_dim_name, nx, x_dim_name, nx, y_dim_name, ny       &
  )

! Define the required vertical levels dimensions
  DO i = 1,file_dims%length
! The dimension names are the dictionary keys, and the sizes are the values
    dim_name = file_dims%keys(i)
    CALL dict_get(file_dims, dim_name, dim_size)

! Define the dimension
    dim_id = file_gridded_def_dim(file, dim_name, dim_size)

! Store the id
    CALL dict_set(file_dim_ids, dim_name, dim_id)
  END DO

!-----------------------------------------------------------------------------
! Define the variables, gathering their ids as we go for future use
!-----------------------------------------------------------------------------
  DO i = 1,nvars
! Get the dimension id to use for the levels dim for this variable
    IF ( LEN_TRIM(var_levs_dims(i)) < 1 ) THEN
! If an empty string was returned as the levels dim name, then the variable
! has no levels dimension
! We exploit the fact that giving -1 as the dimension id is equivalent to not
! giving a dimension id
      levels_dim_id = -1
    ELSE
! Otherwise retrieve the dimension id from the dictionary
      CALL dict_get(file_dim_ids, var_levs_dims(i), levels_dim_id)
    END IF

! Define the variable and store its id
    var_file_ids(i) = file_gridded_def_var(                                   &
      file, sdf_names(i), levels_dim_id, .FALSE.                              &
    )
  END DO

! We have finished defining things!
  CALL file_gridded_enddef(file)

!-----------------------------------------------------------------------------
! Read each variable and fill the appropriate model variable
!-----------------------------------------------------------------------------
  data_nx = nx
  data_ny = ny
  IF ( use_subgrid ) THEN
    data_nx = subgrid_npoints
    data_ny = 1
  END IF

  ALLOCATE(data(data_nx, data_ny, MAXVAL(var_nlevs)))

  DO i = 1,nvars
    CALL file_gridded_read_var(                                               &
      file, var_file_ids(i), data(:,:,1:var_nlevs(i)),                        &
      use_subgrid, points, x_coords, y_coords                                 &
    )

    CALL populate_var(var_ids(i), data(:,:,1:var_nlevs(i)))
  END DO

! Do some cleanup
  DEALLOCATE(data)

  CALL file_gridded_close(file)

  CALL dict_free(file_dims)
  CALL dict_free(file_dim_ids)

  RETURN

END SUBROUTINE fill_variables_from_file
#endif
