#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/output/internal_define_vars.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE internal_define_vars(file, var_ids, var_names, use_time,           &
                                                      var_file_ids, var_nlevs)

  USE io_constants, ONLY : MAX_DIM_FILE

  USE dictionary_mod, ONLY : DICT_KEY_LEN, DICT_CHAR_VAL_LEN, dict,           &
                             dict_create, dict_get, dict_set, dict_free

  USE model_interface_mod, ONLY : get_var_levs_dim, get_var_attrs


  USE file_ts_mod, ONLY : file_ts_def_grid, file_ts_def_dim,                  &
                          file_ts_def_time_dim, file_ts_def_var,              &
                          file_ts_def_attr, file_ts_enddef

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   INTERNAL PROCEDURE TO output_mod
!   Defines the variables identified by identifers on the given file_ts object
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file  ! The file to define variables on
  INTEGER, INTENT(IN) :: var_ids(:)
                                ! The variable ids of the model variables to
                                ! define on the file
                                ! See model_interface_mod
  CHARACTER(len=*), INTENT(IN) :: var_names(:)
                                ! The names to use for each variable in
                                ! output files
  LOGICAL, INTENT(IN) :: use_time(:)
                                ! For each variable:
                                !   T - create the variable as using time
                                !       dimension
                                !   F - create the variable as not using time
                                !       dimension
  INTEGER, INTENT(OUT) :: var_nlevs(:)
                               ! The number of vertical levels for each
                               ! variable
  INTEGER, INTENT(OUT) :: var_file_ids(:)
                               ! The id of each variable in the file


! Work variables
  INTEGER :: nvars  ! The number of variables to define

  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_levs_dims(SIZE(var_ids))
                               ! The name of the levels dimension for each
                               ! variable

  TYPE(dict) :: file_dims  ! Dictionary containing the dimensions to define
                           ! in the file
                           ! Maps dim_name => dim_size

  TYPE(dict) :: file_dim_ids  ! Dictionary containing the dimension ids in file
                              ! Maps dim_name => dim_id

  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name  ! Temporary variables used to
  INTEGER :: dim_size                          ! populate dimension dictionaries
  INTEGER :: dim_id
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name_ignore  ! Used as dummy argument

  INTEGER :: levels_dim_id  ! The id in the file of the vertical levels
                            ! dimension for the variable being processed
                            ! - used when defining variables

  INTEGER :: dummy  ! This is used to store the result from defining the time
                    ! dimension - we don't need to keep it

  TYPE(dict) :: int_attrs, real_attrs, char_attrs  ! Dictionaries containing
                                                   ! attribute values

  CHARACTER(len=DICT_KEY_LEN) :: key  ! Used when iterating over attribute
  INTEGER :: int_val                  ! dictionaries
  REAL :: real_val
  CHARACTER(len=DICT_CHAR_VAL_LEN) :: char_val

  INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  file_dims    = dict_create(MAX_DIM_FILE, INT(1))
  file_dim_ids = dict_create(MAX_DIM_FILE, INT(1))

!-----------------------------------------------------------------------------
! Check arguments
!-----------------------------------------------------------------------------
! Check that we have a value for var_names and use_time for every variable
  nvars = SIZE(var_ids)
  IF ( nvars /= SIZE(var_names) )                                             &
    CALL log_fatal("internal_define_vars",                                    &
                   "var_ids and var_names must have the same number " //      &
                   "of elements")
  IF ( nvars /= SIZE(use_time) )                                              &
    CALL log_fatal("internal_define_vars",                                    &
                   "var_ids and use_time must have the same number " //       &
                   "of elements")
! Check we have enough space to return all the variable ids and levels
  IF ( nvars /= SIZE(var_ids) )                                               &
    CALL log_fatal("internal_define_vars",                                    &
                   "var_ids must have the same number of elements as " //     &
                   "var_ids")
  IF ( nvars /= SIZE(var_nlevs) )                                             &
    CALL log_fatal("internal_define_vars",                                    &
                   "var_nlevs must have the same number of elements as " //   &
                   "var_ids")


!-----------------------------------------------------------------------------
! Define the grid
!-----------------------------------------------------------------------------
  CALL file_ts_def_grid(                                                      &
! Always use a 2d grid, even when ny=1 (for now)
    file, .FALSE., GRID_DIM_NAME, nx, X_DIM_NAME, nx, Y_DIM_NAME, ny          &
  )

!-----------------------------------------------------------------------------
! Define the required levels dimensions
!-----------------------------------------------------------------------------
! Make a pass over the requested variables to work out what dimensions we
! need to have in the file
  DO i = 1,nvars
! We ignore the dimension name for input files and use the name for output files
    CALL get_var_levs_dim(var_ids(i), dim_name_ignore, dim_name, dim_size)

! Store the levels dim name and size for this variable for use later
    var_levs_dims(i) = dim_name
    var_nlevs(i) = dim_size

! For a variable with no levels dimension, the returned dim_name will be the
! empty string - we don't want to add this to the dictionary
    IF ( LEN_TRIM(dim_name) < 1 ) CYCLE

! Otherwise store it in the dictionary
! Duplicate keys will overwritten, so we don't have to worry about duplication
    CALL dict_set(file_dims, dim_name, dim_size)
  END DO

! Define the required dimensions
  DO i = 1,file_dims%length
! The dimension names are the dictionary keys, and the sizes are the values
    dim_name = file_dims%keys(i)
    CALL dict_get(file_dims, dim_name, dim_size)

! Define the dimension
    dim_id = file_ts_def_dim(file, dim_name, dim_size)

! Store the id
    CALL dict_set(file_dim_ids, dim_name, dim_id)
  END DO

!-----------------------------------------------------------------------------
! Define the time dimension
!-----------------------------------------------------------------------------
  dummy = file_ts_def_time_dim(file, TIME_DIM_NAME)

!-----------------------------------------------------------------------------
! Define the variables and save their ids for later
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
! Just use the variable's string identifier as it's name in the file
    var_file_ids(i) = file_ts_def_var(                                        &
      file, var_names(i), levels_dim_id, use_time(i)                          &
    )

! Define attributes
    CALL get_var_attrs(var_ids(i), int_attrs, real_attrs, char_attrs)

! First the integer valued attributes
    DO j = 1,int_attrs%length
      key = int_attrs%keys(j)
      CALL dict_get(int_attrs, key, int_val)

      CALL file_ts_def_attr(file, var_file_ids(i), key, int_val)
    END DO

! Next, real valued attributes
    DO j = 1,real_attrs%length
      key = real_attrs%keys(j)
      CALL dict_get(real_attrs, key, real_val)

      CALL file_ts_def_attr(file, var_file_ids(i), key, real_val)
    END DO

! Lastly, character valued attributes
    DO j = 1,char_attrs%length
      key = char_attrs%keys(j)
      CALL dict_get(char_attrs, key, char_val)

      CALL file_ts_def_attr(file, var_file_ids(i), key, char_val)
    END DO
  END DO

  RETURN

END SUBROUTINE internal_define_vars
#endif
