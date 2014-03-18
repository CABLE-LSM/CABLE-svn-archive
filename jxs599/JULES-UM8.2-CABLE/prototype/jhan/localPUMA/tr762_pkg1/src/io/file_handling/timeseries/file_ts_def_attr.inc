#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_def_attr_real(file, var_id, name, value)

  USE io_constants, ONLY : ATTR_GLOBAL

  USE dictionary_mod, ONLY : dict_create, dict_set

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Define a real valued attribute on the given variable with the given name
!   and value
!   To define a global attribute, specify ATTR_GLOBAL as var_id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  REAL, INTENT(IN) :: value         ! The value of the attribute


!-----------------------------------------------------------------------------


! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ts_def_attr_real",                                   &
                   "Cannot define attribute - file is not in define mode")

!-----------------------------------------------------------------------------
! Work out what dictionary we want to populate
!-----------------------------------------------------------------------------
  IF ( var_id == ATTR_GLOBAL ) THEN
! If it is a global attribute, use the global dictionary
! Create the required dictionary if it has not been used yet
    IF ( file%attrs_real%length == 0 )                                        &
      file%attrs_real = dict_create(MAX_ATTR_FILE, 1.0)

    CALL dict_set(file%attrs_real, name, value)
  ELSE
! Otherwise, use the dictionary on the specified variable
! Create the required dictionary if it has not been used yet
    IF ( file%vars(var_id)%attrs_real%length == 0 )                           &
      file%vars(var_id)%attrs_real = dict_create(MAX_ATTR_FILE, 1.0)

    CALL dict_set(file%vars(var_id)%attrs_real, name, value)
  END IF

  RETURN

END SUBROUTINE file_ts_def_attr_real


SUBROUTINE file_ts_def_attr_int(file, var_id, name, value)

  USE io_constants, ONLY : ATTR_GLOBAL

  USE dictionary_mod, ONLY : dict_create, dict_set

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Define an integer valued attribute on the given variable with the given
!   name and value
!   To define a global attribute, specify ATTR_GLOBAL as var_id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  INTEGER, INTENT(IN) :: value      ! The value of the attribute


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ts_def_attr_int",                                    &
                   "Cannot define attribute - file is not in define mode")

!-----------------------------------------------------------------------------
! Work out what dictionary we want to populate
!-----------------------------------------------------------------------------
  IF ( var_id == ATTR_GLOBAL ) THEN
! If it is a global attribute, use the global dictionary
! Create the required dictionary if it has not been used yet
    IF ( file%attrs_int%length == 0 )                                         &
      file%attrs_int = dict_create(MAX_ATTR_FILE, INT(1))

    CALL dict_set(file%attrs_int, name, value)
  ELSE
! Otherwise, use the dictionary on the specified variable
! Create the required dictionary if it has not been used yet
    IF ( file%vars(var_id)%attrs_int%length == 0 )                            &
      file%vars(var_id)%attrs_int = dict_create(MAX_ATTR_FILE, INT(1))

    CALL dict_set(file%vars(var_id)%attrs_int, name, value)
  END IF

  RETURN

END SUBROUTINE file_ts_def_attr_int


SUBROUTINE file_ts_def_attr_char(file, var_id, name, value)

  USE io_constants, ONLY : ATTR_GLOBAL

  USE dictionary_mod, ONLY : dict_create, dict_set

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Define a character valued attribute on the given variable with the given
!   name and value
!   To define a global attribute, specify ATTR_GLOBAL as var_id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  CHARACTER(len=*), INTENT(IN) :: value
                                    ! The value of the attribute


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ts_def_attr_char",                                   &
                   "Cannot define attribute - file is not in define mode")

!-----------------------------------------------------------------------------
! Work out what dictionary we want to populate
!-----------------------------------------------------------------------------
  IF ( var_id == ATTR_GLOBAL ) THEN
! If it is a global attribute, use the global dictionary
! Create the required dictionary if it has not been used yet
    IF ( file%attrs_char%length == 0 )                                        &
      file%attrs_char = dict_create(MAX_ATTR_FILE, name)

    CALL dict_set(file%attrs_char, name, value)
  ELSE
! Otherwise, use the dictionary on the specified variable
! Create the required dictionary if it has not been used yet
    IF ( file%vars(var_id)%attrs_char%length == 0 )                           &
      file%vars(var_id)%attrs_char = dict_create(MAX_ATTR_FILE, name)

    CALL dict_set(file%vars(var_id)%attrs_char, name, value)
  END IF

  RETURN

END SUBROUTINE file_ts_def_attr_char
#endif
