#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_def_attr_real(file, var_id, name, value)

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
  TYPE(file_ascii), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  REAL, INTENT(IN) :: value         ! The value of the attribute


!-----------------------------------------------------------------------------

! Since all attributes on ASCII files will eventually be char when written
! to file, we just delegate to def_attr_char with a converted value
  CALL file_ascii_def_attr_char(file, var_id, name, TRIM(to_string(value)))

  RETURN

END SUBROUTINE file_ascii_def_attr_real


SUBROUTINE file_ascii_def_attr_int(file, var_id, name, value)

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
  TYPE(file_ascii), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  INTEGER, INTENT(IN) :: value      ! The value of the attribute


!-----------------------------------------------------------------------------

! Since all attributes on ASCII files will eventually be char when written
! to file, we just delegate to def_attr_char with a converted value
  CALL file_ascii_def_attr_char(file, var_id, name, TRIM(to_string(value)))

  RETURN

END SUBROUTINE file_ascii_def_attr_int


SUBROUTINE file_ascii_def_attr_char(file, var_id, name, value)

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
  TYPE(file_ascii), INTENT(INOUT) :: file
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
    CALL log_fatal("file_ascii_def_attr_char",                                &
                   "Cannot define attribute - file is not in define mode")

! If adding another attribute will cause us to have too many attributes,
! error out
  IF ( file%nattrs >= MAX_ATTR_FILE )                                         &
    CALL log_fatal("file_ascii_def_attr_char",                                &
                   "Too many attributes in file - try increasing MAX_ATTR_FILE")

!-----------------------------------------------------------------------------
! Store information about the attribute so that it can be used later
!-----------------------------------------------------------------------------
  file%nattrs = file%nattrs + 1

  file%attr_var_ids(file%nattrs) = var_id
! We trim trailing spaces from name and leading spaces from value so that
! they will be the correct distance apart. Then left align the whole lot
  file%attr_values(file%nattrs)  = ADJUSTL(TRIM(name) // " = " // ADJUSTL(value))

  RETURN

END SUBROUTINE file_ascii_def_attr_char
#endif
