#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_def_attr_real(file, var_id, name, value)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN

  USE driver_ascii_mod, ONLY : file_ascii_def_attr_real
  USE driver_ncdf_mod, ONLY : file_ncdf_def_attr_real

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
  TYPE(file_handle), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  REAL, INTENT(IN) :: value         ! The value of the attribute


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_def_attr_real(file%ascii, var_id, name, value)

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_def_attr_real(file%ncdf, var_id, name, value)

    CASE DEFAULT
      CALL log_fatal("file_def_attr_real",                                    &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_def_attr_real


SUBROUTINE file_def_attr_int(file, var_id, name, value)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN

  USE driver_ascii_mod, ONLY : file_ascii_def_attr_int
  USE driver_ncdf_mod, ONLY : file_ncdf_def_attr_int

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
  TYPE(file_handle), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  INTEGER, INTENT(IN) :: value      ! The value of the attribute


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_def_attr_int(file%ascii, var_id, name, value)

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_def_attr_int(file%ncdf, var_id, name, value)

    CASE DEFAULT
      CALL log_fatal("file_def_attr_int",                                     &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_def_attr_int


SUBROUTINE file_def_attr_char(file, var_id, name, value)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN

  USE driver_ascii_mod, ONLY : file_ascii_def_attr_char
  USE driver_ncdf_mod, ONLY : file_ncdf_def_attr_char

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
  TYPE(file_handle), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  CHARACTER(len=*), INTENT(IN) :: value
                                    ! The value of the attribute


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_def_attr_char(file%ascii, var_id, name, value)

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_def_attr_char(file%ncdf, var_id, name, value)

    CASE DEFAULT
      CALL log_fatal("file_def_attr_char",                                    &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_def_attr_char
#endif
