#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_def_attr_real(file, var_id, name, value)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, ATTR_GLOBAL

  USE file_mod, ONLY : file_def_attr

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
  TYPE(file_gridded), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  REAL, INTENT(IN) :: value         ! The value of the attribute

! Work variables
  INTEGER :: var_id_local


!-----------------------------------------------------------------------------

! Look up the variable id in the underlying file, unless a global attribute
! has been requested
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
  var_id_local = var_id
  IF ( var_id /= ATTR_GLOBAL ) var_id_local = file%vars(var_id)%id

! Now just defer to the underlying file
  CALL file_def_attr(file%fh, var_id_local, name, value)

  RETURN

END SUBROUTINE file_gridded_def_attr_real


SUBROUTINE file_gridded_def_attr_int(file, var_id, name, value)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, ATTR_GLOBAL

  USE file_mod, ONLY : file_def_attr

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
  TYPE(file_gridded), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  INTEGER, INTENT(IN) :: value      ! The value of the attribute

! Work variables
  INTEGER :: var_id_local


!-----------------------------------------------------------------------------

! Look up the variable id in the underlying file, unless a global attribute
! has been requested
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
  var_id_local = var_id
  IF ( var_id /= ATTR_GLOBAL ) var_id_local = file%vars(var_id)%id

! Now just defer to the underlying file
  CALL file_def_attr(file%fh, var_id_local, name, value)

  RETURN

END SUBROUTINE file_gridded_def_attr_int


SUBROUTINE file_gridded_def_attr_char(file, var_id, name, value)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, ATTR_GLOBAL

  USE file_mod, ONLY : file_def_attr

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
  TYPE(file_gridded), INTENT(INOUT) :: file
                                    ! The file to define the attribute in
  INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                    ! attribute on
  CHARACTER(len=*), INTENT(IN) :: name
                                    ! The name of the attribute
  CHARACTER(len=*), INTENT(IN) :: value
                                    ! The value of the attribute

! Work variables
  INTEGER :: var_id_local


!-----------------------------------------------------------------------------

! Look up the variable id in the underlying file, unless a global attribute
! has been requested
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
  var_id_local = var_id
  IF ( var_id /= ATTR_GLOBAL ) var_id_local = file%vars(var_id)%id

! Now just defer to the underlying file
  CALL file_def_attr(file%fh, var_id_local, name, value)

  RETURN

END SUBROUTINE file_gridded_def_attr_char
#endif
