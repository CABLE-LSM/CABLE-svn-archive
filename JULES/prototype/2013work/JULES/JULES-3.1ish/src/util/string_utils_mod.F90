#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/util/string_utils_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

MODULE string_utils_mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: MAX_STR_LEN = 50

  INTERFACE to_string
    MODULE PROCEDURE int_to_string, logical_to_string, real_to_string,        &
                     complex_to_string
  END INTERFACE


CONTAINS


  FUNCTION str_replace(string, search, replace) RESULT(str_repl)

    USE io_constants, ONLY : MAX_FILE_NAME_LEN

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Replaces all occurences of search in string with replace
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: string  ! The string to replace values in
    CHARACTER(len=*), INTENT(IN) :: search  ! The value to search for
    CHARACTER(len=*), INTENT(IN) :: replace ! The value to replace with

! Return type
    CHARACTER(len=MAX_FILE_NAME_LEN) :: str_repl  ! The string with occurences
                                                  ! of search replaced with
                                                  ! replace

! Work variables
    INTEGER :: search_len  ! The length of the search string
    CHARACTER(len=LEN(string)) :: string_local  ! Local copy of string that
                                                ! we can manipulate
    INTEGER :: idx  ! The index of the next occurence of search in string_local

!-----------------------------------------------------------------------------

    search_len = LEN_TRIM(search)

    str_repl = ""
    string_local = string

    DO
      idx = INDEX(string_local, TRIM(search))

      IF ( idx <= 0 ) THEN
! If index is 0, there are no more replacements to do, so just add on
! the rest of the string and exit
        str_repl = TRIM(str_repl) // TRIM(string_local)
        EXIT
      END IF

! If we get to here, then we have a replacement to do
      str_repl = TRIM(str_repl) // string_local(1:idx-1) // TRIM(replace)
      string_local = string_local(idx+search_len:)
    END DO

    RETURN

  END FUNCTION str_replace


  FUNCTION str_starts_with(string, sub_string) RESULT(starts_with)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns true if string starts with sub_string, false otherwise
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: string  ! The string to search in
    CHARACTER(len=*), INTENT(IN) :: sub_string  ! The string to match at the
                                                ! start

! Return type
    LOGICAL :: starts_with  ! T - the string starts with the given substring
                            ! F - otherwise


!-----------------------------------------------------------------------------

! Assume a negative result until we find otherwise
    starts_with = .FALSE.

! If the substring is longer than the string, it can't possibly match
    IF ( LEN_TRIM(sub_string) > LEN_TRIM(string) ) RETURN

    starts_with = ( string(1:LEN_TRIM(sub_string)) == sub_string )

    RETURN

  END FUNCTION str_starts_with


  SUBROUTINE str_split(string, separator, nparts, parts)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Splits the given string into parts around instances of separator. If
!   optional argument nparts is given, then it will contain the number of
!   parts the string is split into. If optional argument parts is given, it
!   will contain the parts that the string is split into - the caller must
!   ensure that there is enough space to do this.
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: string  ! The string to split
    CHARACTER(len=*), INTENT(IN) :: separator  !  The separator to split about
    INTEGER, INTENT(OUT), OPTIONAL :: nparts  ! The number of parts
    CHARACTER(len=*), INTENT(OUT), OPTIONAL :: parts(:)
                                            ! The constituent parts

! Work variables
    CHARACTER(len=LEN(string)) :: string_local  ! Local copy of string that
                                                ! we can manipulate
    INTEGER :: nparts_local  ! Local version of nparts that is always present
    INTEGER :: idx  ! The index of the next occurence of search in string_local


!-----------------------------------------------------------------------------

    nparts_local = 1  ! There is at least one part - the string itself if
                      ! separator is not present
    string_local = string

    DO
! Find the next occurence of separator in the string
      idx = INDEX(string_local, TRIM(separator))

      IF ( idx <= 0 ) THEN
! If index is 0, there are no more occurences , so just put the rest of
! the string in the final place in the parts array and be done with it
        IF ( PRESENT(parts) ) parts(nparts_local) = string_local
        EXIT
      END IF

! If we get to here, then we have a separator to deal with
! The string up to the separator goes into the array of parts
      IF ( PRESENT(parts) ) parts(nparts_local) = string_local(1:idx-1)
! Everything after the separator continues on for further analysis
      string_local = string_local(idx+LEN_TRIM(separator):)

      nparts_local = nparts_local + 1
    END DO

! Set the real nparts if present
    IF ( PRESENT(nparts) ) nparts = nparts_local

    RETURN

  END SUBROUTINE str_split


!------------------------------------------------------------------------------
! Procedures for formatting of variables of intrinsic data types
!------------------------------------------------------------------------------
  FUNCTION int_to_string(variable) RESULT(formatted_string)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: variable

    CHARACTER(len=MAX_STR_LEN) :: formatted_string

    WRITE(formatted_string, *) variable

    formatted_string = ADJUSTL(formatted_string)

  END FUNCTION int_to_string


  FUNCTION logical_to_string(variable) RESULT(formatted_string)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: variable

    CHARACTER(len=MAX_STR_LEN) :: formatted_string

    WRITE(formatted_string, *) variable

    formatted_string = ADJUSTL(formatted_string)

  END FUNCTION logical_to_string


  FUNCTION real_to_string(variable) RESULT(formatted_string)

    IMPLICIT NONE

    REAL, INTENT(IN) :: variable

    CHARACTER(len=MAX_STR_LEN) :: formatted_string

    WRITE(formatted_string, *) variable

    formatted_string = ADJUSTL(formatted_string)

  END FUNCTION real_to_string


  FUNCTION complex_to_string(variable) RESULT(formatted_string)

    IMPLICIT NONE

    COMPLEX, INTENT(IN) :: variable

    CHARACTER(len=MAX_STR_LEN) :: formatted_string

    CHARACTER(len=MAX_STR_LEN) :: real_part, imag_part

    WRITE(real_part, *) REAL(variable)
    WRITE(imag_part, *) AIMAG(variable)

    formatted_string = TRIM(ADJUSTL(real_part)) // ' + ' //                   &
                       TRIM(ADJUSTL(imag_part)) // 'i'

    formatted_string = ADJUSTL(formatted_string)

  END FUNCTION complex_to_string

END MODULE string_utils_mod
#endif
