#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/util/templating_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

MODULE templating_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
! Strings used to indicate templated parts of strings
! These are replaced by functions in this module to create other strings
  CHARACTER(len=3), PARAMETER ::                                              &
    TPL_YR_2DIGIT  = '%y2',                                                   &
               !  2-digit year
    TPL_YR_4DIGIT  = '%y4',                                                   &
               !  4-digit year
    TPL_MON        = '%m1',                                                   &
               !  1- or 2-digit month
    TPL_MON_2DIGIT = '%m2',                                                   &
               !  2-digit month
    TPL_MON_ABBR   = '%mc',                                                   &
               !  3-character month abbreviation
    TPL_VAR_NAME   = '%vv'
               !  Name of a variable

  CHARACTER(len=3) :: MONTH_ABBREVIATIONS(12)

  DATA MONTH_ABBREVIATIONS / 'jan', 'feb', 'mar', 'apr', 'may', 'jun',        &
                             'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /


CONTAINS


!*****************************************************************************
! Variable name templating
!*****************************************************************************
  LOGICAL FUNCTION tpl_has_var_name(tpl)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Checks if a template includes variable name templating
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: tpl  ! The string to check for variable
                                         ! name templating

!-----------------------------------------------------------------------------

! We just check for the variable name templating string in the given template
    tpl_has_var_name = ( INDEX(tpl, TPL_VAR_NAME) > 0 )

    RETURN

  END FUNCTION tpl_has_var_name


  FUNCTION tpl_substitute_var(tpl, var_name) RESULT(file_name)

    USE io_constants, ONLY : MAX_FILE_NAME_LEN

    USE string_utils_mod, ONLY : str_replace

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Replaces all occurences of the variable name template string with the
!   given variable name
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: tpl  ! The string to replace variable
                                         ! name templates in
    CHARACTER(len=*), INTENT(IN) :: var_name  ! The variable name to replace
                                              ! them with

! Return type
    CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name  ! The resulting file name


!-----------------------------------------------------------------------------

    file_name = str_replace(tpl, TPL_VAR_NAME, var_name)

    RETURN

  END FUNCTION tpl_substitute_var


!*****************************************************************************
! Time templating
!*****************************************************************************
  FUNCTION tpl_detect_period(tpl) RESULT(tpl_period)

    USE datetime_mod, ONLY : PERIOD_YEAR, PERIOD_MONTH

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Determines from the template the period of time over which files will
!   apply
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: tpl  ! The template to determine period
                                         ! for

! Return type
    INTEGER :: tpl_period  ! The detected templating period

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Currently, we only have code for handling templating periods of a month
! or a year, so just look for those
!-----------------------------------------------------------------------------

! Return a templating period of 0 if there are no time templating strings
    tpl_period = 0

    IF ( INDEX(tpl, TPL_YR_2DIGIT) > 0 .OR. INDEX(tpl, TPL_YR_4DIGIT) > 0 ) THEN
! If the template contains any year strings, then the templating period is
! at least yearly
      tpl_period = PERIOD_YEAR

! We can't have monthly templating without first having yearly templating...
      IF ( INDEX(tpl, TPL_MON) > 0        .OR.                                &
           INDEX(tpl, TPL_MON_2DIGIT) > 0 .OR.                                &
           INDEX(tpl, TPL_MON_ABBR) > 0 )                                     &
! If the template contains any month strings, then the templating period is
! monthly
        tpl_period = PERIOD_MONTH

    END IF  ! year

    RETURN

  END FUNCTION tpl_detect_period


  FUNCTION tpl_substitute_datetime(tpl, dt) RESULT(file_name)

    USE io_constants, ONLY : MAX_FILE_NAME_LEN

    USE string_utils_mod, ONLY : str_replace

    USE datetime_mod, ONLY : datetime

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Replaces all occurences of the time templating strings with values from
!   the given datetime
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    CHARACTER(len=*), INTENT(IN) :: tpl
                              ! The template to replace templating strings in
    TYPE(datetime), INTENT(IN) :: dt
                              ! The datetime to get values from

! Return type
    CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name
                              ! The resulting file name


! Work variables
    CHARACTER(len=4) :: formatted_str  ! The formatted string to pass to
                                       ! the replace function

!-----------------------------------------------------------------------------

! Just replace each time templating string appropriately

! 4 digit years
    WRITE(formatted_str, "(I4.4)") dt%year
    file_name = str_replace(tpl, TPL_YR_4DIGIT, formatted_str)

! 2 digit years
    WRITE(formatted_str, "(I2.2)") MOD(dt%year, 100)
    file_name = str_replace(file_name, TPL_YR_2DIGIT, TRIM(formatted_str))

! 1 or 2 digit months
    IF ( dt%month < 10 ) THEN
      WRITE(formatted_str, "(I1.1)") dt%month
    ELSE
      WRITE(formatted_str, "(I2.2)") dt%month
    ENDIF
    file_name = str_replace(file_name, TPL_MON, TRIM(formatted_str))

! 2 digit month
    WRITE(formatted_str, '(I2.2)') dt%month
    file_name = str_replace(file_name, TPL_MON_2DIGIT, TRIM(formatted_str))

! 3 character string month
    formatted_str = MONTH_ABBREVIATIONS(dt%month)
    file_name = str_replace(file_name, TPL_MON_ABBR, TRIM(formatted_str))

    RETURN

  END FUNCTION tpl_substitute_datetime

END MODULE templating_mod
#endif
