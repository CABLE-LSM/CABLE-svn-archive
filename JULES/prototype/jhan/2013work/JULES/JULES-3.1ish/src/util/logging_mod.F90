#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/util/logging_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE logging_mod

  USE io_constants, ONLY : UNIT_STDOUT

  IMPLICIT NONE

! Log levels - these can be combined using bitwise operators to indicate
! any combination of log levels
  INTEGER, PARAMETER ::                                                       &
    LOG_LEVEL_INFO  = 1,                                                      &
    LOG_LEVEL_DEBUG = 2,                                                      &
    LOG_LEVEL_WARN  = 4,                                                      &
    LOG_LEVEL_ERROR = 8,                                                      &
    LOG_LEVEL_FATAL = 16


! Determines what log levels are printed to log_unit - this is a bitwise
! combination of values from above.
! The default is to print everything
  INTEGER :: log_print_level = 31

! Determines what log levels cause a program to stop - this is a bitwise
! combination of values from above.
! Fatal errors will always cause the program to stop, by definition.
! The default is to stop only for fatal errors
! Setting this to 15 (i.e. stop for everything, even info) or 14 (i.e. stop
! for everything except info) are useful options for debugging
  INTEGER :: log_stop_level = 0

! The default is to print to stdout
  INTEGER :: log_unit = UNIT_STDOUT

CONTAINS

  SUBROUTINE log_init(unit, print_level, stop_level)

    IMPLICIT NONE

! Argument types
    INTEGER, OPTIONAL, INTENT(IN) :: unit
    INTEGER, OPTIONAL, INTENT(IN) :: print_level
    INTEGER, OPTIONAL, INTENT(IN) :: stop_level


! If unit is supplied, overwrite the default unit
    IF ( PRESENT(unit) ) log_unit = unit

! If print level has been supplied as an argument, overwrite the default
    IF ( PRESENT(print_level) ) log_print_level = print_level

! If stop level has been supplied as an argument, overwrite the default
    IF ( PRESENT(stop_level) ) log_stop_level = stop_level

  END SUBROUTINE log_init


  SUBROUTINE log_shutdown()

    IMPLICIT NONE

! Don't bother catching any errors, as there is nothing sensible we can do
! with them now (like log them...)!!
    CLOSE(log_unit)

  END SUBROUTINE log_shutdown


  SUBROUTINE write_to_log(log_level, message)

    IMPLICIT NONE

! Argument types
    INTEGER, INTENT(IN) :: log_level
      ! The level at which to log the given message
    CHARACTER(len=*), INTENT(IN) :: message
      ! The message to log


!-----------------------------------------------------------------------------

! Use a bitwise and to check if we want to print log messages for the
! given log level
    IF ( IAND(log_print_level, log_level) > 0 ) THEN
      WRITE(log_unit, "(A)") message
    END IF

! Check if we need to stop the program
    IF ( IAND(log_stop_level, log_level) > 0 .OR.                             &
         log_level == LOG_LEVEL_FATAL ) THEN
      CALL log_shutdown()
      STOP
    END IF

  END SUBROUTINE write_to_log


  SUBROUTINE log_info(proc_name, message)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(                                                        &
      LOG_LEVEL_INFO, '[INFO] ' // proc_name // ': ' // message               &
    )

  END SUBROUTINE log_info


  SUBROUTINE log_debug(proc_name, message)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(                                                        &
      LOG_LEVEL_DEBUG, '[DEBUG] ' // proc_name // ': ' // message             &
    )

  END SUBROUTINE log_debug


  SUBROUTINE log_warn(proc_name, message)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(                                                        &
      LOG_LEVEL_WARN, '[WARNING] ' // proc_name // ': ' // message            &
    )

  END SUBROUTINE log_warn


  SUBROUTINE log_error(proc_name, message)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(                                                        &
      LOG_LEVEL_ERROR, '[ERROR] ' // proc_name // ': ' // message             &
    )

  END SUBROUTINE log_error


  SUBROUTINE log_fatal(proc_name, message)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(                                                        &
      LOG_LEVEL_FATAL, '[FATAL ERROR] ' // proc_name // ': ' // message       &
    )

  END SUBROUTINE log_fatal

END MODULE logging_mod
#endif
