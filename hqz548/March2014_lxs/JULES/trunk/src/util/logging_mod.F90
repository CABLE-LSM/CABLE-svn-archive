#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE logging_mod

  USE io_constants, ONLY : UNIT_STDOUT

  IMPLICIT NONE

  INCLUDE "mpif.h"

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
! The default (31) is to print everything
  INTEGER :: log_print_level = 31

! Determines what log levels cause a program to stop - this is a bitwise
! combination of values from above.
! Fatal errors will always cause the program to stop, by definition.
! The default (0) is to stop only for fatal errors
! Setting this to 15 (i.e. stop for everything, even info) or 14 (i.e. stop
! for everything except info) are useful options for debugging
  INTEGER :: log_stop_level = 0

! This is the unit that log files will be opened on when there are multiple
! MPI tasks available
  INTEGER, PARAMETER :: LOG_FILE_UNIT = 90

! The maximum line length for log entries. Any message larger than this is
! truncated
  INTEGER, PARAMETER :: LOG_MAX_LINE_LEN = 1000

! This is the unit that log messages will be written to
! The default is STDOUT - this is only changed if there are multiple MPI
! tasks that might call the logging
  INTEGER :: log_unit = UNIT_STDOUT

! Prefix to attach to log messages indicating the task that generated the
! message
  CHARACTER(len=35) :: task_prefix = ""


! Visibilities
  PRIVATE
  PUBLIC :: log_init, log_shutdown, log_info, log_debug, log_warn, log_error, &
            log_fatal

CONTAINS

  SUBROUTINE log_init(nml_dir)

    USE io_constants, ONLY : MAX_FILE_NAME_LEN, NAMELIST_UNIT

    USE string_utils_mod, ONLY : to_string

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the logging environment
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
    CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                             ! namelist for logging configuration

    CHARACTER(len=MAX_FILE_NAME_LEN) :: log_dir  ! The directory that log files
                                                 ! will be created in

! The namelist also sets the stop and print levels
    NAMELIST /logging/ log_dir, log_print_level, log_stop_level

! Work variables
    CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name
    CHARACTER(len=20) :: task_id_str

    INTEGER :: ntasks, task_id, error, i

    LOGICAL :: dir_exists  ! Used when checking for existence of the logging
                           ! directory


!-----------------------------------------------------------------------------


    log_dir = ""  ! Default log directory is nothing

!-----------------------------------------------------------------------------
! Note that before log_init is called, the logging environment is:
!
!   * All log output goes to STDOUT
!
!   * All log messages are printed (log_print_level = 31)
!
!   * Only fatal messages cause the program to stop (log_stop_level = 1)
!
!   * There is no task prefix (i.e. no indicator of which task generated
!     output)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Start by setting up as if all output from all tasks will go to STDOUT
! (it will to start with!)
!-----------------------------------------------------------------------------
! Find out information about our operating environment
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, error)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, task_id, error)

! Get the task id as a string
! If we ever have a number of tasks that needs more than 20 characters to print,
! then that's awesome!!
    WRITE(task_id_str, "(I20)") task_id
    task_id_str = ADJUSTL(task_id_str)

! Set up the task prefix if there is more than one task
    IF ( ntasks > 1 )                                                         &
      task_prefix = "{MPI Task " // TRIM(task_id_str) // "} "

!-----------------------------------------------------------------------------
! Read the logging configuration from the namelist
!
! Note that any logging we do at this stage will go to the default unit of
! STDOUT
!-----------------------------------------------------------------------------
    OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'logging.nml'),         &
                        STATUS='old', POSITION='rewind', ACTION='read',       &
                        IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("log_init",                                              &
                     "Error opening namelist file logging.nml " //            &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

    READ(NAMELIST_UNIT, nml=logging, IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("log_init",                                              &
                     "Error reading namelist LOGGING " //                     &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

    CLOSE(NAMELIST_UNIT, IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("log_init",                                              &
                     "Error closing namelist file logging.nml " //            &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

! If a log directory has been supplied that doesn't exist, error out
    IF ( LEN_TRIM(log_dir) > 0 ) THEN
! The Intel compiler requires a different form of the statement for directories,
! which we swap in with an ifdef
#if defined(COMPILER_INTEL)
      INQUIRE(DIRECTORY=log_dir, EXIST=dir_exists)
#else
      INQUIRE(FILE=log_dir, EXIST=dir_exists)
#endif

      IF ( .NOT. dir_exists )                                                 &
        CALL log_fatal("log_init",                                            &
                       "Log directory does not exist - " // TRIM(log_dir))
    END IF

!-----------------------------------------------------------------------------
! Set up the logging environment for this task
!
!   * If the user has supplied a log directory, create a file for each task
!     in the log directory
!
!   * If the user has not supplied a log directory, all log output continues
!     to go to STDOUT
!-----------------------------------------------------------------------------
    IF ( LEN_TRIM(log_dir) > 0 ) THEN
!-----------------------------------------------------------------------------
! If a log directory has been specified, open a log file for each MPI task
! (even if there is only one task)
!
! We leave the task prefix for messages as blank, as the task that generated
! a message is obvious from the file it is contained in
!-----------------------------------------------------------------------------
! Open the log file for this task
      file_name = TRIM(log_dir) // '/' //                                     &
                  "task" // TRIM(ADJUSTL(task_id_str)) // ".stdout"
      OPEN(UNIT=LOG_FILE_UNIT, FILE=file_name, IOSTAT=error,                  &
           STATUS='REPLACE', ACTION='WRITE', RECL=LOG_MAX_LINE_LEN)

! Check if we managed to open the file
      IF ( error /= 0 )                                                       &
! If we call log_fatal here, it will write to STDOUT, since we have not changed
! the unit yet
        CALL log_fatal("log_init",                                            &
                       'Error opening log file ' // TRIM(file_name))

! Before we change the log_unit, print info to STDOUT about which log file
! each task is using
! We loop over the available tasks with a barrier before the next iteration starts
! This means each task will print the file it is using to STDOUT in order before
! execution continues
      DO i = 0,ntasks-1
! Print the file name if the loop is on our task
        IF ( task_id == i )                                                   &
          CALL log_info('log_init',                                           &
                        'Output from task ' // TRIM(ADJUSTL(task_id_str)) //  &
                        ' will be written to ' // TRIM(file_name))

        CALL MPI_BARRIER(MPI_COMM_WORLD, error)
      END DO

! Set the log unit
      log_unit = LOG_FILE_UNIT
! Since from now on all log output will go into a task specific file, the
! task prefix is no longer necessary
      task_prefix = ""
    END IF

! We are already set up to use STDOUT if a log directory is not given

    RETURN

  END SUBROUTINE log_init


  SUBROUTINE log_shutdown()

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Does any cleanup required by the logging code
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Don't bother catching any errors, as there is nothing sensible we can do
! with them now (like log them...)!!
    CLOSE(log_unit)

    RETURN

  END SUBROUTINE log_shutdown


  SUBROUTINE write_to_log(log_level, proc_name, message)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Logs the given message at the given level and performs any necessary
!   action
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
    INTEGER, INTENT(IN) :: log_level
      ! The level at which to log the given message
    CHARACTER(len=*), INTENT(IN) :: proc_name
      ! The name of the originating routine/function
    CHARACTER(len=*), INTENT(IN) :: message
      ! The message to log


! Work variables
    CHARACTER(len=LOG_MAX_LINE_LEN) :: full_message
    INTEGER :: error  ! Placeholder for MPI error code


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Construct the full log message
!-----------------------------------------------------------------------------
! Combine the name of the originating procedure with the message
    full_message = TRIM(proc_name) // ': ' // message

! Prepend a fragment to the message depending on the log level
    SELECT CASE ( log_level )
      CASE ( LOG_LEVEL_INFO )
        full_message = "[INFO] " // full_message

      CASE ( LOG_LEVEL_DEBUG )
        full_message = "[DEBUG] " // full_message

      CASE ( LOG_LEVEL_WARN )
        full_message = "[WARNING] " // full_message

      CASE ( LOG_LEVEL_ERROR )
        full_message = "[ERROR] " // full_message

      CASE ( LOG_LEVEL_FATAL )
        full_message = "[FATAL ERROR] " // full_message

      CASE DEFAULT
! This should never happen since the only access to write_to_log is through
! the log_* routines defined below
        CALL log_fatal("write_to_log", "Unknown log level")
    END SELECT

! Prepend the task prefix
    full_message = TRIM(task_prefix) // full_message

!-----------------------------------------------------------------------------
! Use a bitwise and to check if we want to print log messages for the
! given log level
!-----------------------------------------------------------------------------
    IF ( IAND(log_print_level, log_level) > 0 ) THEN
      WRITE(log_unit, "(A)") TRIM(full_message)
    END IF

!-----------------------------------------------------------------------------
! Check if we need to stop the program
!-----------------------------------------------------------------------------
    IF ( IAND(log_stop_level, log_level) > 0 .OR.                             &
         log_level == LOG_LEVEL_FATAL ) THEN
      CALL log_shutdown()
! Abort MPI with a non-zero exit code
      CALL MPI_ABORT(MPI_COMM_WORLD, 1, error)
    END IF

    RETURN

  END SUBROUTINE write_to_log


  SUBROUTINE log_info(proc_name, message)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level info
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(LOG_LEVEL_INFO, proc_name, message)

    RETURN

  END SUBROUTINE log_info


  SUBROUTINE log_debug(proc_name, message)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level debug
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(LOG_LEVEL_DEBUG, proc_name, message)

    RETURN

  END SUBROUTINE log_debug


  SUBROUTINE log_warn(proc_name, message)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level warn
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(LOG_LEVEL_WARN, proc_name, message)

    RETURN

  END SUBROUTINE log_warn


  SUBROUTINE log_error(proc_name, message)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level error
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(LOG_LEVEL_ERROR, proc_name, message)

    RETURN

  END SUBROUTINE log_error


  SUBROUTINE log_fatal(proc_name, message)

    IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Utility function to write a message to the log at level fatal
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

    CHARACTER(len=*), INTENT(IN) :: proc_name
    CHARACTER(len=*), INTENT(IN) :: message


    CALL write_to_log(LOG_LEVEL_FATAL, proc_name, message)

    RETURN

  END SUBROUTINE log_fatal

END MODULE logging_mod
#endif
