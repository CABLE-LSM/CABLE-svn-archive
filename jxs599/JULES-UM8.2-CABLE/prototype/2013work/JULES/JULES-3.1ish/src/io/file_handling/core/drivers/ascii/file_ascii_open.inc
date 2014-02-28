#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/drivers/ascii/file_ascii_open.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_ascii_open(name, mode) RESULT(file)

  USE io_constants, ONLY : MODE_READ, MODE_WRITE, UNIT_STDIN, UNIT_STDOUT

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a file and returns a file_ascii object representing it
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*), INTENT(IN) :: name
                              ! The name of the file
  INTEGER, INTENT(IN) :: mode ! The mode to open the file
                              ! One of MODE_READ or MODE_WRITE

! Return type
  TYPE(file_ascii) :: file


! Work variables
  INTEGER :: unit ! The unit number for the opened file
  INTEGER :: error ! Error code for any errors that occur

  LOGICAL :: unit_in_use ! Used in call to INQUIRE to detect if a unit is in use
  INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Get a spare unit to open the file
!-----------------------------------------------------------------------------
  unit = 0
! Search the allowed unit numbers until we find one with nothing connected
  DO i = 1,MAX_UNIT_NUMBER
! Avoid units for standard i/o
    IF ( i == UNIT_STDIN .OR. i == UNIT_STDOUT ) CYCLE

! Find out if anything is connected to this unit.
    INQUIRE(unit=i, opened=unit_in_use )
    IF ( .NOT. unit_in_use ) THEN
! This unit is free, so exit the loop
      unit = i
      EXIT
    END IF
  END DO

  IF ( unit < 1 ) THEN
    CALL log_fatal("file_ascii_open",                                         &
                   "All allowed units are in use - try increasing " //        &
                   "MAX_UNIT_NUMBER")
  END IF


!-----------------------------------------------------------------------------
! Open the file in the requested mode
!-----------------------------------------------------------------------------
  SELECT CASE ( mode )
    CASE ( MODE_READ )
      CALL log_info("file_ascii_open",                                        &
                    "Opening file " // TRIM(name) // " for reading")
! Open file for reading only - file must already exist
      OPEN(unit, FILE=name, STATUS='old', POSITION='rewind', ACTION='read',   &
                 IOSTAT=error)

    CASE ( MODE_WRITE )
      CALL log_info("file_ascii_open",                                        &
                    "Opening file " // TRIM(name) // " for writing")
! Create an empty file for writing only - if a file with the
! same name exists, overwrite it
      OPEN(unit, FILE=name, STATUS='replace', POSITION='rewind',              &
                 ACTION='write', IOSTAT=error)

    CASE DEFAULT
! Read and write are the only supported modes
      CALL log_fatal("file_ascii_open",                                       &
                     "Unsupported mode - " // TRIM(to_string(mode)))

  END SELECT

!-----------------------------------------------------------------------------
! Report any errors
!-----------------------------------------------------------------------------
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_open",                                         &
                   "Error opening file " // TRIM(name) //                     &
                   " (IOSTAT=" // TRIM(to_string(error)) // ")")


! Initialise the file_ascii object
  file%name = name
  file%mode = mode
  file%unit = unit

  RETURN

END FUNCTION file_ascii_open
#endif
