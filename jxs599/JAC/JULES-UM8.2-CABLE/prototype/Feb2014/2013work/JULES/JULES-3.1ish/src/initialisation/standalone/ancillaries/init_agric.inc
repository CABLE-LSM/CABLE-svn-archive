#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/ancillaries/init_agric.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_agric()

  USE io_constants, ONLY : MDI, MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN,          &
                           NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE input_mod, ONLY : nx, ny

  USE trifctl, ONLY : frac_agr_grid => frac_agr

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the agricultural fraction
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_agric namelist
! For data at a single point, a single fraction is specified in the namelist
! In all other cases, agricultural fraction is read from a file
!-----------------------------------------------------------------------------
  LOGICAL :: zero_agric           ! T - set agr. frac. at all points to 0.0
                                  ! F - read agr. frac. from input
  REAL :: frac_agr                ! Single point fraction

  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                                  ! The file to read fraction from
  CHARACTER(len=MAX_SDF_NAME_LEN) :: agric_name
                                  ! The name of the variable in the file
  NAMELIST /jules_agric/ zero_agric, frac_agr, file, agric_name


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  zero_agric = .TRUE.  ! Default is to set agricultural fraction to 0 everywhere
  frac_agr   = MDI     ! Initialised to missing data so we can tell if it is set
                       ! using the namelist
  file       = ''
  agric_name = ''

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_agric", "Reading JULES_AGRIC namelist...")

! First, we read the namelist
  READ(NAMELIST_UNIT, nml=jules_agric, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_agric",                                              &
                   "Error reading namelist JULES_AGRIC " //                   &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set values derived from namelist and verify for consistency
!-----------------------------------------------------------------------------
! If zero fraction is selected, then that is all we have to do
  IF ( zero_agric ) THEN
    CALL log_info("init_agric", "Zero agricultural fraction indicated")
    frac_agr_grid(:) = 0.0
    RETURN
  END IF

  IF ( nx * ny <= 1 ) THEN
!-----------------------------------------------------------------------------
! If we are running a single point, read fraction from the namelist
!-----------------------------------------------------------------------------
! Check if the namelist lat and lon have been set
    IF ( ABS(frac_agr - MDI) < EPSILON(frac_agr) )                            &
      CALL log_fatal("init_agric",                                            &
                     "For data at a single point, agricultural fraction " //  &
                     "is read from the namelist JULES_AGRIC")

! Now we know they have been set, copy their values into the model arrays
    CALL log_info("init_agric",                                               &
                  "Data is at a single point - reading agricultural " //      &
                  "fraction from namelist JULES_AGRIC")
    frac_agr_grid(1) = frac_agr
  ELSE
!-----------------------------------------------------------------------------
! If we have a grid, set fraction from the specified file
!-----------------------------------------------------------------------------
    CALL log_info("init_agric",                                               &
                  "Data is on a grid - reading agricultural fraction " //     &
                  "from file " // TRIM(file))
    CALL fill_variables_from_file(file, (/ 'frac_agr' /), (/ agric_name /))
  END IF

!-----------------------------------------------------------------------------
! Check that the values seem sensible
!-----------------------------------------------------------------------------
  IF ( ANY(frac_agr_grid < 0.0) .OR. ANY(frac_agr_grid > 1.0) )               &
    CALL log_fatal("init_agric",                                              &
                   "Agricultural fraction should be in range 0.0 to 1.0, " // &
                   "given range is " //                                       &
                   TRIM(to_string(MINVAL(frac_agr_grid))) // " to " //        &
                   TRIM(to_string(MAXVAL(frac_agr_grid))))

  RETURN

END SUBROUTINE init_agric
#endif
