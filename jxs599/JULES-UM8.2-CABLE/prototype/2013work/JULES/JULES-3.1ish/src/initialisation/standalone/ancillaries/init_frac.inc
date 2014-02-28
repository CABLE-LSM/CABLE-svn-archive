#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/ancillaries/init_frac.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_frac()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_SDF_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE switches, ONLY : l_veg_compete

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the tile fractions and checks for consistency
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
! Definition of the jules_frac namelist
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: frac_name
  NAMELIST /jules_frac/ file, frac_name


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  file       = ''
  frac_name = ''

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_frac", "Reading JULES_FRAC namelist...")

! First, we read the namelist
  READ(NAMELIST_UNIT, nml=jules_frac, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_frac",                                               &
                   "Error reading namelist JULES_FRAC " //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set frac using namelist values
!-----------------------------------------------------------------------------
! Frac is prognostic if competing veg is on, so must be read as an initial
! condition
  IF ( l_veg_compete ) THEN
    CALL log_info("init_frac",                                                &
                  "Competing vegetation is enabled - frac will be read " //   &
                  "as an initial condition")
    RETURN
  END IF

  CALL log_info("init_frac",                                                  &
                "Reading tile fractions from file " // TRIM(file))
  CALL fill_variables_from_file(file, (/ 'frac' /), (/ frac_name /))

  RETURN

END SUBROUTINE init_frac
#endif
