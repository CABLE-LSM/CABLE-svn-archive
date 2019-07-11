#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/grid/init_land_frac.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_land_frac()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE coastal, ONLY : flandg

  USE input_mod, ONLY : nx, ny, fill_variables_from_file

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the land fraction on the full grid and checks for consistency
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
! Definition of the jules_land namelist - this combines variables
! from input_mod with some local variables
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: land_frac_name
  NAMELIST /jules_land_frac/ file, land_frac_name


!-----------------------------------------------------------------------------


  CALL log_info("init_land_frac", "Reading JULES_LAND_FRAC namelist...")

! First, we read the namelist
  READ(NAMELIST_UNIT, nml=jules_land_frac, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_land_frac",                                          &
                   "Error reading namelist JULES_LAND_FRAC " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Verify for consistency
!-----------------------------------------------------------------------------
  CALL log_info("init_land_frac",                                             &
                "Getting land fraction for the full input grid...")

! For now, we allocate flandg on the full input grid
  ALLOCATE(flandg(nx,ny), STAT=error)
  IF ( error > 0 )                                                            &
    CALL log_fatal("init_land_frac", "Error allocating flandg")

  IF ( nx * ny <= 1 ) THEN
! For single point data, just assume the single point is 100% land
    CALL log_info("init_land_frac",                                           &
                  "Data is at a single point - setting land fraction to 1.0")
    flandg(1,1) = 1.0
  ELSE
! For any other grid, read land fraction from file
    CALL log_info("init_land_frac",                                           &
                  "Data is on a grid - reading land fraction from file " //   &
                  TRIM(file))
    CALL fill_variables_from_file(file, (/ 'land_fraction' /),                &
                                        (/  land_frac_name   /))
  END IF

! For now, make sure that land fraction is either 1.0 or 0.0 - until JULES
! can deal with coastal tiling
  WHERE ( flandg > EPSILON(1.0) )
    flandg = 1.0
  ELSEWHERE
    flandg = 0.0
  END WHERE

  RETURN

END SUBROUTINE init_land_frac
#endif
