#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/grid/init_surf_hgt.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_surf_hgt()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE input_mod, ONLY : nx, ny, fill_variables_from_file

  USE switches, ONLY : l_aggregate

  USE c_elevate, ONLY : surf_hgt, surf_hgt_io

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the surface height
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
! Definition of the jules_surf_hgt namelist
!-----------------------------------------------------------------------------
  LOGICAL :: zero_height  ! T - set height at all points to 0.0
                          ! F - read height from file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file  ! If input grid has more than one
                                            ! point, read heights from this
                                            ! file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: surf_hgt_name  ! Use the variable with
                                                    ! this name
  NAMELIST /jules_surf_hgt/ zero_height, surf_hgt_io, file, surf_hgt_name


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  zero_height = .TRUE.  ! Default is to use zero height everywhere

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_surf_hgt", "Reading JULES_SURF_HGT namelist...")

! First, we read the namelist
  READ(NAMELIST_UNIT, nml=jules_surf_hgt, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_surf_hgt",                                           &
                   "Error reading namelist JULES_SURF_HGT " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set values derived from namelist and verify for consistency
!-----------------------------------------------------------------------------
! Zero height must be used with the aggregate surface scheme
  IF ( l_aggregate ) THEN
    CALL log_warn("init_surf_hgt",                                            &
                  "Aggregate surface selected, so using zero height")
    zero_height = .TRUE.
  END IF

! If zero height is selected, then that is all we have to do
  IF ( zero_height ) THEN
    CALL log_info("init_surf_hgt",                                            &
                  "Zero height selected - setting all heights to 0.0")
    surf_hgt(:,:) = 0.0
    RETURN
  END IF


  IF ( nx * ny <= 1 ) THEN
!-----------------------------------------------------------------------------
! If we are reading data at a single point, read height from the namelist
!-----------------------------------------------------------------------------
! Copy their values into the model arrays
    CALL log_info("init_surf_hgt",                                            &
                  "Data is at a single point - reading surf_hgt from " //     &
                  "surf_hgt_io in namelist JULES_SURF_HGT")
    surf_hgt(1,:) = surf_hgt_io(:)
  ELSE
!-----------------------------------------------------------------------------
! If we have a grid, set heights from the specified file
!-----------------------------------------------------------------------------
    CALL log_info("init_surf_hgt",                                            &
                  "Data is on a grid - reading surf_hgt from file " //        &
                  TRIM(file))
    CALL fill_variables_from_file(file, (/ 'surf_hgt' /), (/ surf_hgt_name /))
  END IF

  RETURN

END SUBROUTINE init_surf_hgt
#endif
