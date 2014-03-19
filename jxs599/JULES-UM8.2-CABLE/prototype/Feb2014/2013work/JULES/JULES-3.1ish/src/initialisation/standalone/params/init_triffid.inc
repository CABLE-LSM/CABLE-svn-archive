#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/params/init_triffid.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_triffid()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE nstypes, ONLY : npft

  USE switches, ONLY : l_phenol, l_triffid

  USE trif

  USE trif_io

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the TRIFFID parameters
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

! Nothing to do if neither triffid or phenology are selected
  IF ( .NOT. l_triffid .AND. .NOT. l_phenol ) RETURN


  CALL log_info("init_triffid", "Reading JULES_TRIFFID namelist...")

! Open the pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE='triffid_params.nml',                              &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Error opening namelist file triffid_params.nml " //       &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_triffid, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Error reading namelist JULES_TRIFFID " //                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Error closing namelist file triffid_params.nml " //       &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays
!-----------------------------------------------------------------------------
  crop(:)    = crop_io(1:npft)
  g_area(:)  = g_area_io(1:npft)
  g_grow(:)  = g_grow_io(1:npft)
  g_root(:)  = g_root_io(1:npft)
  g_wood(:)  = g_wood_io(1:npft)
  lai_max(:) = lai_max_io(1:npft)
  lai_min(:) = lai_min_io(1:npft)

  RETURN

END SUBROUTINE init_triffid
#endif
