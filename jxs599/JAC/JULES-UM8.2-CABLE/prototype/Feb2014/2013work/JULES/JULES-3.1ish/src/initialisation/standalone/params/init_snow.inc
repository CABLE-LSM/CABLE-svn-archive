#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/params/init_snow.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_snow()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE water_constants_mod, ONLY : tm

  USE ancil_info, ONLY : nsmax

  USE switches, ONLY : l_aggregate, can_model

  USE nstypes, ONLY : npft

  USE rad_param

  USE snow_param

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises snow parameters
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

! Open the snow parameters namelist file
  OPEN(NAMELIST_UNIT, FILE='snow_params.nml',                                 &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_snow",                                               &
                   "Error opening namelist file snow_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! There are two namelists to read from this file
  CALL log_info("init_snow", "Reading JULES_SNOW_PARAM namelist...")
  READ(NAMELIST_UNIT, nml=jules_snow_param, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_snow",                                               &
                   "Error reading namelist JULES_SNOW_PARAM " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_snow", "Reading JULES_RAD_PARAM namelist...")
  READ(NAMELIST_UNIT, nml=jules_rad_param, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_snow",                                               &
                   "Error reading namelist JULES_RAD_PARAM " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_snow",                                               &
                   "Error closing namelist file snow_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values and set derived variables
!-----------------------------------------------------------------------------
! Initialise dzsnow by copying from fixed length IO counterpart
  dzsnow(:) = dzsnow_io(1:nsmax)

! Copy PFT values into cansnowtile
! Set switch for canopy snow model. This can be TRUE only at PFT tiles.
  canSnowTile(:) = .FALSE.
  IF ( .NOT. l_aggregate .AND. can_model==4 ) THEN
    canSnowTile(1:npft) = cansnowpft(1:npft)
  ENDIF

! MOSES code divides kland by dtland. Note this blows up if dtland=0.
  kland = kland_numerator / dtland
  tcland = tm - dtland

  RETURN

END SUBROUTINE init_snow
#endif
