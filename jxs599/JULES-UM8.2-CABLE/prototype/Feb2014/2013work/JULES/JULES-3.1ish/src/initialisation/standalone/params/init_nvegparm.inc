#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/params/init_nvegparm.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_nvegparm()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE switches_urban, ONLY : l_urban2t

  USE nstypes, ONLY : nnvg, npft, urban_canyon, urban_roof

  USE c_z0h_z0m, ONLY : z0h_z0m

  USE nvegparm

  USE nvegparm_io

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the non-vegetation parameters
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


  CALL log_info("init_nvegparm", "Reading JULES_NVEGPARM namelist...")

! Open the pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE='nveg_params.nml',                                 &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Error opening namelist file nveg_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_nvegparm, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Error reading namelist JULES_NVEGPARM " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Error closing namelist file nveg_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays
  nvgname(:)    = nvgname_io(1:nnvg)
  albsnc_nvg(:) = albsnc_nvg_io(1:nnvg)
  albsnf_nvg(:) = albsnf_nvg_io(1:nnvg)
  catch_nvg(:)  = catch_nvg_io(1:nnvg)
  gs_nvg(:)     = gs_nvg_io(1:nnvg)
  infil_nvg(:)  = infil_nvg_io(1:nnvg)
  z0_nvg(:)     = z0_nvg_io(1:nnvg)
  ch_nvg(:)     = ch_nvg_io(1:nnvg)
  vf_nvg(:)     = vf_nvg_io(1:nnvg)
  emis_nvg(:)   = emis_nvg_io(1:nnvg)

  z0h_z0m(npft+1:) = z0hm_nvg_io(1:nnvg)

  RETURN

END SUBROUTINE init_nvegparm
#endif
