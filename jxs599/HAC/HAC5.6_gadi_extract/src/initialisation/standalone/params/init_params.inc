#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_params(nml_dir)

!Module imports
USE lsm_switch_mod,           ONLY: lsm_id, jules, cable



!Common modules
USE ereport_mod,              ONLY: ereport


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises various model parameters
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists


!-----------------------------------------------------------------------------
! Variables

INTEGER :: errorstatus

!-----------------------------------------------------------------------------

! Process the PFT and non-veg parameters

SELECT CASE(lsm_id)

CASE (jules)
  CALL init_pftparm_jules(nml_dir)
  CALL init_nvegparm_jules(nml_dir)

CASE (cable)
  CALL init_pftparm_cable(nml_dir)
  CALL init_soilparm_cable(nml_dir)

CASE DEFAULT
  errorstatus = 101
  CALL ereport('init_params', errorstatus,                                    &
               'Unrecognised surface scheme')
    
END SELECT
  
! Process the crop parameters
CALL init_cropparm(nml_dir)

! Process TRIFFID parameters
CALL init_triffid(nml_dir)

! Process atmospheric deposition species.
CALL init_deposition_species(nml_dir)

RETURN

END SUBROUTINE init_params
#endif
