#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_params(nml_dir)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises various model parameters
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists


!-----------------------------------------------------------------------------


! Process the PFT parameters
  CALL init_pftparm(nml_dir)

! Process the non-veg parameters
  CALL init_nvegparm(nml_dir)

! Process TRIFFID parameters
  CALL init_triffid(nml_dir)

! Process snow parameters
  CALL init_snow(nml_dir)

! Process miscellaneous other parameters
  CALL init_misc(nml_dir)

  RETURN

END SUBROUTINE init_params
#endif
