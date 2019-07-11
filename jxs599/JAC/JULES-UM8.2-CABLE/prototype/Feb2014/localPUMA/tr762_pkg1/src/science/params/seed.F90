! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module sets the values of the variables FRAC_MIN and FRAC_SEED


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE seed

  IMPLICIT NONE

! These are set as parameters in UM include file - we just assign them
! default values.

! Minimum areal fraction for PFTs.
  REAL :: frac_min  = 1.0e-6

! "Seed" fraction for PFTs.
  REAL :: frac_seed = 0.01

  NAMELIST /jules_seed/ frac_min,frac_seed

CONTAINS

  SUBROUTINE print_nlist_jules_seed()
    USE jules_print_mgr, ONLY : jules_print
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer

    CALL jules_print('seed', &
        'Contents of namelist jules_seed')

    WRITE(lineBuffer,*)' frac_min = ',frac_min
    CALL jules_print('seed',lineBuffer)
    WRITE(lineBuffer,*)' frac_seed = ',frac_seed
    CALL jules_print('seed',lineBuffer)

    CALL jules_print('seed', &
        '- - - - - - end of namelist - - - - - -')

  END SUBROUTINE print_nlist_jules_seed

END MODULE seed
