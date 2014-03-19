! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with setting of
! minimum soil carbon

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE csmin

  IMPLICIT NONE

! This is a parameter in the UM include file. In JULES it can be assigned
! in the control routines, so we just give it a default value

! Minimum soil carbon (kg C/m2).
  REAL :: cs_min = 1.0e-6

  NAMELIST /jules_csmin/ cs_min

CONTAINS

  SUBROUTINE print_nlist_jules_csmin()
    USE jules_print_mgr, ONLY : jules_print
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer
    
    CALL jules_print('csmin', &
        'Contents of namelist jules_csmin')
    
    WRITE(lineBuffer,*)' cs_min = ',cs_min
    CALL jules_print('csmin',lineBuffer)
    
    CALL jules_print('csmin', &
        '- - - - - - end of namelist - - - - - -')
    
  END SUBROUTINE print_nlist_jules_csmin
  
END MODULE csmin
