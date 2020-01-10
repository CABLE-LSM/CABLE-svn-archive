! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of
!

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE c_topog

IMPLICIT NONE

!-----------------------------------------------------------------------
! Scalar parameters.
!-----------------------------------------------------------------------
! Topographic index increment:
REAL,PARAMETER :: dti = 0.2

! Standard deviation of LOG(Ksat(0)):
REAL,PARAMETER :: sigma_logk = 0.0

END MODULE c_topog
