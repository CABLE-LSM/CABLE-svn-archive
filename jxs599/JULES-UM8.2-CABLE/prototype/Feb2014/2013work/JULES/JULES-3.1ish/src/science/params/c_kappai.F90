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
    
MODULE c_kappai

IMPLICIT NONE

#include "c_kappai.h"

! Define any variables not in UM include file c_kappai.h

! Effective thickness of sea surface layer (m).
  REAL,PARAMETER:: dzsea = 1.0


END MODULE c_kappai
