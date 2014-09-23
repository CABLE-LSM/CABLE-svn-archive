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

! Use the include file if we are doing a UM run, else define variables

#if defined(UM_JULES)

#include "c_topog.h"

#else

!-----------------------------------------------------------------------
! Scalar parameters.
!-----------------------------------------------------------------------
! Topographic index increment:
  REAL,PARAMETER :: dti = 0.2
  
! Standard deviation of LOG(Ksat(0)):
  REAL,PARAMETER :: sigma_logk = 0.0
  
!-----------------------------------------------------------------------
! Scalar variables.
!-----------------------------------------------------------------------
! Maximum topographic index considered:
  REAL :: ti_max = 10.0
  
! Parameter to remove very high water tables
! from the calculated wetland fraction:
  REAL :: ti_wetl = 2.0
  
! Maximum allowed water table depth (m):
  REAL :: zw_max = 15.0

#endif

END MODULE c_topog
