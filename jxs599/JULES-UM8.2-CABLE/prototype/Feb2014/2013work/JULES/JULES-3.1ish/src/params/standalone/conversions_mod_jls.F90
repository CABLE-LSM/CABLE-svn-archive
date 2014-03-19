#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Global standard conversions 

MODULE conversions_mod

! Description:
!       Model and code section invariant physical constants

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

! Stripped down version for standalone code

IMPLICIT NONE

! 
REAL, PARAMETER :: pi                  = 3.14159265358979323846
      
!Conversion factor degrees to radians   
REAL, PARAMETER :: pi_over_180         = pi/180.0  
!Conversion factor radians to degrees                          
REAL, PARAMETER :: recip_pi_over_180   = 180.0/pi


! zerodegc is a conversion between degrees centigrade and kelvin
REAL, PARAMETER :: zerodegc            = 273.15


END MODULE conversions_mod
#endif
