! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with UM setting of params for the carbon cycle.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

MODULE ccarbon

IMPLICIT NONE

! Carbon cycle and vegetation parameters
REAL, PARAMETER :: m_air   = 28.966      ! molecular weight of dry air
REAL, PARAMETER :: epco2   = 1.5194      ! Ratio molecular weights CO2/dry air
REAL, PARAMETER :: m_co2   = m_air * epco2 ! molecular weight of CO2
REAL, PARAMETER :: m_carbon= 12.0        ! molecular weight of carbon
REAL, PARAMETER :: epo2    = 1.106       ! Ratio molecular weights O2/dry air

END MODULE ccarbon
