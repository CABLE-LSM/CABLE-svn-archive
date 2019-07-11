#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Global atmospheric physical constants

MODULE atmos_constants_mod

! Description:
!       Global atmospheric physical constants

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.


! This is a stripped down version for the JULES standalone code

IMPLICIT NONE

! Von Karman's constant
REAL, PARAMETER :: vkman               = 0.4

! r, the Gas constant for dry air
REAL, PARAMETER :: r                   = 287.05
! cp, specific heat of dry air at constant pressure
REAL, PARAMETER :: cp                  = 1005.0

REAL, PARAMETER :: kappa               = r/cp

! pref, reference surface pressure 
REAL, PARAMETER :: pref                = 100000.0
REAL, PARAMETER :: p_zero              = pref

!scale height h
REAL, PARAMETER :: sclht               = 6.8E+03


! repsilon, ratio of molecular weights of water and dry air
REAL, PARAMETER :: repsilon            = 0.62198
REAL, PARAMETER :: c_virtual           = 1.0/repsilon-1.0

END MODULE atmos_constants_mod
#endif
