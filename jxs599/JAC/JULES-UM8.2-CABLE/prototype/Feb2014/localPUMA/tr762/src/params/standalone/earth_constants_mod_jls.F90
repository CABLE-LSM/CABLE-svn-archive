#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Earth constants

MODULE earth_constants_mod

! Description:
!       Earth's physical constants

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

IMPLICIT NONE

! the Earth's radius
REAL, PARAMETER :: earth_radius        = 6371229.0

! g is the mean acceleration due to gravity at the Earth's surface
REAL, PARAMETER :: g                   = 9.80665


END MODULE earth_constants_mod
#endif
