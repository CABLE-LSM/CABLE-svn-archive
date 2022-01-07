MODULE cable_photo_constants_mod
IMPLICIT NONE

PUBLIC

INTEGER, PARAMETER :: maxiter      = 20      ! max # interations for leaf temperature

REAL, PARAMETER :: a1c4_default    = 4.0
REAL, PARAMETER :: a1c3_default    = 9.0      ! a1c3 is defined incable_canopy.f90
! Empirical coef for vpd sensitivity of stomata
REAL, PARAMETER :: d0c3_hawkesbury = 5000.0
REAL, PARAMETER :: d0c3_default    = 1500.0
REAL, PARAMETER :: d0c4_default    = 1500.0
REAL, PARAMETER :: alpha3          = 0.2

REAL, PARAMETER :: alpha4          = 0.05
REAL, PARAMETER :: cfrd3           = 0.015
REAL, PARAMETER :: cfrd4           = 0.025
REAL, PARAMETER :: conkc0          = 302.0E-6 ! mol mol^-1
REAL, PARAMETER :: conko0          = 256.0E-3 ! mol mol^-1
REAL, PARAMETER :: convx3          = 0.01
REAL, PARAMETER :: convx4          = 0.8
REAL, PARAMETER :: ekc             = 59430.0  ! J mol^-1
REAL, PARAMETER :: eko             = 36000.0  ! J mol^-1
REAL, PARAMETER :: gam0            = 28.0E-6  ! mol mol^-1 @ 20C = 36.9 @ 25C
REAL, PARAMETER :: gam1            = 0.0509
REAL, PARAMETER :: gam2            = 0.0010
REAL, PARAMETER :: gsw03           = 0.01
REAL, PARAMETER :: gsw04           = 0.04
REAL, PARAMETER :: rgbwc           = 1.32
REAL, PARAMETER :: rgswc           = 1.57
REAL, PARAMETER :: tmaxj           = 45.0
REAL, PARAMETER :: tmaxv           = 45.0
REAL, PARAMETER :: tminj           = -5.0
REAL, PARAMETER :: tminv           = -5.0
REAL, PARAMETER :: toptj           = 20.0
REAL, PARAMETER :: toptv           = 20.0
REAL, PARAMETER :: trefk           = 298.2    ! reference temperature K

END module cable_photo_constants_mod

