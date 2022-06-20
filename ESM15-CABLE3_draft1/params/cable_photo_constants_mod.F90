!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Defines constants for CABLE
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Combines cable_*_constants from earlier versions
!          Will include define_types in future version.
!
!
! ==============================================================================

module cable_photo_constants_mod
   implicit none

   public

   ! definition of major types of constants

      integer:: maxiter=20 ! max # interations for leaf temperature
      real :: gam0 = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
      real :: gam1 = 0.0509
      real :: gam2 = 0.0010
      real :: rgbwc  = 1.32
      real :: rgswc  = 1.57
      real :: tmaxj  = 45.0
      real :: tmaxv  = 45.0
      real :: tminj  = -5.0
      real :: tminv  = -5.0
      real :: toptj  = 20.0
      real :: toptv  = 20.0
      real :: trefk= 298.2  !reference temperature K

!ESM15 uses - replaced formulation in trunk
      real :: gsw03  = 0.01
      real :: gsw04  = 0.04
      real :: conkc0 = 302.e-6  !mol mol^-1
      real :: conko0 = 256.e-3  !mol mol^-1
      real :: ekc = 59430.0  !J mol^-1
      real :: eko = 36000.0  !J mol^-1
      real :: d0c3 = 1500.0
      real :: d0c4 = 1500.0
      real :: a1c3 = 9.0
      real :: a1c4 = 4.0
      real :: cfrd3  = 0.010
      real :: cfrd4  = 0.025
      
      real :: alpha3 = 0.200
      real :: alpha4  = 0.05
      real :: convx3 = 1.0E-2
      real :: convx4 = 0.8
End module cable_photo_constants_mod

