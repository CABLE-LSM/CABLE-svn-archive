#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing constants/parameters used in the dust scheme
!
MODULE dust_parameters_mod

!
! Description:
!   This module contains declarations for constants and tunable 
!   parameters used to diagnose the emission and deposition of 
!   mineral dust
!  
! Owner: 
!   See Unified Model Code Owner's HTML page
!   This file belongs in section 17: CLASSIC Aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards version 8.2.
!

IMPLICIT NONE

!
! Parameters fundamental to the dust scheme 
! (Details of dust properties/size divisions etc.)
! =======================================================================
!
! Number of discrete particle size divisions (bins)
INTEGER, PARAMETER :: ndiv=6 ! number of divisions that can be lifted 
                             !  from the surface
INTEGER, PARAMETER :: ndivh=9! number of divisions that can be blown 
                             !  horizontally along the surface and contribute 
                             !  to the lifting of the 1st NDIV divisions
!
! The size of particles included in each division 
! Note that by using two arrays for the max/min diameter here we can set
! up overlapping divisions. We must take care to make them consistent

REAL, PARAMETER ::  dmax(ndiv) =                                        &
       (/ 2.0E-7, 6.32456E-7, 2.0E-6, 6.32456E-6, 2.0E-5, 6.32456E-5 /)
                ! max diameter of particles in each div.
REAL, PARAMETER ::  dmin(ndiv) =                                        &
       (/ 6.32456E-8, 2.0E-7, 6.32456E-7, 2.0E-6, 6.32456E-6, 2.0E-5 /)
                ! min diameter of particles in each div.
REAL, PARAMETER ::  drep(ndiv) =                                        &
                     (/ 0.112468E-06, 0.355656E-06, 0.112468E-05,       &   
                        0.355656E-05, 0.112468E-04, 0.355656E-04 /)
                ! representative particle diameter
!
! Physical properties of the dust
REAL, PARAMETER :: rhop = 2.65E+3  ! density of a dust particle (quartz)

!
! Parameters used during dust emissions calculations
! =======================================================================
!
! Parameters based on observations/published research
REAL, PARAMETER, DIMENSION(ndivh) :: ustd_bas =                       &
       (/ 0.85, 0.72, 0.59, 0.46, 0.33, 0.16, 0.14, 0.18, 0.28 /)
                                      ! impact U*t derived from Bagnold (1941)
REAL, PARAMETER :: horiz_c = 2.61     ! C in horizontal flux calc (White 1979)
REAL, PARAMETER :: vert_a = 13.4      ! A in vertical flux calc(Gillette 1979)
REAL, PARAMETER :: vert_b = -6.       ! B in vertical flux calc(Gillette 1979)
REAL, PARAMETER :: vert_c = 0.01      ! cgs to si conversion

!
! Input variables needed when using the dust lifting scheme with 1 tile
!
REAL            :: z0_soil            ! Roughness length over bare soil

!
! Tuning parameters. Note that these can be redefined in the dust_params 
!                    namelist below
REAL            :: us_aa = 0.0        ! ustar correction (additional)
REAL            :: us_am = 1.0        ! ustar correction (multiplic.)
REAL            :: ust_aa = 0.0       ! ustar_t correction (add)
REAL            :: ust_am = 1.0       ! ustar_t correction (multi.)
REAL            :: sm_corr = 1.0      ! soil moist. correction factor
REAL            :: horiz_d = 1.0      ! Global tuning param for horizontal
                                      !  (and hence vertical) flux

!
! Limits used in dust flux calculations
REAL            :: u_s_min = 1.e-5    ! Minimum val of u_s_std,
                                      !  below which no dust is produced
                                      !  (avoids divide by zero problem)
REAL            :: clay_max = 0.1     ! Max clay fraction.
REAL            :: snowmin = 1.e-6    ! Min snow depth for dust
REAL            :: h_orog_limit = 150.! 1/2pk-trough height above which
                                      !  no dust is produced
REAL            :: fland_lim =0.99999 ! No ems if fland<lim as windspeed
                                      !  too high at coastal points
!
! Switch to diagnose vertical flux using a fixed (user definable) size
!  distribution. If set to false, the vertical flux in each bin at a point 
!  is poportional to the horizontal flux in that bin at that point.
LOGICAL         :: l_fix_size_dist = .false.
REAL            :: size_dist(ndiv) =                                    &
                    (/ 0.0005, 0.0049, 0.0299, 0.2329, 0.4839, 0.2479 /)
                                      ! The proportion of flux from each 
                                      !  bin if using the fixed
                                      !  size distribution

!
! Parameters used during the gravitational settling of dust
! =======================================================================
!
REAL, PARAMETER :: accf = 1.257 ! Cunningham correction factor term A
REAL, PARAMETER :: bccf = 0.4   ! Cunningham correction factor term B
REAL, PARAMETER :: cccf = -1.1  ! Cunningham correction factor term C

LOGICAL         :: l_cam_hotspot_fix = .false.
! Dust hotspot fix previously used with the CAM dust scheme. 
!  This allows tracer to settle directly from level 1 rather 
!  than treating this process in parallel with BL mixing in 
!  bl_trmix_dd, which is a more scientifically valid approach. 
!  Dust hotspots have now been fixed by performing
!  gravitational settling before BL mixing, so this should no 
!  longer be necessary. This option will probably be removed in the
!  near future.

!
! Parameters used during the scavenging of dust
! =======================================================================
!
REAL, PARAMETER :: krain_dust(ndiv) =                                   &
       (/ 2.e-5, 2.e-5, 3.e-5, 6.e-5, 4.e-4, 4.e-4 /)
                                       ! scav. coeff. for rain
REAL, PARAMETER :: ksnow_dust(ndiv) =                                   &
       (/ 2.e-5, 2.e-5, 3.e-5, 6.e-5, 4.e-4, 4.e-4 /)
                                       ! scav. coeff. for snow

!
! RUN_Dust namelist via which non-parameter values herein can be set
! =======================================================================
! 
  NAMELIST /RUN_Dust/                                                   &
       us_aa, us_am, ust_aa, ust_am, sm_corr, horiz_d,                  &
       clay_max, fland_lim,                                             &
       l_fix_size_dist, size_dist,                                      &
       l_cam_hotspot_fix

  CONTAINS
!
! Internal subroutine to check that entries to RUN_dust are consistent
! =======================================================================
! 
  SUBROUTINE dust_parameters_check( )
!
!   Check that a user-defined emissions size distribution is normalised
!
    IF (l_fix_size_dist) THEN
      size_dist(:)=size_dist(:)/SUM(size_dist)
    END IF
    
  END SUBROUTINE dust_parameters_check

END MODULE dust_parameters_mod

#endif
