#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/update/update_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE update_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
  INTEGER ::                                                                  &
    io_precip_type,                                                           &
      !  Flag indicating how precipitation is input
      !      1 = total precipitation is read in
      !      2 = values for total rainfall and total snowfall are read in
      !      3 = values for large-scale rainfall, convective rainfall and
      !          total snowfall are read in
      !      4 = values for convective rainfall, large-scale rainfall,
      !          convective snowfall and large-scale snowfall are read in
    io_rad_type
      !  Flag indicating how radiation is input
      !      1 downward fluxes provided
      !      2 net (all wavelength) longwave flux and downward shortwave flux
      !        are provided
      !      3 net downward fluxes are provided

  LOGICAL ::                                                                  &
    io_wind_speed,                                                            &
      !   T means that the windspeed is input
      !   F means 2 components of wind are input
    use_diff_rad
      !   T means diffuse radiation is input
      !   F means diffuse radiation is set to a constant fraction (diffFracConst)

  REAL ::                                                                     &
    t_for_snow      = 274.0,                                                  &
      !   air temperature (K) at or below which precipitation is assumed to
      !   be snow
    t_for_con_rain  = 373.15,                                                 &
      !   air temperature (K) at or above which rainfall is taken to
      !   be convective (rather than large-scale) in nature
    diff_frac_const = 0.0
      !   a constant value for fraction of radiation that is diffuse

  LOGICAL ::                                                                  &
    have_prescribed_veg = .FALSE.
      !   T - a vegetation variable has been prescribed that requires a call
      !       to sparm after update
      !   F - no such variables have been prescribed


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "update_derived_variables.inc"

END MODULE update_mod
#endif
