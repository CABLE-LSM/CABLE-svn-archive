#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE interpolation_mod

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
  CHARACTER(len=2), PARAMETER ::                                              &
    INTERP_AVE_BACKWARD = 'b',                                                &
       ! backward time average, i.e. time average ending at given time
       ! (in GSWP2 this is L)
    INTERP_AVE_CENTRED = 'c',                                                 &
       ! centred time average, i.e. time average centred on given time
       ! (GSWP2 C)
    INTERP_AVE_FORWARD = 'f',                                                 &
       ! forward time average, i.e. time average starting at given time
       ! (GSWP2 N)
    INTERP_INSTANT = 'i',                                                     &
       ! instantaneous value at given time (interpolation will be
       ! linear in time)(GSWP2 I)
    NO_INTERP_END = 'nb',                                                     &
       ! no interpolation, value is valid over time interval ending
       ! at given time (not in GSWP2)
    NO_INTERP_CENTRED = 'nc',                                                 &
       ! no interpolation, value is valid over time interval centred
       ! on given time  (GSWP2 0)
    NO_INTERP_START = 'nf'
       ! no interpolation, value is valid over time interval starting
       ! at given time (GSWP2 default)

! Note that interpolation for b,c and f will generate values that, generally,
! lie outside the range of the input averages


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC get_required_time_bounds, interpolate


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "get_required_time_bounds.inc"
#include "interpolate.inc"

END MODULE interpolation_mod
#endif
