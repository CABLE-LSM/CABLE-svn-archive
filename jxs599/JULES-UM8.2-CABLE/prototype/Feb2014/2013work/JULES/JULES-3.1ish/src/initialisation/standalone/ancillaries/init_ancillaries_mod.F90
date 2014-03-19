#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/ancillaries/init_ancillaries_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE init_ancillaries_mod

  USE input_mod, ONLY : fill_variables_from_file

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "init_ancillaries.inc"
#include "init_frac.inc"
#include "init_soil.inc"
#include "init_top.inc"
#include "init_agric.inc"

#include "calc_fit_fsat.inc"

! CABLE
#include "init_cable.inc"

END MODULE init_ancillaries_mod

#endif
