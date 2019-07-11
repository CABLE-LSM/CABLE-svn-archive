#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/var/latlon_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

MODULE latlon_mod

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
  REAL, ALLOCATABLE :: latitude(:,:)  ! The latitude of model points
  REAL, ALLOCATABLE :: longitude(:,:) ! The longitude of model points

END MODULE latlon_mod
#endif
