#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
! ******************************** COPYRIGHT *********************************

!-----------------------------------------------------------------------------
! Module containing boundary layer-related variables.
!-----------------------------------------------------------------------------

MODULE boundary_layer_var_mod

IMPLICIT NONE

REAL, ALLOCATABLE ::                                                          &
  dzl(:,:,:),                                                                 &
    ! Separation of boundary layer levels (m).
    ! The levels are listed starting at the surface and working up.
  zh(:,:)
    ! Height above surface of top of boundary layer (m).

END MODULE boundary_layer_var_mod
#endif
