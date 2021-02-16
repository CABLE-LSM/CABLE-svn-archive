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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BOUNDARY_LAYER_VAR_MOD'

CONTAINS
  
SUBROUTINE boundary_layer_var_alloc(t_i_length,t_j_length,                    &
                                    bl_levels,                                &
                                    l_deposition)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: t_i_length,t_j_length,                                 &
                       bl_levels

LOGICAL, INTENT(IN) :: l_deposition

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BOUNDARY_LAYER_VAR_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( zh(t_i_length,t_j_length))
zh(:,:) = 0.0

IF ( l_deposition ) THEN
  ALLOCATE( dzl(t_i_length,t_j_length,bl_levels))
  dzl(:,:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE boundary_layer_var_alloc

END MODULE boundary_layer_var_mod
#endif
