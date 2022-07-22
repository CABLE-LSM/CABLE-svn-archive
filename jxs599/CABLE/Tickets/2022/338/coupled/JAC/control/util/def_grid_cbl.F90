MODULE def_cable_grid_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: def_cable_grid

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='def_cable_grid_mod'

CONTAINS

!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------
SUBROUTINE def_cable_grid(mp, nsurft, surft_pts)

IMPLICIT NONE

INTEGER, INTENT(IN)   :: nsurft             !# surface types
INTEGER, INTENT(IN)   :: surft_pts(nsurft)  !# land points per PFT
INTEGER, INTENT(OUT)  :: mp

! Determine the number of active tiles
mp = SUM(surft_pts)

RETURN
END SUBROUTINE def_cable_grid

END MODULE def_cable_grid_mod
