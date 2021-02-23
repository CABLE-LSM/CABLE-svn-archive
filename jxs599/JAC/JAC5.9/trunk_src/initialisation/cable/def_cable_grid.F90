MODULE def_cable_grid_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: def_cable_grid

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='def_cable_grid_mod'

CONTAINS

SUBROUTINE def_cable_grid(mp, surft_pts, frac_surft)

USE init_active_tile_mask_mod,        ONLY: init_active_tile_mask_cbl
USE ancil_info,           ONLY: nsurft, land_pts

IMPLICIT NONE

INTEGER :: mp
integer:: surft_pts(nsurft)  
REAL:: frac_surft(land_pts,nsurft)  
!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Determine the number of active tiles
mp = SUM(surft_pts)

CALL init_active_tile_mask_cbl(frac_surft)

END SUBROUTINE def_cable_grid

END MODULE def_cable_grid_mod
