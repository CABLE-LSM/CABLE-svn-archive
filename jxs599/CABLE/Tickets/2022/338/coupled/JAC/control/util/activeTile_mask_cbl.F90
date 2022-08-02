MODULE init_active_tile_mask_mod
IMPLICIT NONE
PUBLIC :: init_active_tile_mask_cbl
PRIVATE

CONTAINS

SUBROUTINE init_active_tile_mask_cbl(l_tile_pts, land_pts, nsurft, tile_frac )
!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT(OUT), ALLOCATABLE :: L_tile_pts(:,:)
INTEGER, INTENT(IN) :: land_pts, nsurft
REAL,    INTENT(IN) :: tile_frac(land_pts,nsurft)        !fraction of each surf type

!Local vars:
INTEGER :: i, j

! Determine active tiles map
IF ( .NOT. ALLOCATED(l_tile_pts)) ALLOCATE( l_tile_pts(land_pts, nsurft) )

l_tile_pts(:,:) = .FALSE.

DO j = 1, nsurft
  DO i = 1, land_pts
    IF ( tile_frac(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO

RETURN

END SUBROUTINE init_active_tile_mask_cbl

END MODULE init_active_tile_mask_mod
