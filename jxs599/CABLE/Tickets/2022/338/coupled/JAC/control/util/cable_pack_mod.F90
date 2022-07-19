MODULE cable_pack_mod

CONTAINS

!-----------------------------------------------------------------------------
! Description:
!   JULES met forcing vars needed by CABLE commonly have JULES dimensions
!   (row_length,rows), which are no good for CABLE. These have to be
!   re-packed in a single vector of active tiles. Hence we use
!   conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!   if the land point is/has an active tile
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
!-----------------------------------------------------------------------------
SUBROUTINE cable_pack_rr( cable_var, jules_var, mp, l_tile_pts, row_length,    &
                          rows, ntype, land_pts, land_index, surft_pts,        &
                          surft_index )
IMPLICIT NONE

REAL, INTENT(OUT)   :: cable_var(mp)
INTEGER, INTENT(IN) :: mp, row_length, rows, ntype, land_pts
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, ntype)
INTEGER, INTENT(IN) :: land_index(land_pts)           !Index in (x,y) array
INTEGER, INTENT(IN) :: surft_pts(ntype)              !# land points per PFT
INTEGER, INTENT(IN) :: surft_index(land_pts,ntype)   !Index in land_pts array
REAL, INTENT(IN)    :: jules_var(row_length, rows)
!local vars
REAL    :: fvar(land_pts, ntype)
INTEGER :: n, k, l, j, i

fvar(:, :) = 0.0
DO n = 1, ntype
  ! loop over number of points per tile
  DO k = 1, surft_pts(n)
    l = surft_index(k, n)
    j = (land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length
    fvar(l, n) = jules_var(i, j)
  END DO
END DO
cable_var =  PACK(fvar, l_tile_pts)

RETURN
END SUBROUTINE cable_pack_rr

!-----------------------------------------------------------------------------
! Description:
!   UM met forcing vars needed by CABLE which have UM dimensions
!   (land_points)[_lp], which is no good to cable. These have to be
!   re-packed in a single vector of active tiles. Hence we use
!   conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!   if the land point is/has an active tile
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
!-----------------------------------------------------------------------------
SUBROUTINE cable_pack_soil( cablevar, umvar, defaultin, mp, l_tile_pts,        &
                            ntype, land_pts, ICE_type, ICE_soilType,           &
                            SurfaceType, surft_pts, surft_index )

IMPLICIT NONE

REAL, INTENT(OUT)   :: cablevar(mp)
INTEGER, INTENT(IN) :: mp, ntype, land_pts
LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, ntype)
INTEGER, INTENT(IN) :: ICE_type, ICE_soilType
INTEGER, INTENT(IN) :: surft_pts(ntype)              !# land points per PFT
INTEGER, INTENT(IN) :: surft_index(land_pts,ntype)   !Index in land_pts array
INTEGER, INTENT(IN) :: SurfaceType(mp)
REAL, INTENT(IN)    :: umvar(land_pts)
REAL, INTENT(IN)    :: defaultin(:)
!local vars:
REAL    :: fvar(land_pts, ntype)
INTEGER :: n,k,l

fvar(:, :) = 0.0

!set all tiles to spatial UM var per land_point
DO n = 1, ntype
  DO k = 1, surft_pts(n)
    ! index of each point per tile in an array of dim=(land_pts,ntype)
    l = surft_index(k,n)
    fvar(l,n) = umvar(l)
  END DO
END DO

!PACK all to spatial UM var
cablevar = PACK(fvar, l_tile_pts)
DO n = 1, mp
  !*EXCEPT* land-ICE points use per CABLE's per PFT=ICE values
  IF (SurfaceType(n) == ICE_type)                                              &
    cablevar(n) =  defaultin(ICE_soilType)
END DO

RETURN
END SUBROUTINE cable_pack_soil

END MODULE cable_pack_mod
