MODULE cbl_LAI_canopy_height_mod

IMPLICIT NONE
PUBLIC :: limit_HGT_LAI
PRIVATE

! subroutine limit_HGT_LAI irestricts the range of canopy height and LAI as
!inherited from JULES/UM spatial maps

CONTAINS

SUBROUTINE limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, ntiles,      &
                          npft, tile_pts, tile_index, tile_frac,L_tile_pts,    &
                          LAI_pft, HGT_pft, CLAI_thresh )

IMPLICIT NONE
INTEGER, INTENT(IN) :: mp
REAL, INTENT(OUT) :: LAI_pft_cbl(mp)
REAL, INTENT(OUT) :: HGT_pft_cbl(mp)
INTEGER, INTENT(IN) :: land_pts, ntiles, npft
REAL, INTENT(IN) :: LAI_pft(land_pts, ntiles)
REAL, INTENT(IN) :: HGT_pft(land_pts, ntiles)
REAL, INTENT(IN):: tile_frac(land_pts,ntiles)
REAL, INTENT(IN) :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
INTEGER, INTENT(IN) :: tile_pts(ntiles)
INTEGER, INTENT(IN):: tile_index(land_pts,ntiles)
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts,ntiles)

!local vars
REAL :: HGT_pft_temp(land_pts,ntiles) ! needed to filter spatail map
REAL :: LAI_pft_temp(land_pts,ntiles) ! needed to filter spatail map
INTEGER :: i,j, n

! init everywhere, even where tile_frac=0
LAI_pft_temp = 0.0
HGT_pft_temp = 0.0

DO n=1,ntiles
  DO j=1,tile_pts(n)

    i = tile_index(j,n)  ! It must be landpt index

    IF ( tile_frac(i,n)  >   0.0 ) THEN
      ! LAI set either just below threshold OR from INput field
      LAI_pft_temp(i,n) = MAX(CLAI_thresh*.99,LAI_pft(i,n))
      IF (n>npft)  LAI_pft_temp(i,n) = 0.0 !to match offline Loobos
       ! hard-wired vegetation type numbers need to be removed
      IF (n < 5 ) THEN ! set trees min. height
        HGT_pft_temp(i,n) = MAX(1.0,HGT_pft(i,n))
      ELSE  ! set shrubs/grass min. height
        HGT_pft_temp(i,n) = MAX(0.1, HGT_pft(i,n))
      END IF

    END IF

  END DO
END DO

! pack filtered JULE/UM maps to CABLE variables
LAI_pft_cbl  = PACK(LAI_pft_temp, l_tile_pts)
HGT_pft_cbl  = PACK(HGT_pft_temp, l_tile_pts)

END SUBROUTINE limit_HGT_LAI

END MODULE cbl_LAI_canopy_height_mod

