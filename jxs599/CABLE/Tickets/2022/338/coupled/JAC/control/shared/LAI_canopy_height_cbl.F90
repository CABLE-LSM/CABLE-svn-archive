MODULE cbl_LAI_canopy_height_mod

CONTAINS

SUBROUTINE limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, ntiles,      &
                          tile_pts, tile_index, tile_frac,L_tile_pts,          &
                          LAI_pft, HGT_pft, CLAI_thresh )

IMPLICIT NONE
REAL :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
INTEGER :: land_pts, ntiles
INTEGER :: tile_pts(ntiles)
INTEGER:: tile_index(land_pts,ntiles)
REAL:: tile_frac(land_pts,ntiles)
LOGICAL :: L_tile_pts(land_pts,ntiles)
INTEGER :: mp
REAL :: LAI_pft(land_pts, ntiles)
REAL :: HGT_pft(land_pts, ntiles)

!return vars
REAL :: LAI_pft_cbl(mp)
REAL :: HGT_pft_cbl(mp)

!local vars
REAL :: HGT_pft_temp(land_pts,ntiles)
REAL :: LAI_pft_temp(land_pts,ntiles)
INTEGER :: i,j, n

!Retain init where tile_frac=0
LAI_pft_temp = 0.0
HGT_pft_temp = 0.0

DO n=1,ntiles
  DO j=1,tile_pts(n)

    i = tile_index(j,n)  ! It must be landpt index

    IF ( tile_frac(i,n)  >   0.0 ) THEN

      LAI_pft_temp(i,n) = MAX(CLAI_thresh*.99,LAI_pft(i,n))
      IF (n>13)  LAI_pft_temp(i,n) = 0.0 !to match offline Loobos
       ! hard-wired vegetation type numbers need to be removed
      IF (n < 5 ) THEN ! trees
        HGT_pft_temp(i,n) = MAX(1.0,HGT_pft(i,n))
      ELSE ! shrubs/grass
        HGT_pft_temp(i,n) = MAX(0.1, HGT_pft(i,n))
      END IF

    END IF

  END DO
END DO

  !surface_type = PACK(surface_type_temp, um1%L_TILE_PTS)
LAI_pft_cbl  = PACK(LAI_pft_temp, l_tile_pts)
HGT_pft_cbl  = PACK(HGT_pft_temp, l_tile_pts)

END SUBROUTINE limit_HGT_LAI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cbl_LAI_canopy_height_mod

