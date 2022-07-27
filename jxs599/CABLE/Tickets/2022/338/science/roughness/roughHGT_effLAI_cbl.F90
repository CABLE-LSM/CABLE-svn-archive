MODULE eff_LAI_mod_cbl

IMPLICIT NONE

PUBLIC HgtAboveSnow
PUBLIC LAI_eff

CONTAINS

!Computes height above surface given that therre is snow present
!variable formerly known as rough%hruff
SUBROUTINE HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft,             &
                         SnowDepth, SnowDensity )
IMPLICIT NONE

!re-decl input args
INTEGER  :: mp
!result to return
REAL :: HeightAboveSnow(mp)

!re-decl input args
REAL :: z0surf_min
REAL :: HGT_pft(mp)
REAL :: SnowDepth(mp) !snow depth (liquid water )
REAL :: SnowDensity(mp)


!local_vars:
REAL, PARAMETER  :: fmin = 10.0 ! [meters]?
REAL, PARAMETER  :: SnowDensity_min = 100.0 ! [meters]?

REAL :: SnowDensity_eff(mp)
REAL :: HgtAboveSnow_min(mp)
REAL :: HgtAboveSnow_comp(mp)

! restrict Effective snow density to be .GE. "a" minimum
SnowDensity_eff= MAX( SnowDensity_min, SnowDensity )

! min. allowed canopy height (fixed @ 10* min. surface roughness)
HgtAboveSnow_min =  fmin * z0surf_min

!Canopy Hgt above snow given computed snow depth & PFT height
HgtAboveSnow_comp =   HGT_pft - ( 1.2 * SnowDepth / SnowDensity_eff )

! Finally Set Effective canopy height above snow level
HeightAboveSnow = MAX( HgtAboveSnow_min,  HgtAboveSnow_comp )

RETURN

END SUBROUTINE HgtAboveSnow

!Computes Effective LAI of exposed canopy given effect of snow present
!variable formerly known as canopy%vlaiw
SUBROUTINE LAI_eff( mp, lai_pft, Hgt_PFT, HgtAboveSnow,                        &
                    reducedLAIdue2snow )

  !re-decl input args
INTEGER  :: mp
REAL :: lai_pft(mp)
REAL :: Hgt_PFT(mp)
REAL :: HgtAboveSnow(mp)
REAL :: reducedLAIdue2snow(mp)

!local_vars:
REAL :: Hgt_eff(mp)
REAL :: FracOfCanopyAboveSnow(mp)

!Fraction Of Canopy Above Snow
FracOfCanopyAboveSnow = HgtAboveSnow/ MAX( 0.01, Hgt_PFT)

! LAI decreases due to snow:
reducedLAIdue2snow = lai_pft * FracOfCanopyAboveSnow

END SUBROUTINE LAI_eff

END MODULE eff_LAI_mod_cbl
