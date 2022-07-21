MODULE cable_land_albedo_mod

CONTAINS

SUBROUTINE cable_land_albedo (                                                 &
  !OUT: JULES (per rad band) albedos [GridBoxMean & per tile albedo]
  land_albedo , alb_surft,                                                     &
  !IN: JULES dimensions and associated
  row_length, rows, land_pts, nsurft, npft,                                    &
  surft_pts, surft_index, land_index,                                          &
  !IN: JULES Surface descriptions generally parametrized
  dzsoil, tile_frac, LAI_pft_um, HGT_pft_um,                                   &
  soil_alb,                                                                    &
  !IN: JULES  timestep varying fields
  cosine_zenith_angle, snow_tile,                                              &
  !IN:CABLE dimensions from grid_constants_cbl
  nsl, nsnl, nrb, nrs, mp,                                                     &
  !IN: CABLE constants
  Cz0surf_min, Clai_thresh, Ccoszen_tols, Cgauss_w, Cpi, Cpi180,               &
  !IN: CABLE Vegetation/Soil parameters. decl in params_io_cbl.F90
  VeginXfang, VeginTaul, VeginRefl, ICE_type, ICE_soiltype,                    &
  !IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90
  SoilTemp_CABLE, OneLyrSnowDensity_CABLE, SnowAge_CABLE                       &
)

!USE subroutines
! Route into CABLE & UNPACKING what we have computed for JULES
USE cable_rad_driv_mod,         ONLY: cable_rad_driver
USE cable_rad_unpack_mod,       ONLY: cable_rad_unpack

! PACKING from JULES array dims to CABLE active tile 1-D vector
USE init_active_tile_mask_mod,  ONLY: init_active_tile_mask_cbl
USE cable_pack_mod,             ONLY: cable_pack_soil,cable_pack_rr

! Define CABLE grid, sunlit/veg masks & initialize surface type params
USE def_cable_grid_mod,         ONLY: def_cable_grid
USE init_cable_parms_mod,       ONLY: init_cable_parms_rad
USE cbl_masks_mod,              ONLY: fveg_mask, fsunlit_mask,                 &
                                      fsunlit_veg_mask, L_tile_pts

!Compute canopy exposed above (potential) snow
USE cbl_LAI_canopy_height_mod,  ONLY: limit_HGT_LAI
USE cbl_hruff_mod,              ONLY: HgtAboveSnow
USE cbl_LAI_eff_mod,            ONLY: LAI_eff

IMPLICIT NONE

!-- IN: JULES model dimensions
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) :: row_length, rows               !# grid cell x, y
INTEGER, INTENT(IN) :: land_pts                       !# land points on x,y grid
INTEGER, INTENT(IN) :: nsurft, npft                   !# surface types, PFTS
INTEGER, INTENT(IN) :: nrs                            !# rad streams
                                                      !(:,:,1) direct beam VIS
                                                      !(:,:,2) diffuse visible
                                                      !(:,:,3) direct beam NIR
                                                      !(:,:,4) diffuse NIR
!-------------------------------------------------------------------------------

!--- OUT: JULES (per rad band) albedos [GridBoxMean & per tile albedo]
!-------------------------------------------------------------------------------
REAL,    INTENT(OUT) :: land_albedo(row_length,rows,nrs)    ! [land_albedo_ij]
REAL,    INTENT(OUT) :: alb_surft(Land_pts,nsurft,nrs)      ! [alb_tile]
INTEGER, INTENT(OUT) :: mp  !curr. NOT requirted OUT, however it likely will
!-------------------------------------------------------------------------------

!---IN: JULES model associated
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) :: surft_pts(nsurft)            ! # land points per PFT
INTEGER, INTENT(IN) :: surft_index(land_pts,nsurft) ! Index in land_pts array
INTEGER, INTENT(IN) :: land_index(land_pts)         ! Index in (x,y) array
!-------------------------------------------------------------------------------

!--- IN: declared in grid_cell_constants_cbl
INTEGER, INTENT(IN) :: nsl                            !# soil layers

!-- IN: JULES Surface descriptions generally parametrized
REAL, INTENT(IN) :: dzsoil(nsl)
REAL, INTENT(IN) :: tile_frac(land_pts,nsurft)      ! fraction of each surf type
REAL, INTENT(IN) :: LAI_pft_um(land_pts, npft)      ! Leaf area index.
REAL, INTENT(IN) :: HGT_pft_um(land_pts, npft)      ! Canopy height
REAL, INTENT(IN) :: soil_alb(land_pts)              ! Snow-free, soil albedo

!---IN: JULES  timestep varying fields
!-------------------------------------------------------------------------------
REAL, INTENT(IN) :: cosine_zenith_angle(row_length,rows)  ! zenith angle of sun
REAL, INTENT(IN) :: snow_tile(land_pts,nsurft)            ! snow depth (units?)
!-------------------------------------------------------------------------------

!--- IN: CABLE  declared in grid_cell_constants_cbl
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) :: nsnl                           ! max # snow layers(3)
INTEGER, INTENT(IN) :: nrb                            ! # radiation bands
!-------------------------------------------------------------------------------

!---IN: CABLE constants
!-------------------------------------------------------------------------------
REAL, INTENT(IN) :: Cz0surf_min     ! the minimum roughness of bare soil
REAL, INTENT(IN) :: Clai_thresh     ! min. LAI signalling a cell is vegetated
REAL, INTENT(IN) :: Cgauss_w(nrb)   ! Gaussian integration weights
REAL, INTENT(IN) :: Cpi             ! PI
REAL, INTENT(IN) :: Cpi180          ! PI in radians
REAL, INTENT(IN) :: Ccoszen_tols    ! sun rise/set threshold for zenith angle
                                    ! signals daylit
!---IN: CABLE Vegetation/Soil parameters. decl in params_io_cbl.F90
INTEGER, INTENT(IN) :: ICE_type, ICE_soiltype
REAL, INTENT(IN)    :: VeginXfang(nsurft)               ! Leaf Angle
REAL, INTENT(IN)    :: VeginTaul(nrb, nsurft )          ! Leaf Transmisivity
REAL, INTENT(IN)    :: VeginRefl(nrb, nsurft )          ! Leaf Reflectivity
!-------------------------------------------------------------------------------

!---IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90
REAL, INTENT(IN) :: SoilTemp_CABLE(land_pts, nsurft, nsl )
REAL, INTENT(IN) :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: SnowAge_CABLE(land_pts, nsurft )

!--- local vars (passed to subrs)

! Albedos req'd by JULES - Effective Surface Relectance as seen by atmosphere
REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)

!masks
LOGICAL, ALLOCATABLE :: veg_mask(:),  sunlit_mask(:),  sunlit_veg_mask(:)

!Vegetation/soil parameters
INTEGER, ALLOCATABLE :: SurfaceType(:)     ! veg%iveg
INTEGER, ALLOCATABLE :: SoilType(:)        ! soil%isoilm
REAL, ALLOCATABLE :: VegXfang(:)           ! Leaf Angle [veg%xfang]
REAL, ALLOCATABLE :: VegTaul(:,:)          ! Leaf Transmisivity [veg%taul]
REAL, ALLOCATABLE :: VegRefl(:,:)          ! Leaf Reflectivity [veg%refl]

REAL, ALLOCATABLE :: reducedLAIdue2snow(:) ! Eff. LAI IF snow [canopy%vlaiw]
REAL, ALLOCATABLE :: HeightAboveSnow(:)    ! Canopy Hgt above snow (rough%hruff)

! arrays to map IN progs to CABLE vector length
REAL, ALLOCATABLE :: SnowDepth(:)     ! Total Snow depth - water eqivalent -
                                      ! ssnow%snowd
REAL, ALLOCATABLE :: SnowDensity(:)   ! Total Snow density (assumes 1 layer
REAL, ALLOCATABLE :: SoilTemp(:)      ! Soil Temperature of top layer (soil%tgg)
REAL, ALLOCATABLE :: SnowAge(:)       ! Snow age (assumes 1 layer describes snow
                                      ! ssnow%isflag

REAL, ALLOCATABLE :: AlbSoil(:,:)     ! BareSoil Albedo [soill%albsoil]
REAL, ALLOCATABLE :: coszen(:)
!computed from UM HT(LAI)_PFT passed in explicit call - need at rad call
REAL, ALLOCATABLE  :: LAI_pft_cbl(:)           !Formerly: ~veg%vlai
REAL, ALLOCATABLE  :: HGT_pft_cbl(:)           !Formerly:  ~veg%hc

LOGICAL :: um_online = .FALSE.
LOGICAL :: jls_standalone = .TRUE.
LOGICAL :: jls_radiation = .TRUE.     !um_radiation = jls_radiation

INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: subr_name = "cable_rad_main"
! End header

! initialise INTENT(OUT) fields
land_albedo = 0.0; alb_surft = 0.0

! Define mapping mask. i.e. l_tile_pts =TRUE (active) , where tile_frac > 0
CALL init_active_tile_mask_cbl(l_tile_pts, land_pts, nsurft, tile_frac )

! CABLE works  on 1D vector (length=mp) of active points (ONLY)
CALL def_cable_grid(mp, nsurft, surft_pts)

! alloc/zero each timestep
CALL alloc_local_vars( EffSurfRefl_dif, EffSurfRefl_beam,                      &
                       reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,           &
                       HeightAboveSnow, coszen, mp, nrb )

CALL alloc_local_progs( SnowDepth,             SnowDensity,                    &
                        SoilTemp,             SnowAge, AlbSoil, mp, nrb )

! -----------------------------------------------------------------------------
! PACK CABLE fields
! -----------------------------------------------------------------------------
! map PFT/soil parameters to mp format
CALL init_cable_parms_rad( VegXfang, VegTaul, VegRefl, SurfaceType, SoilType,  &
                           mp, nrb, l_tile_pts, ICE_Type, ICE_SoilType,        &
                           VeginXfang, VeginTaul, VeginRefl,                   &
                           land_pts, nsurft, tile_frac )

!JULES doesn't have albedos per nrb: assume same: i.e (:,2) = (:,1) & (:,3)=0
albsoil(:,3)  = 0.0
CALL cable_pack_soil(  albsoil(:,1), soil_alb, soil_alb, mp, l_tile_pts,       &
                       nsurft, land_pts, ICE_type, ICE_soilType,               &
                       SurfaceType, surft_pts, surft_index )

CALL cable_pack_soil(  albsoil(:,2), soil_alb, soil_alb, mp, l_tile_pts,       &
                       nsurft, land_pts, ICE_type, ICE_soilType,               &
                       SurfaceType, surft_pts, surft_index )

CALL cable_pack_rr(    coszen, cosine_zenith_angle,                            &
                       mp, l_tile_pts, row_length, rows, nsurft, land_pts,     &
                       land_index, surft_pts, surft_index )

CALL cable_pack_progs( SnowDepth,             SnowDensity,                     &
                       SoilTemp,           SnowAge, mp, land_pts, nsurft, nsl, &
                       nsnl, l_tile_pts,             snow_tile,                &
                       OneLyrSnowDensity_CABLE,                                &
                       SoilTemp_CABLE,                 SnowAge_CABLE )
! -----------------------------------------------------------------------------

! limit IN height, LAI  and initialize existing cable % types
CALL limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, nsurft,            &
                    surft_pts, surft_index, tile_frac, l_tile_pts,             &
                    LAI_pft_um, HGT_pft_um, CLAI_thresh )

!set Height of Canopy above snow (rough%hruff)
CALL HgtAboveSnow( HeightAboveSnow, mp, Cz0surf_min, HGT_pft_cbl,              &
                   SnowDepth, SnowDensity )

!set Effective LAI considering ipotentail snow coverage
CALL LAI_eff( mp, LAI_pft_cbl, HGT_pft_cbl, HeightAboveSnow,                   &
                reducedLAIdue2snow)

!Define logical masks for vegetated,sunlit, both
CALL fveg_mask( veg_mask, mp, Clai_thresh, reducedLAIdue2snow )
CALL fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, coszen )
CALL fsunlit_veg_mask( sunlit_veg_mask, veg_mask, sunlit_mask, mp )

!------------------------------------------------------------------------------
! Call CABLE_rad_driver to run specific and necessary components of CABLE
!------------------------------------------------------------------------------
CALL cable_rad_driver( EffSurfRefl_dif, EffSurfRefl_beam,                      &
                       mp, nrb, Clai_thresh, Ccoszen_tols,                     &
                       CGauss_w, Cpi, Cpi180, Cz0surf_min,                     &
                       veg_mask, sunlit_mask, sunlit_veg_mask,                 &
                       jls_standalone, jls_radiation, SurfaceType,             &
                       SoilType, LAI_pft_cbl, HGT_pft_cbl,                     &
                       SnowDepth, SnowDensity, SoilTemp, SnowAge,              &
                       AlbSoil ,coszen, VegXfang, VegTaul, VegRefl,            &
                       HeightAboveSnow, reducedLAIdue2snow )

! Unpack variables (CABLE computed albedos) to JULES
!------------------------------------------------------------------------------
CALL cable_rad_unpack( land_albedo, alb_surft, mp, nrb, row_length, rows,      &
                       land_pts, nsurft, nsl, surft_pts, surft_index,          &
                       land_index, tile_frac, l_tile_pts,                      &
                       EffSurfRefl_dif, EffSurfRefl_beam  )

CALL flush_local( EffSurfRefl_dif, EffSurfRefl_beam,                           &
                        SnowDepth,             SnowDensity,                    &
                        SoilTemp,           SnowAge, AlbSoil,                  &
                        reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,          &
                        HeightAboveSnow, coszen )

!flick switches before leaving
jls_radiation= .FALSE.

RETURN

END SUBROUTINE cable_land_albedo

!==============================================================================

SUBROUTINE alloc_local_vars( EffSurfRefl_dif, EffSurfRefl_beam,                &
                             reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,     &
                             HeightAboveSnow, coszen, mp, nrb )
IMPLICIT NONE

INTEGER :: mp, nrb
REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)
REAL, ALLOCATABLE :: reducedLAIdue2snow(:) !make also again in cbl_model_driver
REAL, ALLOCATABLE :: HeightAboveSnow(:)  ! Canopy hgt above snow (rough%hruff)
REAL, ALLOCATABLE :: LAI_pft_cbl(:)      ! Formerly: ~veg%vlai
REAL, ALLOCATABLE :: HGT_pft_cbl(:)      ! Formerly: ~veg%hc
REAL, ALLOCATABLE :: coszen(:)

IF ( .NOT. ALLOCATED(reducedLAIdue2snow)) ALLOCATE( reducedLAIdue2snow(mp) )
IF ( .NOT. ALLOCATED(HeightAboveSnow) )   ALLOCATE( HeightAboveSnow(mp) )
IF ( .NOT. ALLOCATED(LAI_pft_cbl) )       ALLOCATE( LAI_pft_cbl(mp) )
IF ( .NOT. ALLOCATED(HGT_pft_cbl) )       ALLOCATE( HGT_pft_cbl(mp) )
IF (.NOT. ALLOCATED(coszen) )            ALLOCATE (coszen(mp) )
IF (.NOT. ALLOCATED(EffSurfRefl_dif) )   ALLOCATE (EffSurfRefl_dif(mp, nrb) )
IF (.NOT. ALLOCATED(EffSurfRefl_beam) )  ALLOCATE (EffSurfRefl_beam(mp, nrb) )

RETURN

END SUBROUTINE alloc_local_vars

!==============================================================================

SUBROUTINE alloc_local_progs( SnowDepth,             SnowDensity,              &
                              SoilTemp,           SnowAge, AlbSoil, mp, nrb )

IMPLICIT NONE

INTEGER :: mp, nrb
!local to CABLE and can be flushed every timestep
REAL, ALLOCATABLE :: SnowDepth(:)        ! Total Snow depth - water eqivalent
REAL, ALLOCATABLE :: SnowDensity(:)      ! Total Snow density (assumes 1 layer
REAL, ALLOCATABLE :: SoilTemp(:)         ! Soil Temp. of top layer (soil%tgg)
REAL, ALLOCATABLE :: SnowAge( :)         ! assumes 1 layer describes snow
REAL, ALLOCATABLE :: AlbSoil(:,:)        ! bare soil albedo parametrized

IF ( .NOT. ALLOCATED(SnowDepth) )         ALLOCATE( SnowDepth(mp) )
IF ( .NOT. ALLOCATED(SnowDensity) )       ALLOCATE( SnowDensity(mp) )
IF ( .NOT. ALLOCATED(SoilTemp) )          ALLOCATE( SoilTemp(mp) )
IF ( .NOT. ALLOCATED(SnowAge) )           ALLOCATE( SnowAge(mp) )
IF (.NOT. ALLOCATED(AlbSoil) )           ALLOCATE (AlbSoil(mp, nrb) )

AlbSoil(:,:) = 0.0

RETURN

END SUBROUTINE alloc_local_progs

!==============================================================================

SUBROUTINE cable_pack_progs( SnowDepth,             SnowDensity,               &
                       SoilTemp,           SnowAge, mp, land_pts, nsurft, nsl, &
                       nsnl, l_tile_pts,             snow_tile,                &
                       OneLyrSnowDensity_CABLE,                                &
                       SoilTemp_CABLE,                 SnowAge_CABLE )

IMPLICIT NONE

INTEGER, INTENT(IN)  :: land_pts, nsurft, nsl, nsnl,mp

! map IN progs to CABLE veector length
REAL, INTENT(OUT) :: SnowDepth(mp)     ! Total Snow depth - water eqivalent -
REAL, INTENT(OUT) :: SnowDensity(mp)   ! Total Snow density (assumes 1 layer
REAL, INTENT(OUT) :: SoilTemp(mp)      ! Soil Temperature of top layer (soil%tgg)
REAL, INTENT(OUT) :: SnowAge(mp)      ! Snow age (assumes 1 layer describes snow

LOGICAL, INTENT(IN) :: l_tile_pts(land_pts, nsurft )

!---IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90
REAL, INTENT(IN) :: SoilTemp_CABLE(land_pts, nsurft, nsl )
REAL, INTENT(IN) :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: SnowAge_CABLE(land_pts, nsurft )
REAL, INTENT(IN) :: snow_tile(land_pts,nsurft)            ! snow depth (units?)

!Store Snow Depth from previous timestep. Treat differently on 1st timestep
SnowDepth   = PACK( snow_tile, l_tile_pts )
SnowDensity = PACK( OneLyrSnowDensity_CABLE, l_tile_pts )

!Surface skin/top layer Soil/Snow temperature
SoilTemp =   PACK( SoilTemp_CABLE(:,:,1), l_tile_pts )
SnowAge  =   PACK( SnowAge_CABLE(:,:), l_tile_pts )

RETURN

END SUBROUTINE cable_pack_progs

!==============================================================================

SUBROUTINE flush_local( EffSurfRefl_dif, EffSurfRefl_beam,                     &
                        SnowDepth,             SnowDensity,                    &
                        SoilTemp,           SnowAge, AlbSoil,                  &
                        reducedLAIdue2snow, LAI_pft_cbl, HGT_pft_cbl,          &
                        HeightAboveSnow, coszen )

IMPLICIT NONE

INTEGER, ALLOCATABLE :: SnowFlag_3L(:)   ! treat snow as 3 layer - if enough
REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)
REAL, ALLOCATABLE :: SnowDepth(:)        ! Total Snow depth - water eqivalent
REAL, ALLOCATABLE :: SnowODepth(:)       ! Total Snow depth before any update
REAL, ALLOCATABLE :: SnowDensity(:)      ! Total Snow density (assumes 1 layer
REAL, ALLOCATABLE :: SoilTemp(:)         ! Soil Temp. of top layer (soil%tgg)
REAL, ALLOCATABLE :: SnowTemp(:)         ! Snow Temp. of top layer (soil%tss)
REAL, ALLOCATABLE :: SnowAge( :)         ! assumes 1 layer describes snow
REAL, ALLOCATABLE :: AlbSoil(:,:)        ! bare soil albedo parametrized
REAL, ALLOCATABLE :: reducedLAIdue2snow(:) !make also again in cbl_model_driver
REAL, ALLOCATABLE :: HeightAboveSnow(:)  ! Canopy hgt above snow (rough%hruff)
REAL, ALLOCATABLE :: LAI_pft_cbl(:)      ! Formerly: ~veg%vlai
REAL, ALLOCATABLE :: HGT_pft_cbl(:)      ! Formerly: ~veg%hc
REAL, ALLOCATABLE :: coszen(:)

IF ( ALLOCATED (EffSurfRefl_dif)  ) DEALLOCATE (EffSurfRefl_dif )
IF ( ALLOCATED (EffSurfRefl_beam) ) DEALLOCATE (EffSurfRefl_beam )
IF ( ALLOCATED(SnowDepth) )         DEALLOCATE( SnowDepth )
IF ( ALLOCATED(SnowODepth) )        DEALLOCATE( SnowODepth )
IF ( ALLOCATED(SnowDensity) )       DEALLOCATE( SnowDensity )
IF ( ALLOCATED(SnowFlag_3L) )       DEALLOCATE( SnowFlag_3L )
IF ( ALLOCATED(SoilTemp) )          DEALLOCATE( SoilTemp )
IF ( ALLOCATED(SnowTemp) )          DEALLOCATE( SnowTemp )
IF ( ALLOCATED(SnowAge) )           DEALLOCATE( SnowAge )
IF ( ALLOCATED(AlbSoil) )           DEALLOCATE( AlbSoil )
IF ( ALLOCATED(reducedLAIdue2snow)) DEALLOCATE( reducedLAIdue2snow)
IF ( ALLOCATED(HeightAboveSnow) )   DEALLOCATE( HeightAboveSnow )
IF ( ALLOCATED(LAI_pft_cbl) )       DEALLOCATE( LAI_pft_cbl )
IF ( ALLOCATED(HGT_pft_cbl) )       DEALLOCATE( HGT_pft_cbl )
IF ( ALLOCATED(coszen) )            DEALLOCATE( coszen )

RETURN
END SUBROUTINE flush_local

!==============================================================================

END MODULE cable_land_albedo_mod

