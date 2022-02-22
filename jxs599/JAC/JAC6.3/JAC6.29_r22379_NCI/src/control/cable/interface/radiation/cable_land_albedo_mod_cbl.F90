MODULE cable_land_albedo_mod 
  
CONTAINS

SUBROUTINE cable_land_albedo (                                                 &
  !OUT: JULES (per rad band) albedos [GridBoxMean & per tile albedo]
  land_albedo , alb_surft,                                                     &
  !IN: JULES dimensions and associated
  row_length, rows, land_pts, nsurft, npft,                                    &
  tile_pts, tile_index, land_index,                                            &
  !IN: JULES Surface descriptions generally parametrized
  dzsoil, tile_frac, LAI_pft_um, HGT_pft_um,                                   &
  soil_alb,                                                                    &
  !IN: JULES  timestep varying fields
  cosine_zenith_angle, snow_tile,                                              &
  !IN:CABLE dimensions from grid_constants_cbl 
  nsl, nsnl, nrb,                                                              &
  !IN: CABLE constants
  Cz0surf_min, Clai_thresh, Ccoszen_tols, Cgauss_w, Cpi, Cpi180,               &
  !IN: CABLE Vegetation/Soil parameters. decl in params_io_cbl.F90 
  VegXfang, VegTaul, VegRefl,                                                  &
  !IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90 
  SoilTemp_CABLE, SnowTemp_CABLE,                                              &
  OneLyrSnowDensity_CABLE                                                      &
) 
 
IMPLICIT NONE

!-- IN: JULES model dimensions
INTEGER, INTENT(IN) :: row_length, rows               !# grid cell x, y
INTEGER, INTENT(IN) :: land_pts                       !# land points on x,y grid
INTEGER, INTENT(IN) :: nsurft, npft                   !# surface types, PFTS  

!--- OUT: JULES (per rad band) albedos [GridBoxMean & per tile albedo]
REAL, INTENT(OUT) :: land_albedo(row_length,rows,4)    ! [land_albedo_ij] 
REAL, INTENT(OUT) :: alb_surft(Land_pts,nsurft,4)      ! [alb_tile] 

!---IN: JULES model associated
INTEGER, INTENT(IN) :: tile_pts(nsurft)               !# land points per PFT 
INTEGER, INTENT(IN) :: tile_index(land_pts,nsurft)    !Index in land_pts array
INTEGER, INTENT(IN) :: land_index(land_pts)           !Index in (x,y) array

!--- IN: declared in grid_cell_constants_cbl
INTEGER, INTENT(IN) :: nsl                            !# soil layers

!-- IN: JULES Surface descriptions generally parametrized
REAL, INTENT(IN) :: dzsoil(nsl)
REAL, INTENT(IN) :: tile_frac(land_pts,nsurft)        !fraction of each surf type
REAL, INTENT(IN) :: LAI_pft_um(land_pts, npft)        !Leaf area index.
REAL, INTENT(IN) :: HGT_pft_um(land_pts, npft)        !Canopy height
REAL, INTENT(IN) :: soil_alb(land_pts)                !Snow-free, bare soil albedo

!---IN: JULES  timestep varying fields
REAL, INTENT(IN) :: cosine_zenith_angle(row_length,rows)  !zenith angle of sun
REAL, INTENT(IN) :: snow_tile(land_pts,nsurft)            ! snow depth (units?)

!--- IN: CABLE  declared in grid_cell_constants_cbl
INTEGER, INTENT(IN) :: nsnl                           !max # snow layers(3)
INTEGER, INTENT(IN) :: nrb                            !# radiation bands

!---IN: CABLE constants
REAL, INTENT(IN) :: Cz0surf_min     !the minimum roughness of bare soil
REAL, INTENT(IN) :: Clai_thresh     !min. LAI signalling a cell is vegetated
REAL, INTENT(IN) :: Cgauss_w(nrb)   !Gaussian integration weights
REAL, INTENT(IN) :: Cpi             !PI - describing the ratio of circumference to diameter
REAL, INTENT(IN) :: Cpi180          !PI in radians
REAL, INTENT(IN) :: Ccoszen_tols    !sun rise/set threshold for zenith angle 
                                    !signals daylit
!---IN: CABLE Vegetation/Soil parameters. decl in params_io_cbl.F90 
REAL, INTENT(IN) :: VegXfang(nsurft)                !Leaf Angle
REAL, INTENT(IN) :: VegTaul(nsurft, nrb)            !Leaf Transmisivity 
REAL, INTENT(IN) :: VegRefl(nsurft, nrb)            !Leaf Reflectivity 

!---IN: CABLE prognostics. decl in progs_cbl_vars_mod.F90 
REAL, INTENT(IN) :: SoilTemp_CABLE(land_pts, nsurft, nsl )
REAL, INTENT(IN) :: SnowTemp_CABLE(land_pts, nsurft, nsnl)
REAL, INTENT(IN) :: OneLyrSnowDensity_CABLE(land_pts, nsurft )

! End header

! initialise INTENT(OUT) fields for now until CABLE is implemented
land_albedo = 0.0 
alb_surft = 0.0 

write(6,*) "Currently CABLE rad/albedo @6.3 not implemented"

RETURN

END SUBROUTINE cable_land_albedo

END MODULE cable_land_albedo_mod 

