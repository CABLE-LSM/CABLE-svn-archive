!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

!Variables to be calculated and returned by CABLE
!------------------------------------------------------------------------------
land_albedo,   & ! GridBoxMean albedo per rad band (row_length,rows,4) [land_albedo_ij]
alb_surft,     & ! albedo per rad band per tile (land_pts, ntiles, 4) [alb_tile] 
!------------------------------------------------------------------------------

!This is an "outlier" and possibly misleading. CABLE does not actually calculate radiation:
!Generally we speak of 4-band radiation. This is actually only 2-bands VIS/NIR in the SW.
!We further split each of these into Direct Beam and Diffuse components. 
!Offline CABLE splits the total SW forcing into VIS/NIR using a Spitter() fuction
!We include this variable here to connect back to JULES toplevel routines because:
!Online the UM radiation scheme DOES compute surf_down_sw using a more sophisticated model 
!than that which we use Offline, however not until AFTER the surface albedos have been 
!calculated which IS what is done AND here and technically does not require knowledge of 
!the downward SW. JULES aggregates this SW and threads this to the LSM as it is called 
!explicitly. 
!------------------------------------------------------------------------------
surf_down_sw,     & ! ShortWave radiation per rad band (row_length,rows,4) 
!------------------------------------------------------------------------------

!Mostly model dimensions and associated
!------------------------------------------------------------------------------
row_length,          & !grid cell x
rows,                & !grid cell y
land_pts,            & !grid cell land points on the x,y grid
ntiles,              & !grid cell number of surface types [nsurft] 
sm_levels,           & !grid cell number of soil levels 
npft,                & !grid cell number of PFTs 
tile_pts,            & !Number of land points per PFT [surft_pts] 
tile_index,          & !Index of land point in (land_pts) array[surft_index] 
land_index,          & !Index of land points in (x,y) array - see  corresponding *decs.inc
!------------------------------------------------------------------------------

!Surface descriptions generally parametrized
!------------------------------------------------------------------------------
dzsoil,              & !soil thicknesses in each layer  
tile_frac,           & !fraction of each surface type per land point [frac_surft] 
LAI_pft_um,          & !Leaf area index. [LAI_pft]
HGT_pft_um,          & !Canopy height [canht_pft]
soil_alb,            & !(albsoil)Snow-free, bare soil albedo [albsoil_soilt(:,1)]
!------------------------------------------------------------------------------

!Variables passed from JULES/UM
!------------------------------------------------------------------------------
snow_tile,           & !snow depth equivalent (in water?) [snow_surft]
                       !This is the total snow depth per tile. CABLE also has depth per layer
cosine_zenith_angle, & ! cosine_zenith_angle [cosz_ij]
!------------------------------------------------------------------------------

!CABLE prognostics
!------------------------------------------------------------------------------
soil_temp_cable,     & ! Soil temperature (3 layer) [SoilTemp_cable] 
snow_temp_cable,     & ! Snow temperature (3 layer) [SnowTemp_cable]
snow_flag_cable,     & ! flag indicating whether enough snow to treat as 3 layers
                       ! [ThreeLayerSnowFlag_CABLE]
snow_avg_rho_cable   & ! density considering snow as 1 layer [OneLyrSnowDensity_CABLE] 
!------------------------------------------------------------------------------
