!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

!Variables to be calculated and returned by CABLE
!------------------------------------------------------------------------------
FTL_TILE,           & !surface sensible heat flux  [W/m2]? (up=+ve?) -"FTL_tile" in CABLE
FQW_TILE,           & !surface moisture flux flux  [kg/m^2/s/m2](up=+ve?) units could be changed?
TSTAR_TILE,         & !surface temperature [K] per tile
U_S,                & ! land point surface friction velocity [m/s]
U_S_STD_TILE,       & ! per tile surface friction velocity [m/s] canopy%us
CD_TILE,            &
CH_TILE,            &
RADNET_TILE,        & !Net radiation at surface [W/m2]
FRACA,              & !Fraction of surface moisture flux with only aerodynamic resistance for
                      !snow-free land tiles.
RESFS,              & !Combined soil, stomatal and aerodynamic resistance factor for fraction
                      !(1-FRACA) of snow-free land tiles.
RESFT,              & !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for snow.
Z0H_tile,           & ! Tile roughness lengths for heat and moisture (m).
Z0M_tile,           & ! OUT Tile roughness lengths for momentum.
RECIP_L_MO_tile,    & ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
EPOT_tile,          & ! Potential evaporation from surface, per tile
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
row_length,         & !grid cell x
rows,               & !grid cell y
land_pts,           & !grid cell land points on the x,y grid
ntiles,             & !grid cell number of surface types [] 
sm_levels,          & !grid cell number of soil levels 
npft,               & !grid cell number of PFTs 
tile_pts,           & !Number of land points per PFT [] 
tile_index,         & !Index of land point in (land_pts) array[] 
land_index,         & !Index of land points in (x,y) array - see  corresponding *decs.inc
timestep_width,     & !bin width in seconds of timestep
endstep,            & !last timestep of experiment
timestep_number,    &
doy,                &  
mp,                 &
nrb,                &
!------------------------------------------------------------------------------

!Surface descriptions generally parametrized
!------------------------------------------------------------------------------
Fland,              & !fraction of land per land point (could be coastal) 
tile_frac,          & !fraction of each surface type per land point [frac_surft] 
L_tile_pts,         & !Logical mask TRUE where tile frac > 0. used to PACK/UNPACK
LAI_ft,             & !Leaf area index. [LAI_pft/LAI_pft_um in radiation]
canht_ft,           & !Canopy height [canht_pft/HGT_pft_um in radiation]
albsoil,            & !(albsoil)Snow-free, bare soil albedo [albsoil_soilt(:,1) in um ]
z0surf_min,         &
dzsoil,             & !soil thicknesses in each layer  
bexp,               &
hcon,               &
satcon,             &
sathh,              &
smvcst,             &
smvcwt,             &
smvccl,             &
!------------------------------------------------------------------------------

!Variables passed from JULES/UM
!------------------------------------------------------------------------------
latitude,           & !latitude
longitude,          & !longitude
cosine_zenith_angle,& ! cosine_zenith_angle [cosz_ij]
sin_theta_latitude, &
sthu,               &
lw_down,            &
ls_rain,            &
ls_snow,            &
tl_1,               &
qw_1,               &
vshr_land,          &
pstar,              &
z1_tq,              &
z1_uv,              &
snow_tile,          & !snow depth equivalent (in water?) [snow_surft]
                      !This is the total snow depth per tile. CABLE also has depth per layer
canopy_tile,        &
CO2_MMR,            &
!------------------------------------------------------------------------------

!CABLE prognostics
!------------------------------------------------------------------------------
iThreeLayerSnowFlag,& ! flag indicating whether enough snow to treat as 3 layers
                      ! [real(ThreeLayerSnowFlag_CABLE)]
OneLyrSnowDensity,  & ! density considering snow as 1 layer [OneLyrSnowDensity_CABLE 
SnowAge,            & 
SnowDensity,        & 
SnowDepth,          &
SnowMass,           & 
SnowTemp,           & 
SoilMoisture,       & 
FrozenSoilFrac,     & 
SoilTemp            & 
!------------------------------------------------------------------------------




