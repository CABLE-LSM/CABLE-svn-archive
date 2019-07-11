integer :: row_length                       !grid cell x
integer :: rows                             !grid cell y
integer :: land_pts                         !grid cell land points on the x,y grid
integer :: ntiles                           !grid cell number of surface types 
integer :: sm_levels                        !grid cell number of soil levels 
integer :: npft                             !grid cell number of PFTs 
real :: dzsoil(sm_levels)                   !soil thicknesses in each layer  
real :: tile_frac(land_pts,ntiles)          !fraction of each surface type per land point 
integer :: tile_pts(ntiles)                 !Number of land points per PFT 
integer :: tile_index(land_pts,ntiles)      !Index of land point in (land_pts) array
integer :: land_index(land_pts)             !Index of land points in (x,y) array - see below
real :: LAI_pft_um(land_pts, npft)          !Leaf area index.
real :: HGT_pft_um(land_pts, npft)          !Canopy height

real :: snow_tile(land_pts,ntiles)          ! snow depth equivalent (in water?)
real :: land_albedo(row_length,rows,4)      
real :: alb_surft(Land_pts,ntiles,4)        
real :: soil_alb(land_pts)                  !(albsoil)Snow-free, bare soil albedo

real :: soil_temp_cable( land_pts, ntiles, sm_levels )   !Mean snow density
real :: snow_temp_cable(land_pts, ntiles, 3)             !Snow temperature (3 layer)
real :: snow_flag_cable(land_pts,ntiles)                 !surface type fraction 
real :: snow_avg_rho_cable(LAND_PTS,NTILES)              !snow cover in the ice tile.

real :: surf_down_sw(row_length,rows,4)                  !4-band ShortWave forcing
real :: cosine_zenith_angle(row_length,rows)             !cosine_zenith_angle          

!Convoluted mapping using land_index(array) to get back to the row_length,rows co-ordinate
! J = ( LAND_INDEX(L)-1 ) / ROW_LENGTH + 1
! I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
! FTL_1(I,J) = where ftl_1(row_length,rows) 
