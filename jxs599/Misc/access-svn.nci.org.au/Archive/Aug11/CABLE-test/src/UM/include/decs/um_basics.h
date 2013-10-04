         integer:: row_length, rows, land_pts, ntiles, npft, &
                     sm_levels, timestep 
         integer, dimension(:) :: tile_pts, land_index
         integer, dimension(:,:) :: tile_index
         real, dimension(:,:) :: tile_frac
         real, dimension(:,:) :: latitude, longitude 
         !___true IF vegetation (tile) fraction is greater than 0
         logical,dimension(:,:) :: l_tile_pts
         real, intent(in) :: rho_water 
