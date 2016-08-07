
module cable_rrm_nc_names
    implicit none
    public
    character(len=*), parameter :: mask_name     = "land_mask"
    character(len=*), parameter :: length_name   = "river_distance"
    character(len=*), parameter :: slope_name    = "stddev_elevation"  !"slope"
    character(len=*), parameter :: elev_name     = "outlet_elevation"   !"elevation"
    character(len=*), parameter :: rdir_name     = "river_direction"
    character(len=*), parameter :: src_area_name = "source_area"
    character(len=*), parameter :: rr_lat_dim_name = "lat"
    character(len=*), parameter :: rr_lon_dim_name = "lon"
    character(len=*), parameter :: rr_lat_var_name = "latitude"
    character(len=*), parameter :: rr_lon_var_name = "longitude"
    
end module cable_rrm_nc_names
