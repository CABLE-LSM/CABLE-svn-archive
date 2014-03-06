

  ! CABLE-JULES vars to be initialized from file
  ALLOCATE( snow_depth_CABLE(land_pts, ntiles, nsmax), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( soil_temp_CABLE(land_pts, ntiles, sm_levels), STAT=error ) 
  error_sum = error_sum + error

