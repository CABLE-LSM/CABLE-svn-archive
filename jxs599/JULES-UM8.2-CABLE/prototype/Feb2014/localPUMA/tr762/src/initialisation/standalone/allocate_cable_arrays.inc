

  ! CABLE-JULES vars to be initialized from file
  ALLOCATE( snow_depth_CABLE(land_pts, ntiles, nsmax), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( soil_temp_CABLE(land_pts, ntiles, sm_levels), STAT=error ) 
  error_sum = error_sum + error

  ALLOCATE( smcl_CABLE(land_pts, ntiles, sm_levels), STAT=error ) 
  error_sum = error_sum + error

  ALLOCATE( sthf_CABLE(land_pts, ntiles, sm_levels), STAT=error ) 
  error_sum = error_sum + error

  ALLOCATE( SNOW_MASS_CABLE(land_pts, ntiles, nsmax), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( snow_temp_CABLE(land_pts, ntiles, nsmax), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( SNOW_RHO_CABLE(land_pts, ntiles, nsmax), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( SNOW_FLG3L_CABLE(land_pts, ntiles), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( SNOW_RHO1L_CABLE(land_pts, ntiles), STAT=error )
  error_sum = error_sum + error

  ALLOCATE( SNOW_AGE_CABLE(land_pts, ntiles), STAT=error )
  error_sum = error_sum + error

