      integer, intent(in), dimension(um1%land_pts, um1%ntiles) :: isnow_flg3l 
      real, intent(inout), dimension(um1%land_pts, um1%ntiles) :: snow_tile
      real, intent(in), dimension(um1%land_pts, um1%ntiles) :: snow_rho1l, &
                                                  snage_tile
      real, intent(inout), dimension(um1%land_pts, um1%ntiles,3) :: snow_cond
      real, intent(in), dimension(um1%land_pts, um1%ntiles,3) :: snow_rho3l, &
                                                   snow_depth3l, snow_mass3l,  &
                                                   snow_tmp3l
      real, intent(in), dimension(um1%land_pts) :: fland 
