      real, intent(inout), dimension(um1%row_length, um1%rows) :: sw_down
      real, intent(in), dimension(um1%row_length, um1%rows) :: lw_down, &
                                        sin_theta_latitude
      real, intent(inout), dimension(um1%row_length, um1%rows) :: cos_zenith_angle
      real, intent(in), dimension(um1%row_length, um1%rows, 4) :: surf_down_sw 
           
