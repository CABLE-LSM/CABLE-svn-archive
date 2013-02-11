      !___IN: UM dimensions, array indexes, flags
      integer, intent(in) :: row_length, rows, land_pts, ntiles, npft, sm_levels
      integer, intent(in), dimension(land_pts) :: land_index 
      integer, intent(in), dimension(ntiles) :: tile_pts 
      integer, intent(in), dimension(land_pts, ntiles) :: tile_index ,isnow_flg3l 

      !___UM parameters 
      integer, intent(in) :: itimestep
      real, intent(in) :: rho_water 
      real, intent(in), dimension(sm_levels) :: dzsoil

      !___UM soil/snow/radiation/met vars
      real, intent(in), dimension(land_pts) :: bexp, hcon, satcon, sathh, smvcst,&
            smvcwt, smvccl, albsoil, fland 
      real, intent(inout), dimension(row_length,rows) :: sw_down, &
               cos_zenith_angle
      real, intent(in), dimension(row_length,rows) :: latitude, longitude,       &
            lw_down, ls_rain, ls_snow, tl_1, qw_1,    &
            vshr_land, pstar, z1_tq, z1_uv
      real, intent(inout), dimension(land_pts, ntiles) :: snow_tile
      real, intent(in), dimension(land_pts, ntiles) :: tile_frac, &   
            snow_rho1l, snage_tile
      real, intent(in), dimension(row_length, rows, 4) :: surf_down_sw 
      real, intent(in), dimension(land_pts, npft) :: canht_ft, lai_ft 
      real, intent(in),dimension(land_pts, ntiles) :: canopy_tile
      real, intent(inout), dimension(land_pts, ntiles,3) :: snow_cond
      real, intent(in), dimension(land_pts, ntiles,3) :: &
            snow_rho3l, snow_depth3l, snow_mass3l, snow_tmp3l
      real, intent(in), dimension(land_pts, sm_levels) :: sthu 
      real, intent(in), dimension(land_pts, ntiles, sm_levels) :: sthu_tile,     &
            sthf_tile,smcl_tile, tsoil_tile
      real, intent(in) :: co2_mmr

      !___true IF vegetation (tile) fraction is greater than 0
      logical, intent(inout),dimension(land_pts, ntiles) :: L_tile_pts
     
      real :: sin_theta_latitude(row_length,rows) 
