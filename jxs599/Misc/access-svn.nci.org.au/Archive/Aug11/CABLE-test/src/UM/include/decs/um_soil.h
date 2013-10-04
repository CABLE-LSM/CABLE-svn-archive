      real, intent(in), dimension(um1%land_pts) :: bexp, hcon, satcon, sathh, smvcst, &
            smvcwt, smvccl, albsoil 
      real, intent(in), dimension(um1%land_pts, um1%sm_levels) :: sthu
      real, intent(in), dimension(um1%land_pts, um1%ntiles, um1%sm_levels) :: sthu_tile,     &
            sthf_tile, smcl_tile, tsoil_tile
      real, intent(in), dimension(um1%sm_levels) :: dzsoil
