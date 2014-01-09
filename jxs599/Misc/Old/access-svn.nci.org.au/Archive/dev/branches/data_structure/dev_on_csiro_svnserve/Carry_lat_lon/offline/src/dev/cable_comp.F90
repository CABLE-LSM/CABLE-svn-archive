

program cable_comp
   use air_module
   implicit none
   integer :: mp = 1          !total # of points
   integer :: kend = 17520    !toal # of timesteps

      call define_air(mp,kend)
   stop
end program cable_comp
