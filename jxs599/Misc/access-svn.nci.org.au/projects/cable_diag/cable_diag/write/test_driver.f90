
program test_driver
   use cable_diag_mod 
   implicit none
      integer,parameter :: n =100
      integer ,parameter :: timesteps =10
      integer ,parameter :: node = 0
      integer :: i, this_timestep 
      real :: ar(n)

      do i=1,n
         ar(i) = real(i )      
      enddo        

      do i=1,timesteps
         this_timestep = i     
         call cable_diag( 1, 'test', n, timesteps, this_timestep, node, &
                                    'test', ar )      
      enddo        

   stop
end program test_driver

