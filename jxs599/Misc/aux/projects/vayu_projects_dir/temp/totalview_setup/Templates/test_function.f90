


program test_driver
   use generic_pars
   implicit none
   real :: var
   var = test_fn()
   stop
end program test_driver


   real function test_fn() 
      use generic_pars
      implicit none
      real :: val
         test_fn = val
      return 
   subroutine test_fn 

