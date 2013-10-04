
program test_driver
   implicit none
      character(len=222) :: plum
      call get_environment_variable(NAME="cmarsvn",VALUE=plum)
      print *,plum 
      call system('./wr_plum')
      !print *, 'hello again tiger'
   stop
end program test_driver

