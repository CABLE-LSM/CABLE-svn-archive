
program test_driver
   implicit none
      character(len=222) :: plum
      call get_environment_variable(NAME="cmarsvn",VALUE=plum)
      print *,plum 
      !call system('/short/p66/jxs599/UM_ROUTDIR/xafsy/script')
      !call system('/short/p66/jxs599/UM_ROUTDIR/xafsy/xafsy.exe')
      !print *, 'hello again tiger'
   stop
end program test_driver

