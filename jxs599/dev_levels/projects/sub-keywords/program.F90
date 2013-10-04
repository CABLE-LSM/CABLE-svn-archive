module m2
   implicit none
   real, allocatable, dimension(:) :: myvn

   contains

   subroutine afunc(fint,uvar)
      implicit none
      integer :: fint
      real, allocatable, dimension(:) :: uvar 
         allocate(uvar(fint))
      return
   end subroutine afunc

end module m2

!#define svnLCD '$LastChangedDate: 2012-10-12 13:24:23 +1100 (Fri, 12 Oct 2012) $' 
#define svnLCR '$Rev$' 
!#define svnLCA '$LastChangedAuthor$' 

program test_driver
   use m2
   implicit none
   integer, parameter :: fint = 3
   real, allocatable, dimension(:) :: heater 
      allocate(heater(fint))
      heater = 3440.99
      print *,svnLCR 
      !print *,svnLCA 
      !print *,'jhan: ', svnLCD 
      call afunc(fint, myvn)
      myvn = heater
   stop
end program test_driver


