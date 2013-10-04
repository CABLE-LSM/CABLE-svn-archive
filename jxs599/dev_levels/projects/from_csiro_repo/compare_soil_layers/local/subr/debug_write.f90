
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_write_mod
   use debug_common
   implicit none
   contains

   !==========================================================================!
   !==========================================================================!
   subroutine write_txt_file( filename )
      implicit none
      character(len=*), intent(in) :: filename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      integer :: i,j,k
         open(unit=1,file=filename//'.txt', status="unknown",action="write", iostat=gopenstatus )
            if(gopenstatus==gok) then
               do i=1,Nvars
                  write(1,*), ar_varname(i)
                  write(1,*),'dimensions:  ', dimx, dimy
                  do j=1, dimy
                     do k=1, dimx
                        write(1,*), ar_data(i,j,k) 
                     enddo 
                     write(1,*)
                  enddo 
                  write(1,*), 'END ', ar_varname(i) 
                  write(1,*)
               enddo
            else
               write (*,*), filename//'.txt',' : unable to write'
            endif
         close(1)
      return 
   end subroutine write_txt_file 

   !==========================================================================!
   !==========================================================================!

end module debug_write_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
