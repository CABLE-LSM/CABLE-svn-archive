
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module cat_nodes_mod 
   implicit none
   private
   public :: cat_nodes

   contains

   subroutine cat_nodes() 
      use debug_common                                                            
      use debug_read_mod
      use debug_write_mod
      implicit none                                                               
      integer :: i
      integer, save :: ndimx=0, pndimx=0
      character(len=30) :: nfilename, chnodes
      real(gs), dimension(:), allocatable :: ar_x
        allocate( dimx_i(n_nodes) )
        do i=1, n_nodes
            write(chnodes,10), i-1 
  10        format(I2.2)   
            nfilename = trim( trim(filename)//trim(chnodes) )
            call read_txt_ifile( trim( nfilename) )
         end do   
         do i=1, n_nodes
            write(chnodes,10), i-1 
            nfilename = trim( trim(filename)//trim(chnodes) )
            pndimx = ndimx + pndimx 
            call read_txt_file( trim( nfilename) , ndimx )
            call read_dat_file( trim( nfilename ), pndimx+1, pndimx+ndimx )
         end do   
         call rewrite_txt_file( trim(filename), pndimx+ndimx )
         call rewrite_dat_file( trim(filename) )
      return
   end subroutine cat_nodes 

end module cat_nodes_mod 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


