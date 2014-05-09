
module cat_nodes_mod 
   implicit none
   private
   public :: cat_nodes

   contains

   subroutine cat_nodes() 
      use cat_common                                                            
      use cat_read_mod
      use cat_write_mod
      implicit none                                                               
      integer :: i
      integer, save :: ndimx=0, pndimx=0
      character(len=30) :: nfilename, chnodes
      real, dimension(:), allocatable :: ar_x
         
         ! array  holding the "sub-mp" per node file
         allocate( dimx_i(i_nodes) )
 
         ! loop over descriptive text.dat node files
         do i=1, i_nodes
            
            write(chnodes,10), i-1 
  10        format(I3.3)   
            nfilename = trim( trim(filename)//trim(chnodes) )

            ! read descriptive .dat node files to set up dimensions
            call read_txt_ifile( trim( nfilename) )

         end do   

         ! loop over node files
         do i=1, i_nodes

            write(chnodes,10), i-1 
            nfilename = trim( trim(filename)//trim(chnodes) )
            pndimx = ndimx + pndimx 

            ! read descriptive .dat node files to set per node file dimensions
            call read_txt_file( trim( nfilename) , ndimx )
!print *, ""
!print *, "filename ", trim( nfilename)
!print *, "dimx ", ndimx 

            ! read binary (data) node files 
            call read_dat_file( trim( nfilename ), pndimx+1, pndimx+ndimx )
         end do   
         
!print *, ""
!print *, "filename ", trim( filename)
!print *, "dimx ", pndimx + ndimx 

         call rewrite_txt_file( trim(filename), pndimx+ndimx )
         call rewrite_dat_file( trim(filename) )

      return
   end subroutine cat_nodes 

end module cat_nodes_mod 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


