
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
 
         print *, "Reading N files to allocate super-array:" 
         ! loop over descriptive text.dat node files
         do i=1, i_nodes
            
            write(chnodes,10), i-1 
  10        format(I3.3)   
            nfilename = trim( trim(filename)//trim(chnodes) )

            ! read descriptive .dat node files to set up dimensions
            call read_txt_ifile( trim( nfilename) )

         end do   
         print *, "Allocated super-array" 

         ! loop over node files
         do i=1, i_nodes

            write(chnodes,10), i-1 
            nfilename = trim( trim(filename)//trim(chnodes) )
            pndimx = ndimx + pndimx 

            ! read descriptive .dat node files to set per node file dimensions
            call read_txt_file( trim( nfilename) , ndimx )

            !print *, "Text file read for file: ", nfilename 
            ! read binary (data) node files 
            !print *, "Data file: ",  nfilename
            !print *, "x limits: ",  pndimx+1, pndimx+ndimx 
            call read_dat_file( trim( nfilename ), pndimx+1, pndimx+ndimx )
            print *, "Data read for file: ",  nfilename
         end do   
         
         print *, "Writing to a singlefile:" 
         call rewrite_txt_file( trim(filename), pndimx+ndimx )
         call rewrite_dat_file( trim(filename) )

      return
   end subroutine cat_nodes 

end module cat_nodes_mod 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


