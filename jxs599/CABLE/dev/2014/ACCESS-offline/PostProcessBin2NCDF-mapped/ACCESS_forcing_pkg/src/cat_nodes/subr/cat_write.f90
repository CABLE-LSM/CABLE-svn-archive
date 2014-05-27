

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module cat_write_mod
   implicit none
   contains

   !==========================================================================!
   !==========================================================================!
   
   subroutine rewrite_txt_file( Lfilename, Ldimx )
      use cat_common
      implicit none
      character(len=*), intent(in) :: Lfilename
      integer, intent(in) :: Ldimx
      integer, parameter :: gok=0
      integer :: gopenstatus
      integer :: i,j,k
         print *, "Writing single file: ",  Lfilename
         open(unit=1,file=Lfilename//'.dat', status="unknown",action="write",  &
              iostat=gopenstatus )
            
            if(gopenstatus==gok) then
               write (1,*), 'Number of var(s): '
               write (1,*) Nvars
               write (1,*), 'Name of var(s): '
   7139        format(a)            
               do i=1,Nvars
                  write (1,7139), trim(ar_varname(i))
               enddo
               write (1,*), 'dimension of var(s) in x: '
               write (1,*), Ldimx 
               write (1,*), 'dimension of var(s) in y: '
               write (1,*), dimy 
            else
               write (*,*), Lfilename//'.dat',' NOT found to write'
            endif
         
         close(1)
         print *, "Written single file: ",  Lfilename
      
      return 
   end subroutine rewrite_txt_file 

   !==========================================================================!
   !==========================================================================!

   subroutine rewrite_dat_file( Lfilename )
      use cat_common
      implicit none
      character(len=*), intent(in) :: Lfilename
      integer, parameter :: gok=0
      integer :: gopenstatus
      integer :: i,j
      real, dimension(:,:), allocatable :: Lar
!      integer :: frecl
!      frecl = Nvars * dimx*r_1
      
      print *, "Writing field: "
      
      allocate( Lar(dimy, Nvars*dimx_tot) ) 
      print *, "Allocated typesetting array: "

         do i=1,dimy
            do j=1,Nvars
               Lar(i, (j-1)*dimx_tot+1 : j*dimx_tot )  = ar_data(j,i,:)
               !print *, "typeset field for timestep: ", i
           enddo   
         enddo   
         print *, "Filled typesetting array: "

         !open(unit=2,file=Lfilename//'.bin',status="unknown",action="write", &
         open(unit=2,file='/home/599/jxs599/test.bin',status="unknown",action="write", &
                  iostat=gopenstatus, form="unformatted" )

            !if(gopenstatus==gok) then
               do i=1,dimy
                  write (2) Lar(i,:) 
                  !print *, "Written field for timestep: ", i
               enddo
            !else
               write (*,*), Lfilename//'.bin',' NOT found for write'
            !endif
         close(2)
      return 
   end subroutine rewrite_dat_file 
   !==========================================================================!
   !==========================================================================!



end module cat_write_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
