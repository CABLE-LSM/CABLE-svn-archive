
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_read_mod
   implicit none
   contains

   subroutine read_txt_ifile( Lfilename )
      use debug_common
      implicit none
      character(len=*), intent(in) :: Lfilename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      integer :: i
      integer, save :: first_time_caller=1
         open(unit=1,file=Lfilename//'.dat', status="old",action="read", iostat=gopenstatus )
            if(gopenstatus==gok) then
               read(1,*)
               read (1,*) Nvars; read (1,*)
               do i=1,Nvars
                  read (1,*)
               enddo
               read (1,*); read (1,*), dimx
               read (1,*); read (1,*), dimy 
            else
               write (*,*), Lfilename//'.dat',' NOT found to read ini'
            endif
         close(1)
         dimx_i(first_time_caller) = dimx
         dimx_tot = dimx_tot + dimx
!         print *, 'dim_tot, nodes  ',dimx_tot, n_nodes, first_time_caller
         if (first_time_caller==n_nodes) then
            allocate( ar_Nvars( Nvars*dimx_tot ) )
            allocate( ar_data(Nvars,dimy,dimx_tot) )
         endif
         first_time_caller = first_time_caller + 1
      return 
   end subroutine read_txt_ifile 

   !==========================================================================!
   !==========================================================================!


   !==========================================================================!
   !==========================================================================!

   subroutine read_txt_file( Lfilename, ndimx )
      use debug_common
      implicit none
      integer, intent(out) :: ndimx 
      character(len=*), intent(in) :: Lfilename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      character(len=99) :: trash 
      integer :: i
      integer, save :: first_time_caller=1
         open(unit=1,file=Lfilename//'.dat', status="old",action="read", iostat=gopenstatus )
            if(gopenstatus==gok) then
               read(1,*)
               read (1,*) Nvars; read (1,*), trash
               if (first_time_caller==1)  allocate( ar_varname(Nvars) )
               do i=1,Nvars
                  read (1,*), ar_varname(i)
               enddo
               read (1,*), trash
               read (1,*), dimx 
               read (1,*), trash
               read (1,*), dimy 
            else
               write (*,*), Lfilename//'.dat',' NOT found to read'
            endif
         close(1)
         first_time_caller = first_time_caller + 1
         ndimx=dimx
      return 
   end subroutine read_txt_file 

   !==========================================================================!
   !==========================================================================!

   subroutine read_dat_file( Lfilename, Lfrom,Lto )
      use debug_common
      implicit none
      integer, intent(in) :: Lfrom, Lto
      character(len=*), intent(in) :: Lfilename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      integer :: i,j
!      integer(i_d) :: frecl
!      frecl = Nvars * dimx*r_1
         open(unit=2,file=Lfilename//'.bin',status="unknown",action="read", &
                  iostat=gopenstatus, form="unformatted" )
            if(gopenstatus==gok) then
               do i=1,dimy
                  !read(2), ar_Nvars(1:Lto-Lfrom )
                  read(2), ar_Nvars(1:Nvars*dimx )
                  do j=1,Nvars
                     ar_data(j,i,Lfrom:Lto) = ar_Nvars( ( (j-1)*dimx )+1 : j*dimx )
!                     print *, ar_data(j,i,Lfrom:Lto)
!                     print *, 'Lfrom:Lto ',Lfrom,Lto
                  enddo 
               enddo
            else
               write (*,*), Lfilename//'.bin',' NOT found for read'
            endif
         close(2)
      return 
   end subroutine read_dat_file 
   !==========================================================================!
   !==========================================================================!


end module debug_read_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
