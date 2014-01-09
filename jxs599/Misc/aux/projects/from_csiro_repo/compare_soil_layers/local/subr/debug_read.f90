
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_read_mod
   use debug_common
   implicit none
   contains
      
   !==========================================================================!
   !==========================================================================!

   subroutine read_args
      use debug_common
      implicit none
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
         open(unit=1,file='input.dat', status="unknown",action="read", iostat=gopenstatus )
            if(gopenstatus==gok) then
               read(1,*), filename1
               read(1,*), filename2
               read(1,*), filename3
               read(1,*), filename4
               read(1,*), write_flag
               read(1,*), plot_flag
               read(1,*), t_window
               read(1,*), x_window
            else
               write (*,*), 'input.dat',' NOT found to read'
            endif
         close(1)
      return
   end subroutine read_args
      
 
   !==========================================================================!
   !==========================================================================!

   subroutine read_txt_file( filename )
      implicit none
      character(len=*), intent(in) :: filename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      character(len=99) :: trash 
      integer, save :: statalloc_read=1, statalloc_cmp=1, statalloc_name=1, statalloc_nvars=1
      integer :: i
         open(unit=1,file=filename//'.dat', status="old",action="read", iostat=gopenstatus )
            if(gopenstatus==gok) then
               read(1,*)
               read (1,*) Nvars; read (1,*), trash
               if(statalloc_name >0) then                                         
                 allocate( ar_varname(Nvars), stat=statalloc_name )
               endif
               do i=1,Nvars
                  read (1,*), ar_varname(i)
               enddo
               read (1,*), trash; read (1,*), dimx 
               read (1,*), trash; read (1,*), dimy 
            else
               write (*,*), filename//'.dat',' NOT found to read'
            endif
         close(1)
      return 
   end subroutine read_txt_file 

   !==========================================================================!
   !==========================================================================!

   subroutine read_dat_file( filename )
      implicit none
      character(len=*), intent(in) :: filename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      integer :: i,j
!      integer(i_d) :: frecl
!      frecl = Nvars * dimx*r_1
         open(unit=2,file=filename//'.bin',status="unknown",action="read", &
                  iostat=gopenstatus, form="unformatted" )
            if(gopenstatus==gok) then
               do i=1,dimy
                  read(2), ar_Nvars(1:(Nvars*dimx))
                  do j=1,Nvars
                     ar_data(j,i,:) = ar_Nvars( ( (j-1)*dimx )+1 : j*dimx )
!if(j==1)  print *, ar_data(j,i,:) 
                  enddo 
               enddo
            else
               write (*,*), filename//'.bin',' NOT found for read'
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
