
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_read_mod
   implicit none
   contains

   !==========================================================================!
   !==========================================================================!

   subroutine read_txt_file( Lfilename )
      use debug_common
      implicit none
      character(len=*), intent(in) :: Lfilename
      integer(i_d), parameter :: gok=0
      integer(i_d) :: gopenstatus
      character(len=99) :: trash 
      integer :: i
         open(unit=1,file=Lfilename//'.dat', status="old",action="read", iostat=gopenstatus )
            if(gopenstatus==gok) then
               read(1,*)
               read (1,*) Nvars; read (1,*), trash
               allocate( ar_varname(Nvars) )
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
         allocate( ar_Nvars( Nvars*dimx ) )
         allocate( ar_data(Nvars,dimy,dimx) )
      return 
   end subroutine read_txt_file 

   !==========================================================================!
   !==========================================================================!

   subroutine read_dat_file( Lfilename )
      use debug_common
      implicit none
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
                  read(2), ar_Nvars
                  do j=1,Nvars
                     ar_data(j,i,:) = ar_Nvars( ( (j-1)*dimx )+1 : j*dimx )
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
