

module debug_read_mod
   implicit none
   
   interface read_dat_file
      module procedure read_dat_file1, read_dat_file2
   end interface 
   
   contains

   !==========================================================================!
   !==========================================================================!

   subroutine read_file( filename, varname, dimx, dimy  )
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(out) :: varname
      integer,intent(out) :: dimx, dimy 
      
      integer, parameter :: gok=0
      integer :: gopenstatus
      character(len=99) :: trash
     
      open(unit=1,file=filename//'.dat', status="old",action="read", iostat=gopenstatus )
         if(gopenstatus==gok) then
            read(1,*) trash 
            read (1,*) trash 
            read (1,*), trash
            read (1,*), varname
            read (1,*), trash
            read (1,*), dimx 
            read (1,*), trash
            read (1,*), dimy 
         else
            write (*,*), filename//'.dat',' NOT found to read'
         endif
      close(1)

      return 
   end subroutine read_file 

   !==========================================================================!
   !==========================================================================!


   subroutine read_dat_file1( filename,ar_data, dimx, dimy )
      implicit none
      character(len=*), intent(in) :: filename
      !interfacing var
      real*8, intent(out), dimension(:) :: ar_data
      integer,intent(in) :: dimx, dimy 
      
      integer, parameter :: gok=0
      integer :: gopenstatus
      integer :: i,j
     
      open(unit=2,file=filename//'.bin',status="unknown",action="read", &
               iostat=gopenstatus, form="unformatted" )
         if(gopenstatus==gok) then
            do i=1,dimy
               read(2), ar_data
            enddo
         else
            write (*,*), filename//'.bin',' NOT found for read'
         endif
      close(2)

      return 
   end subroutine read_dat_file1 

   !==========================================================================!
   !==========================================================================!



   subroutine read_dat_file2( filename, ar_data, dimx, dimy )
      implicit none
      character(len=*), intent(in) :: filename
      integer,intent(in) :: dimx, dimy 
      !interfacing var
      real*8, intent(out), dimension(:,:) :: ar_data
      real*8, dimension(dimx) :: temp
      
      integer, parameter :: gok=0
      integer :: gopenstatus
      integer :: i,j
     
      open(unit=2,file=filename//'.bin',status="unknown",action="read", &
               iostat=gopenstatus, form="unformatted" )
         if(gopenstatus==gok) then
            do i=1,dimy
print *, 'jhan',i,dimy,temp
               read(2), temp
               ar_data(i,:) = temp 
            enddo
         else
            write (*,*), filename//'.bin',' NOT found for read'
         endif
      close(2)

      return 
   end subroutine read_dat_file2 

   !==========================================================================!
   !==========================================================================!















!   subroutine read_dat_file2( filename,ar_data )
!      implicit none
!      real, intent(out), dimension(:,:) :: ar_data
!      character(len=*), intent(in) :: filename
!      integer, parameter :: gok=0
!      integer :: gopenstatus
!      integer :: i,j
!!      integer :: frecl
!!      frecl = Nvars * dimx*r_1
!         open(unit=2,file=filename//'00.bin',status="unknown",action="read", &
!                  iostat=gopenstatus, form="unformatted" )
!            if(gopenstatus==gok) then
!               do i=1,dimy
!                  read(2), ar_Nvars
!                     ar_data(i,:) = ar_Nvars( 1 : dimx )
!               enddo
!            else
!               write (*,*), filename//00'.bin',' NOT found for read'
!            endif
!         close(2)
!      return 
!   end subroutine read_dat_file2 

   !==========================================================================!
   !==========================================================================!

end module debug_read_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
