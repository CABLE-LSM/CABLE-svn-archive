

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
      integer, parameter :: gok=0
      integer :: gopenstatus
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
!         print *, 'dim_tot, nodes  ',dimx_tot, i_nodes, first_time_caller
         if (first_time_caller==i_nodes) then
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
      integer, parameter :: gok=0
      integer :: gopenstatus
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
      integer, parameter :: gok=0
      integer :: gopenstatus
      integer :: i,j

         open(unit=2,file=Lfilename//'.bin',status="unknown",action="read", &
                  iostat=gopenstatus, form="unformatted" )
            if(gopenstatus==gok) then

               do i=1,dimy
   
                  read(2), ar_Nvars(1:dimx )
                  ar_data(1,i,Lfrom:Lto) = ar_Nvars( 1 : dimx )

               enddo

            else
               write (*,*), Lfilename//'.bin',' NOT found for read'
            endif

         close(2)

   end subroutine read_dat_file 
   !==========================================================================!
   !==========================================================================!


end module debug_read_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
