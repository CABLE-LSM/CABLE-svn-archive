

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!jhan: reuse reading .dat code more wisely

module debug_read_mod
   implicit none
   contains

   !subr called from loop over "node" files to set global dimensions
   subroutine read_txt_ifile( Lfilename )
      use debug_common
      implicit none
      character(len=*), intent(in) :: Lfilename
      integer, parameter :: gok=0
      integer :: gopenstatus
      integer :: i
      integer, save :: first_time_caller=1

         ! .dat text file contains basic info describing the dimensions of the binary file to be read
         open(unit=1,file=Lfilename//'.dat', status="old",action="read", iostat=gopenstatus )
            if(gopenstatus==gok) then
               read(1,*)
               read (1,*) Nvars; read (1,*)
               do i=1,Nvars
                  read (1,*)
               enddo
               !dimx = spatial points
               read (1,*); read (1,*), dimx
               !dimy = temporal points
               read (1,*); read (1,*), dimy 
            else
               write (*,*), Lfilename//'.dat',' NOT found to read ini'
            endif
         close(1)

         ! spatial points per node file 
         dimx_i(first_time_caller) = dimx
         ! Sum spatial points per node file. e.g dimx_tot = mp
         dimx_tot = dimx_tot + dimx

         ! we now have enough data to alloc arrays
         ! e.g assume Nvars=1, dimx_tot = mp 
         if (first_time_caller==i_nodes) then
            ! e.g. ar_Nvars = mp, where mp= sum of "dimx-es" in each
            ! node file
            allocate( ar_Nvars( Nvars*dimx_tot ) )
            ! array to hold same data as above, per timestep

            ! force only 3 months worth of data
            dimy = dimy/4
            allocate( ar_data(Nvars,dimy,dimx_tot) )
         endif

         first_time_caller = first_time_caller + 1

      return 
   end subroutine read_txt_ifile 

   !==========================================================================!
   !==========================================================================!


   !==========================================================================!
   !==========================================================================!

   !subr called from loop over "node" files to read node file dimensions
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
         !force this as 1 year is too big
         dimy=dimy/4
      
      return 
   end subroutine read_txt_file 

   !==========================================================================!
   !==========================================================================!

   ! read binary (data) node files 
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

               ! loop over binary (data) node files PER timstep
               do i=1,dimy
   
                  ! read all spatial data per node file @ this timstep
                  read(2), ar_Nvars(1:dimx )
                  ! put spatial data per timstep into global array
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
