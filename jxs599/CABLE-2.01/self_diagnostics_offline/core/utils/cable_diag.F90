!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ USE this module in any subr. you wish to write vars from.             +++!
!+++ can write up to three vars of the form foo(x), per call               +++!
!+++ where x is typically the number of landpoints(tiles). binary file is  +++!
!+++ then appended every timestep with the new foo(x)                      +++!
!+++                                                                       +++! 
!+++ CALL syntax:                                                          +++!  
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, var1 )     +++!
!+++                                                                       +++! 
!+++    OR                                                                 +++! 
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, vname2, &  +++!
!+++              var1, var2 )                                             +++!      
!+++                                                                       +++! 
!+++    OR                                                                 +++! 
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, vname2, &  +++!
!+++              vname3, var1, var2, var3 )                               +++!
!+++                                                                       +++! 
!+++ where:    Nvars       number of vars being written in call (1,2 or 3) +++!
!+++           filename    base of preferred filename.dat, filename.bin    +++!
!+++           dimx        length of x-dimension (# landpoints)            +++!      
!+++           dimy        length of y-dimension (# time steps )           +++!     
!+++           timestep    # of the timestep                               +++!     
!+++           vnameX      preferred/recognizable  name of var             +++!
!+++           varX        vsr to output                                   +++!   
!+++                                                                       +++! 
!+++ to plot the temp. of the first three soil layers in subr. bar()       +++!
!+++                                                                       +++! 
!+++ e.g.   in subr. bar()                                                 +++!
!+++     .                                                                 +++! 
!+++     .                                                                 +++! 
!+++      use cable_diag_mod                                               +++! 
!+++     .                                                                 +++! 
!+++     .                                                                 +++! 
!+++      call cable_diag( 3, 'bar_T', 2950, 48, ktau, 'tsoil1',           +++! 
!+++                'tsoil2','tsoil3',tsoil(:,1), tsoil(:,2), tsoil(:,3))  +++!
!+++     .                                                                 +++! 
!+++     .                                                                 +++!
!+++ following the run binaries can be interpreted from the command line   +++!
!+++ using the cable_diag.pl command (see cable_diag.pl for details)       +++!  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


module cable_diag_module
   implicit none
   integer, parameter :: gok=0
   integer :: galloctest=1
  
   !--- subrs overloaded to respond to call cable_diag 
   interface cable_diag
      module procedure cable_diag1
   end interface cable_diag
  
   contains

   !==========================================================================!
   ! cable_diag1/2/3 call subrs to write filename.dat which contains description
   ! of data and format etc., and filename.bin containing the data   
   !==========================================================================!

   subroutine cable_diag1( Nvars, basename, dimx, dimy, timestep, node, &
                           vname1, var1 )
      integer, intent(in) :: Nvars,dimx, dimy, timestep,node
      real, intent(in), dimension(:) :: var1
      integer :: i=0
      character(len=*), intent(in) :: basename, vname1
      character(len=30) :: filename, chnode
     
         write(chnode,10) node
      10 format(i2.2)   
         filename=trim(trim(basename)//trim(chnode))
         
         if (timestep == 1) & 
            call cable_diag_desc1( Nvars, trim(filename), dimx, dimy, vname1 )
         
         call cable_diag_data1( Nvars, trim(filename), dimx, timestep, dimy, &
                                var1 )
   end subroutine cable_diag1

!=============================================================================!
!=============================================================================!

   subroutine cable_diag_desc1( Nvars, filename, dimx, dimy, vname1 )
      implicit none
      integer, intent(in) :: Nvars,dimx,dimy 
      character(len=*), intent(in) :: filename, vname1
      integer, save :: gopenstatus = 1

         open(unit=713941,file=filename//'.dat', status="replace", &
              action="write", iostat=gopenstatus )
         
         if(gopenstatus==gok) then
               write (713941,*) 'Number of var(s): '
               write (713941,*) Nvars
               write (713941,*) 'Name of var(s): '
               write (713941,7139) vname1 
   7139        format(a)            
               write (713941,*) 'dimension of var(s) in x: '
               write (713941,*) dimx 
               write (713941,*) 'dimension of var(s) in y: '
               write (713941,*) dimy 
         else
            write (*,*), filename//'.dat',' Error: unable to write'
         endif
         
         close(713941)
     
   end subroutine cable_diag_desc1


   subroutine cable_diag_data1( Nvars, filename, dimx, timestep, kend, var1  )
      implicit none
      integer, intent(in) :: Nvars, dimx, timestep, kend
      real, intent(in), dimension(:) :: var1
      character(len=*), intent(in) :: filename
      integer, save :: gopenstatus = 1

         if (timestep == 1)  then 
            open(unit=713942,file=filename//'.bin',status="unknown", &
                 action="write", iostat=gopenstatus, form="unformatted", &
                 position='append' )
         endif   
    
         if(gopenstatus==gok) then
               write (713942) var1
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif

         if (timestep == kend) & 
            close(713942)

      return 
   end subroutine cable_diag_data1

   !==========================================================================!
   !--- cable generic print status
   !==========================================================================!

  subroutine cable_stat( routname)
      use cable_common_module, only : ktau_gl, knode_gl
      implicit none
      character(len=*), intent(in) :: routname
         if(knode_gl==1) & 
            write(6,*) 'CABLE@  ', routname, ktau_gl
      return 
   end subroutine cable_stat



end module cable_diag_module



