 
   subroutine eg_pars 
      implicit none
      integer :: n
      integer :: nstates
         do i=1, N
            x(i) = real(i/10) 
         enddo
      return 
   end subroutine eg_pars 
 
   subroutine test_plot( N, x, y ) 
      use generic_pars
      implicit none
      integer(gs), intent(in) :: N
      real(gs), intent(in), dimension(N) :: x, y 
      real(gs), parameter :: xmin=-3.0, xmax=3.0, ymin = 0.0, ymax = 2.0 
      integer(gs), parameter :: sci=1, sls = 1, slw = 3 
         call pgbeg(0,"/xs",1, 1)
         call pgsci(sci); call pgsls(sls); call pgslw(slw) 
         call pgmtxt("B",2.0,0.5,0.5,"\gl\\dMFP\\u (Mpc)")
         call pgmtxt('L',2.0,0.5,0.5,"\gG\\d12\\u")
         call pgswin(xmin,xmax,ymin,ymax)
         call pgline( N, x, y )
         call pgbox( 'BNCTS', 0.0, 0, 'BCNTS', 0.0, 0 ) 
         call pgend()
      return 
   end subroutine test_plot 
!
