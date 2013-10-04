
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_plot_mod                                                          
   use debug_common                                                            
   use debug_comp_mod, only : ar_x, ar_y, ar_t, ary_twindow
   implicit none                                                               
   real(gs) :: xmin, xmax, ymin, ymax,xmean,ymean
  contains                                                                    
                                                                              
  !==========================================================================!
  !==========================================================================!
                                                                             
  subroutine plot_data()                                                 
     implicit none
         xmean = (sum(ar_x/dimx))
         xmin = minval(ar_x)-.1*xmean; xmax = maxval(ar_x) +0.1*xmean
          call pgbeg(0,'/cps',3,2)

          call fplot(1, 7, 'layer 1' )
          call fplot(2, 8, 'layer 2' )
          call fplot(3,9, 'layer 3' )
          call fplot(4,10, 'layer 4' )
          call fplot(5,10, 'layer 5' )
          call fplot(6, 10, 'layer 6' )

          call pgclos()

     return 
  end subroutine plot_data

   !==========================================================================!
   !==========================================================================!
   
    subroutine fplot(n1, n2, layer)
     implicit none
     integer, intent(in) :: n1, n2
     character(len=*), intent(in) :: layer 
     integer :: p
     integer(gi) :: sci, sls, slw
         p=n1
         ymean = (sum(cmp_data(p,:)/dimx))
         ymin = minval(cmp_data(p,:))-.1*ymean; ymax = maxval(cmp_data(p,:))+0.1*ymean
print *,n1
print *,cmp_data(n1,1:10)
print *,n2
print *,cmp_data(n2,1:10)

         call pgpage()
         call pgswin(xmin,xmax,ymin,ymax)
         p=n2; sci=15; sls=1; slw=10
         !p=n2; sci=15; sls=2; slw=20
         call pgsci(sci); call pgsls(sls); call pgslw(slw) 
         call pgline(dimx,ar_x,cmp_data(p,:) ) 

         p=n1
         sci=1; sls=1; slw=3
         call pgsci(sci); call pgsls(sls); call pgslw(slw) 
         call pgline(dimx,ar_x,cmp_data(p,:) ) 
         call pgmtxt('T',2.5,0.5,0.5, layer  )
         call pgmtxt('L',2.5,0.5,0.5, 'T (K)' )
         call pgmtxt("B",2.5,0.5,0.5,"lanpoint")
         call pgbox('BNCST', 0.0, 0, 'BSTNC', 0.0,0)  !draw box art end
         return
   end subroutine fplot

  !==========================================================================!
  !==========================================================================!

end module debug_plot_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!#ifdef manual_extrema
!         ymin = fymin
!         ymax = fymax
!#else
!#endif
