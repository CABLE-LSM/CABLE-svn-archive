
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_plot_mod                                                          
   use debug_common                                                            
   use debug_comp_mod, only : ar_x, ar_y, ar_t, ary_twindow
   implicit none                                                               
  contains                                                                    
                                                                              
  !==========================================================================!
  !==========================================================================!
                                                                             
  subroutine plot_data(j, n_call )                                                 
     implicit none
     integer,intent(in) :: j, n_call
     character(len=30), save :: pn1, plotname, tstr, pn2
     integer :: i
     integer(gi), save :: sci, sls, slw
     real(gs), save :: xmin, xmax, ymin, ymax,xmean,ymean
         if(n_call==1) then
            xmean = (sum(ar_x/dimx))
            xmin = minval(ar_x)-.1*xmean; xmax = maxval(ar_x) +0.1*xmean
         endif
#ifdef manual_extrema
     ymin = fymin
     ymax = fymax
#else
         if(n_call==1) then
            ymean = (sum(ar_y/dimx))
            ymin = minval(ar_y)-.1*ymean; ymax = maxval(ar_y)+0.1*ymean
         endif
#endif
         sci=1; sls=1; slw=3
         pn1 = '.png'; pn2 = '/png'
         write(tstr,10), n_call
   10    format(I5.5)   !in the form ofIw.m where w is width and m is min. places 
#ifdef testxs
         call pgopen( '/xs' )
#else
         call pgopen( trim('scratch/'// trim( ar_varname(j) )//trim(tstr)//pn1 )//pn2 )
#endif
         call pgpage()
         call pgsci(sci); call pgsls(sls); call pgslw(slw) 
         call pgswin(xmin,xmax,ymin,ymax)
        !call pgsvp(0.1,0.6,0.1,0.9)
         call pgmtxt("T",2.0,0.5,0.5,tstr)
         call pgmtxt("B",2.5,0.5,0.5,"lanpoint")
         call pgmtxt('L',2.5,0.5,0.5, trim( ar_varname(j) ) )
         call pgline(dimx,ar_x,ar_y ) 
         call pgbox('BNCST', 0.0, 0, 'BSTNC', 0.0,0)  !draw box art end
         call pgclos()
     return 
  end subroutine plot_data

   !==========================================================================!
   !==========================================================================!

end module debug_plot_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

