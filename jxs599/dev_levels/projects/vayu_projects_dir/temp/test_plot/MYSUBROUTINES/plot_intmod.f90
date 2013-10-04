!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module plot_interface
      interface plotl
         module procedure plottop,plotbottom,plotfesc_repeat,plotBHpoints,plotmypoints,plotFanpoints 
      end interface 
      contains 
         subroutine plottop(lambda,casex) 
            use genericread_module
            use plot_types_module
            implicit none
            type (fiducial_case) :: casex
            real, dimension (casex%N) :: lambda
            real ::  xmin=0.0,xmax=120,ymin=0.02,ymax=0.55   
!            real ::  xmin=10.0,xmax=1100.0,ymin=0.05,ymax=2.2   
               call readfile(casex%filename,casex%N,casex%x)
               !call pgsvp(0.1,0.6,0.1,0.9)
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               call pgmtxt("B",2.0,0.5,0.5,"\gl\\dMFP\\u (Mpc)")
               call pgmtxt('L',2.0,0.5,0.5,"\gG\\d12\\u")
!               call pgswin(0.2,log10(xmax),log10(ymin),log10(ymax))
               call pgswin(xmin,xmax,ymin,ymax)
!               call pgline(casex%N,log10(lambda),log10(casex%x)) 
               call pgline(casex%N,lambda,casex%x) 
            return
         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine plotbottom(lambda,casex) 
            use genericread_module
            use plot_types_module
            implicit none
            type (fesc_case) :: casex
            real, dimension (casex%N) :: lambda
            real ::  xmin=0.0,xmax=80,ymin=0.02,ymax=0.48   
               call readfile(casex%filename,casex%N,casex%x)
               !call pgsvp(0.1,0.6,0.1,0.9)
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               call pgmtxt("B",2.0,0.5,0.5,"\gl\\dMFP\\u (Mpc)")
               call pgmtxt('L',2.0,0.5,0.5,"\gG\\d12\\u")
               call pgswin(xmin,xmax,ymin,ymax)
               call pgline(casex%N,lambda,casex%fesc*casex%x) 
            return
         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine plotfesc_repeat(lambda,casex,ifesc) 
            use genericread_module
            use plot_types_module
            implicit none
            type (fesc_case) :: casex
            real, dimension (casex%N) :: lambda
            real :: ifesc
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               call pgline(casex%N,lambda,ifesc*casex%x) 
               call pgsci(1); call pgsls(1); call pgslw(3) 
            return
         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine plotmypoints(casex) 
            use genericread_module
            use plot_types_module
            implicit none
            type (mypoints) :: casex
            real, dimension (casex%N) :: mfp
            real :: correction_factor = (0.38/0.7)   !!correction for SED integral  
               mfp=(/28.7,18.9,4.2/)   !!these are in Mpc
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               if(casex%uniform) then
                  call readfile(casex%filename,casex%N,casex%x,casex%uniform)
                  call pgpt(casex%N,mfp,correction_factor*casex%x,casex%marker)
               else
                  call readfile(casex%filename,casex%N,casex%x,casex%lower,casex%upper,casex%errdirn)
                  call pgerry(casex%N,mfp,correction_factor*casex%upper,correction_factor*casex%lower,6.0)
                  call pgslw(10)              
                  call pgpt(casex%N,mfp,correction_factor*casex%x,casex%marker)
               endif
               call pgsci(1); call pgsls(1); call pgslw(3) 
               casex%uniform=.false.
            return
         end subroutine 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine plotFanpoints(lambda,casex) 
            use genericread_module
            use plot_types_module
            implicit none
            type (Fanpoints) :: casex
            real, dimension (casex%N) :: lambda
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               call readfile(casex%filename,casex%N,casex%lower,casex%x,casex%upper,casex%errdirn)
               call pgerry(casex%N,lambda,casex%upper,casex%lower,6.0)
               call pgpt(casex%N,lambda,casex%x,casex%marker)
               call pgsci(1); call pgsls(1); call pgslw(3) 
            return
         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine plotBHpoints(lambda,casex) 
            use genericread_module
            use plot_types_module
            implicit none
            type (BHpoints) :: casex
            real, dimension (casex%N) ::lambda 
            real, dimension (2) :: mfp=(/37.0,114.0/) 
            real, dimension (2) :: Gammay=(/0.19,0.52/) 
            real, dimension (2) :: upperx=(/56.0,163.0/) 
            real, dimension (2) :: lowerx=(/22.0,70.0/) 
            real, dimension (2) :: uppery=(/0.34,0.87/) 
            real, dimension (2) :: lowery=(/0.09,0.31/) 
            integer :: i,j
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               call pgerry(casex%N,mfp,upperx,lowerx,6.0)
               call pgerry(casex%N,mfp,uppery,lowery,6.0)
               call pgpt(casex%N,mfp,Gammay,casex%marker)
               call pgsci(1); call pgsls(1); call pgslw(3) 
            return
         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         subroutine plotBHpoints(lambda,casex) 
!            use genericread_module
!            use plot_types_module
!            implicit none
!            type (BHpoints) :: casex
!            integer :: findxintercept 
!            real, dimension (casex%N) :: mfp 
!            real, dimension (casex%Nlambda) :: lambda,Gammax
!            integer :: i,j
!               call readfile(casex%filename,casex%N,casex%lower,casex%x,casex%upper,casex%errdirn)
!               do i=1,casex%N
!                  if (i==1) then
!                     Gammax=casex%analytica%x
!                  else if (i==2) then 
!                     Gammax=casex%analyticb%x
!                  endif
!                  j=findxintercept(casex%x(i),casex%Nlambda,Gammax)
!                  mfp(i)=lambda(j)
!               enddo
!               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
!               call pgerry(casex%N,mfp,casex%upper,casex%lower,6.0)
!               call pgpt(casex%N,mfp,casex%x,casex%marker)
!               call pgsci(1); call pgsls(1); call pgslw(3) 
!               do i=1,casex%N
!                  print *, mfp(i),casex%x(i)
!                  print *, 
!               enddo
!            return
!         end subroutine 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer function findxintercept(x,N,y)
            implicit none
            integer :: N,i
            real :: min_diff,x
            real :: cmin=10000.0
            real, dimension (N) :: y
               do i=1,N
                  min_diff=abs(x-y(i))
write (*,10), x,y(i),min_diff
10 format(F6.3xxx,f6.3xxx,f6.3) 
                     if (min_diff<cmin) then
                        cmin=min_diff
                        findxintercept=i
                     endif
               enddo
print *,
               cmin=10000.0
            return
         end function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         subroutine findintercept(xintercept,yintercept,N,x,y)
!            use genericread_module
!            use plot_types_module
!            implicit none
!            integer :: i
!            real :: min_diff
!            real :: cmin=10000.0
!            real :: xintercept,yintercept
!            real, dimension (N):: BHx 
!               do i=1,N
!                  min_diff=abs(Gamma_F(i)-BHx(i))
!                     if (min_diff<cmin) then
!                        cmin=min_diff
!                        xintercept=x(i)
!                        yintercept=Gamma_F(i)
!                     endif
!               enddo
!               cmin=10000.0
!         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         subroutine plotpoints(casex) 
!            use genericread_module
!            use plot_types_module
!            implicit none
!            type (datapoints) :: casex
!            real, dimension (casex%N) :: mfp
!            real :: correction_factor = (0.38/0.7)   !!correction for SED integral  
!               mfp=(/28.7,18.9,4.2/)   !!these are in Mpc
!               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
!               if(casex%uniform) then
!                  call readfile(casex%filename,casex%N,casex%x,casex%uniform)
!                  call pgpt(casex%N,mfp,correction_factor*casex%x,casex%marker)
!               else
!                  call readfile(casex%filename,casex%N,casex%x,casex%lower,casex%upper,casex%errdirn)
!                  call pgerry(casex%N,mfp,correction_factor*casex%upper,correction_factor*casex%lower,3.0)
!                  call pgslw(10)              
!                  call pgpt(casex%N,mfp,correction_factor*casex%x,casex%marker)
!               endif
!               call pgsci(1); call pgsls(1); call pgslw(3) 
!               casex%uniform=.false.
!            return
!         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
