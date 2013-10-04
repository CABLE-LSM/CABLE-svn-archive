!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module plottypes_module
      type fiducial_case 
         character (len=99) :: filename 
         integer :: N
         real, pointer :: x(:)
         integer :: sci,slw,sls
      end type
      type fesc_case 
         character (len=99) :: filename 
         integer :: N
         real, pointer :: x(:)
         real :: fesc
         integer :: sci,slw,sls
      end type
      type mypoints
         character (len=99) :: filename 
         integer :: N,errdirn=1
         logical :: uniform=.false.
         real, pointer :: x(:),upper(:),lower(:)
         integer :: sci,slw,sls,marker
         real :: y,z
      end type
      type BHpoints
         character (len=99) :: filename 
         integer :: N,errdirn=1,Nlambda
         logical :: uniform=.false.
         real, pointer :: x(:),upper(:),lower(:)
         integer :: sci,slw,sls,marker
         real :: y,z
         type (fiducial_case) :: analytica,analyticb
      end type
      type Fanpoints
         character (len=99) :: filename 
         integer :: N,errdirn=1
         logical :: uniform=.false.
         real, pointer :: x(:),upper(:),lower(:)
         integer :: sci,slw,sls,marker
         real :: y,z
      end type
   end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module plotcases_module
      interface plotcase
         module procedure plot_fiducial,plotfesc,plotBHpoints,plotmypoints,plotFanpoints 
      end interface 
      contains 
         subroutine plot_casex1(lambda,ptype,N,sci,slw,sls,filename)
            use genericread_module
            use plot_types_module
            implicit none
            type (fiducial_case) :: ptype
            integer :: N,sci,sls,slw
            character (len=*) :: filename
            real, dimension(N): lambda
!               ptype%filename=filename
!               ptype%N=N;  ptype%sci=sci;ptype%slw=slw; ptype%sls=sls
!               call plotcase(lambda%x,ptype)
!

               call readfile(casex%filename,casex%N,casex%x)
               call pgsci(casex%sci); call pgsls(casex%sls); call pgslw(casex%slw) 
               call pgline(casex%N,lambda%x,casex%x) 
            return
         end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
          subroutine plotfesc(lambda,casex,ifesc) 
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
!      interface plotcase
!         module procedure plot_casex1,plot_casex2 
!      end interface
!      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         subroutine plot_casex1(lambda,ptype,N,sci,slw,sls,filename)
!            use plot_types_module
!            use plot_interface
!            implicit none
!            integer :: N,sci,sls,slw
!            character (len=*) :: filename
!            type (fiducial_case) :: lambda,ptype
!               ptype%filename=filename
!               ptype%N=N;  ptype%sci=sci;ptype%slw=slw; ptype%sls=sls
!               call plotcase(lambda%x,ptype)
!            return
!         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         subroutine plot_casex2(ptype,N,sci,slw,sls,filename,uniform,marker)
!            use plot_types_module
!            use plot_interface
!            implicit none
!            logical :: uniform
!            integer :: N,sci,sls,slw,marker
!            character (len=*) :: filename
!            type (mypoints) :: ptype
!               ptype%filename=filename;ptype%uniform=uniform;ptype%marker=marker
!               ptype%N=N;  ptype%sci=sci;ptype%slw=slw; ptype%sls=sls
!               call plotcase(ptype)
!            return
!         end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   end module    
!



