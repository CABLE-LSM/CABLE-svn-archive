!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module global_module
      integer :: N=29
      real, pointer :: x(:), Gamma_BH5(:),Gamma_BH6(:),Gamma_Srb(:),Gamma_F(:)
   end module
   PROGRAM plot
      use global_module 
      use genericread_module
      implicit none
         integer :: i
         real :: fesc
            real ::  xmin=0.0,xmax=78,ymin=0.02,ymax=0.48
            call readGamma("Input/lambda.dat",N,x)
            call readGamma("Input/BHGamma.5.5",N,Gamma_BH5)
            call readGamma("Input/BHGamma.6.0",N,Gamma_BH6)
!            call plottop(xmin,xmax,ymin,ymax,fesc)
!            call plotbottom(xmin,xmax,ymin,ymax,fesc)
!            call pgend
      stop
   end program plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine plottop(xmin,xmax,ymin,ymax,fesc)
      use global_module 
      use genericread_module
      implicit none
         real :: fesc,xmin,xmax,ymin,ymax
            call plot_setup
            !call pgsvp(0.1,0.6,0.1,0.9)
               call pgmtxt("B",2.0,0.5,0.5,"\gl\\dMFP\\u (Mpc)")
               call pgmtxt('L',2.0,0.5,0.5,"\gG\\d12\\u")
            call pgswin(xmin,xmax,ymin,ymax)  
            call pgsci(15); call pgslw(5); call pgsls(4)
            call pgline(N,x,Gamma_BH5)
            call pgline(N,x,Gamma_BH6)
   !         call plotJ21
            call plotGammatop
            !call pgsvp(0.6,0.9,0.1,0.9)
            !call pgbox('BC', 0.0, 0, 'BC', 0.0,0)  
            !  call pgsvp(0.1,0.9,0.1,0.9)
            call pgsci(1); call pgslw(3); call pgsls(1)
            call pgbox('BNCST', 0.0, 0, 'BSTNC', 0.0,0)  
      return
   end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine plotbottom(xmin,xmax,ymin,ymax,fesc)
      use global_module 
      use genericread_module
      implicit none
         real :: fesc,xmin,xmax,ymin,ymax
            call pgpage
            !call pgsvp(0.1,0.6,0.1,0.9)
            call pgswin(xmin,xmax,ymin,ymax)  
            call pgsci(15); call pgslw(5); call pgsls(4)
            call plotJ21
            call plotGammabottom
            !call pgsvp(0.6,0.9,0.1,0.9)
            !call pgbox('BC', 0.0, 0, 'BC', 0.0,0)  
            !  call pgsvp(0.1,0.9,0.1,0.9)
            call pgsci(1); call pgslw(3); call pgsls(1)
            call pgbox('BNCST', 0.0, 0, 'BSTNC', 0.0,0)  
      return
      end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine plot_setup
         implicit none
               call pgbeg(0,"/cps",1,2)  
               call pgpage()
     !          call pgsvp(0.1,0.9,0.1,0.9)
               call pgsch(1.5)
         return
      end SUBROUTINE plot_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine plotGammatop
         use global_module 
         use genericread_module
         implicit none
            real :: xintercept,yintercept
            call readGamma("Input/myGamma.005",N,Gamma_F)
            call pgsci(1); call pgslw(3); call pgsls(1)
            call pgline(N,x,Gamma_F)
            call readGamma("Input/myGamma.010",N,Gamma_F)
            call pgsci(1); call pgslw(3); call pgsls(1)
            call pgline(N,x,Gamma_F)
            call readGamma("Input/myGamma.040",N,Gamma_F)
            call pgsci(1); call pgslw(3); call pgsls(1)
            call pgline(N,x,Gamma_F)
            call readGamma("Input/myGamma.080",N,Gamma_F)
            call pgsci(1); call pgslw(3); call pgsls(1)
            call pgline(N,x,Gamma_F)
      end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine plotGammabottom
         use global_module 
         use genericread_module
         implicit none
            real :: errx1,errx2,erry1,erry2,datax,datay

      subroutine plot_data()
            call readGamma("Input/myGamma.005",N,Gamma_F)
            call findintercept(errx1,erry1,Gamma_BH6)
            call findintercept(errx2,erry2,Gamma_BH5)
            datax=errx1+(0.5*abs(errx2-errx1))
            datay=erry1+(0.5*abs(erry2-erry1))
            call pgsci(3); call pgslw(3); call pgsls(1)
            call pgpt1(datax,datay,14)
            call pgerrx(1,errx1,errx2,datay,4.0)
 
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine findintercept(xintercept,yintercept,BHx)
         use global_module 
         use genericread_module
         implicit none
            integer :: i
            real :: min_diff
            real :: cmin=10000.0
            real :: xintercept,yintercept
            real, dimension (N):: BHx 
               do i=1,N
                  min_diff=abs(Gamma_F(i)-BHx(i))
                     if (min_diff<cmin) then
                        cmin=min_diff
                        xintercept=x(i)
                        yintercept=Gamma_F(i)
                     endif
               enddo
            cmin=10000.0
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine plotJ21
         use global_module 
         use genericread_module
         implicit none
            real :: fesc
            fesc=0.05
            call readGamma("Input/J21_mfp.dat",N,Gamma_Srb,fesc)
            call pgsci(15); call pgslw(5); call pgsls(1)
            call pgline(N,x,Gamma_Srb)
            fesc=0.1
            call readGamma("Input/J21_mfp.dat",N,Gamma_Srb,fesc)
            call pgsci(15); 
            call pgline(N,x,Gamma_Srb)
            fesc=0.2
            call readGamma("Input/J21_mfp.dat",N,Gamma_Srb,fesc)
            call pgsci(15); 
            call pgline(N,x,Gamma_Srb)
      end subroutine plotJ21













!
