!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine top_panel
      use genericread_module
      use plottypes_module
      use plotcases_module
      implicit none
      type (fiducial_case) :: lambda,bh5,bh6,srb005,srb010,srb040,srb080
      real ::  xmin=0.0,xmax=120,ymin=0.02,ymax=0.55   
         call setup_plot(xmin,xmax,ymin,ymax)  
         call readfile("Input/lambda.dat",29,lambda%x)   !filename,N,whereto
         call plotcase(lambda%x,srb005,N,1,3,1,"Input/myGamma.005") 
         call pgbox('BNCST', 0.0, 0, 'BSTNC', 0.0,0)  !draw box art end
      return
   end subroutine         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






















   subroutine bottom_panel
      implicit none
      real ::  xmin=0.0,xmax=120,ymin=0.02,ymax=0.55   
         call setup_plot(xmin,xmax,ymin,ymax)    
         call pgbox('BNCST', 0.0, 0, 'BSTNC', 0.0,0)  !draw box art end
      return
   end subroutine         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine setup_plot(xmin,xmax,ymin,ymax)    
      implicit none
      real ::  xmin,xmax,ymin,ymax   
         call pgpage()
         call pgmtxt("B",2.0,0.5,0.5,"\gl\\dMFP\\u (Mpc)")
         call pgmtxt('L',2.0,0.5,0.5,"\gG\\d12\\u")
         call pgswin(xmin,xmax,ymin,ymax)
      return
   end subroutine 
