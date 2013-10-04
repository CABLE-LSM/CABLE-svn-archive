      subroutine pgcontour(af,n_x,n_y,level,x_min,x_max,y_min,y_max,& 
     & x_axis,y_axis,string1)
      integer :: n_x,n_y
      real, dimension(n_x,n_y) :: af
      real :: level, x_min,x_max,y_min,y_max
      character(len=10) :: x_axis,y_axis,string1

      integer i,x,y
      real, dimension(6) :: atransform
      real :: delta_x, delta_y
      real, dimension(n_x*n_y) :: alocal 

      do i=1,n_x*n_y
        alocal[i]=0.0
      enddo
      if(x_axis=='l') then
        x_min = log10(x_min)
        x_max = log10(x_max)
      if(y_axis=='l') then 
        y_min = log10(y_min)
        y_max = log10(y_max)
      endif
      delta_x = (x_max - x_min) / (n_x-1)
      delta_y = (y_max - y_min) / (n_y-1)
      atransform(0) = x_min - delta_x
      atransform(1) = delta_x
      atransform(2) = 0.0
      atransform(3) = y_min - delta_y
      atransform(4) = 0.0
      atransform(5) = delta_y
      
      do x=1,n_x
      do y=1,n_y
        alocal(x+(n_x*y))=af(x)(y)
      enddo
      enddo

      call pgcont(alocal, n_x, n_y, 1, n_x, 1, n_y, &level,& 
     & - 1, atransform)
      call pgsch(1.9)
!     cpgconl(alocal, n_x, n_y, 1, n_x, 1, n_y, level, atransform, string1,100,6);
      call pgconl(alocal, n_x, n_y, 1, n_x, 1, n_y, level, atransform, string1,100,1000)

      return

      end subroutine pgcontour



    

