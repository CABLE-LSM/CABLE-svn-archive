
program test_driver
   implicit none
   external subr
   !testing
   !integer, parameter ::   n=8,nsegments=8,row=73,row_length=96, levels=38
   integer, parameter ::   nsegments=2,row=2,row_length=2, levels=2
   integer :: i,j,k, first_point,step   
   integer :: isegment_x(nsegments), isegment_y(nsegments)
   integer :: ar1(row,row_length,levels) 

      do i=1,levels 
         do j=1, row_length
            do k=1,row 
               !ar1(i,j,k) = ( (i-1)*nsegments*nsegments ) + ( (j-1)*nsegments )  + k
               ar1(i,j,k) = ( (i-1)*row_length*row ) + ( (j-1)*row )  + k
               print *,'level, row_length, row',i,j,k,ar1(i,j,k) 
            end do              
         end do              
      end do   
       
      first_point = 1
      step = row*row_length/ nsegments
      Do i = 1, nsegments
         isegment_y(i) = (first_point-1)/row_length + 1
         isegment_x(i) = first_point - ( isegment_y(i) - 1 ) * row_length
         first_point = first_point+step
      end do    

      Do i = 1, nsegments
               print *, isegment_y(i) 
               print *, isegment_x(i) 
      end do    

       
               i=1 
      !call test_subr( ar1(isegment_x(i),isegment_y(i),1), levels, row,row_length)
      call test_subr( ar1(isegment_x(i):isegment_x(i+1),isegment_y(i):isegment_y(i+1),:), levels, row/nsegments,row_length/nsegments)

      print *,'' 
      do i=1,levels 
         do j=1, row_length
            do k=1,row 
               !print *,ar1(i,j,k) 
            end do              
         end do              
      end do   
   stop
end program test_driver

