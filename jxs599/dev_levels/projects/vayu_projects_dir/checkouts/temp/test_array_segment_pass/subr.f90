 
   subroutine test_subr(ar2, levels, row, row_length) 
      implicit none
      integer :: i,j,k
      integer :: levels,row,row_length 
      integer :: ar2(row, row_length,levels)
      integer :: ar3(row*row_length,levels)

         print *,'' 
         do i=1,levels 
            do j=1, row_length 
               do k=1, row 
                  !ar2(i,j) = ar2(i,j) * 2
                  print *, 'level, row,row_length', i,j,k, ar2(i,j,k) 
               end do
            end do
         end do

         !ar3 = reshape(ar2,(/row*row_length,levels/)) 
     
         do i=1,levels 
            do j=1, row_length 
               do k=1, row 
                  ar3(( (j-1)*row ) +k,i ) = ar2(i,j,k) 
                  !print *, 'i,j,ar3', ( (j-1)*row ) +k, ar2(i,j,k)
               end do
            end do
         end do

      
         print *,'' 
         do i=1,levels 
            do j=1, row*row_length 
                  print *, 'level, row*row_length', i,j, ar3(j,i) 
            end do
         end do

         do i=1,levels 
            do j=1, row*row_length 
                  
            end do
         end do
      
      return 
   end subroutine test_subr 


