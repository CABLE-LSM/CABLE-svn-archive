
#include "debug_directives.f90" 

module debug_comp_mod                                                          
   use debug_common                                                            
   implicit none                                                               
   real(gs), dimension(:), allocatable :: ar_x, ar_t, ar_y, ary_twindow
   contains                                                                    
                                                                               
   !==========================================================================!
   !==========================================================================!
                                                                               
   subroutine comp_data(j, t )                                                 
      implicit none                                                            
      integer,intent(in) :: j,t 
      integer :: i                                                             
      integer, save :: statalloc_ar_x=1, statalloc_ar_t=1, statalloc_ar_y=1    
            if(statalloc_ar_y >0) then                                         
               allocate( ar_y(dimx), stat=statalloc_ar_y )                     
               allocate( ary_twindow(dimx), stat=statalloc_ar_y )                     
            endif                                                              
            if(statalloc_ar_t >0) then                                         
               allocate( ar_t(dimy), stat=statalloc_ar_t )                          
            endif                                                              
            if(statalloc_ar_x >0) then
               allocate( ar_x(dimx), stat=statalloc_ar_x )                     
               do i=1,dimx                                                     
                  ar_x(i) = real(i)                                            
               enddo      
            endif       
            !smoother in here                                                       
            ar_y = real( ar_data(j,t,:) )    
            do i=1,dimx                                                     
               if(isnan( ar_y(i) ) ) then
                   ar_y(i) =  -1.0
               endif                                    
            enddo
      return                                                                   
   end subroutine comp_data
                                                                               
   !==========================================================================!
   !==========================================================================!

      SUBROUTINE smooth_av(nNj,z,F, var)
      implicit none
      integer(i_d), intent(in) :: nNj
      integer :: j=0,nnNj
      integer,intent(in) :: var
      real, dimension(nNj),intent(in)  :: z
      real, dimension(nNj),intent(inout)  :: F
      real, dimension(nNj) :: tempy,tempx
      !real :: smooth(),zx,mean,sub,real_var
      real :: zx,mean,sub,real_var
      real_var = real(var)
      nnNj = int(nNj)
      do j=1,nnNj
        tempy(j) = F(j)
      end do
      do j=1,nnNj
        zx=z(j)
        F(j) = smooth(zx,z,tempy,nnNj,real_var)
!        print *,F(j) 
      end do
      sub=0.0
      do j=1,nnNj
        sub =sub+F(j) 
      end do
      mean = sub/nnNj 
  !    print *,'mean  ',mean 
      return

      END SUBROUTINE smooth_av
   
      !==========================================================================!
      !==========================================================================!

      real FUNCTION smooth(x_0,ax,af,n,var)
      implicit none
      
      integer,intent(in) ::  n
      real,intent(in) :: x_0, var
      real, dimension(n),intent(in)  :: ax,af
      integer :: x
      real :: delta=0.0,weight=0.0,f=0.0
      real, dimension(:),allocatable :: aweight
      allocate(aweight(n))
      if(n<= 0) then
        print *,'ff_x_ax: n_x <= 0'
        stop
      else if(var<=0.0) then 
        print*,'smooth: var <= 0.0'
        stop
      else 
        do x=1,n
          delta = sqrt(((ax(x)-x_0)*(ax(x)-x_0)/var/var))   
!      print *,delta,ax(x), x_0 
          if(delta<20.0) then 
            aweight(x) = exp(-1*delta*delta);
          else 
            aweight(x) = 0.0
          endif
          weight = weight + aweight(x);
        end do
        if(weight==0.0) then
          print *,"ff_x_ax: interp too small"
          stop
        endif
        do x=1,n
          aweight(x) = aweight(x) / weight
          f = f + (aweight(x) * af(x))
!      print *,weight,aweight(x),af(x) 
        end do 
      end if
      smooth=f
      weight=0.0
      x=0
      delta =0.0
      f=0.0
      deallocate (aweight)
      return
      END FUNCTION smooth

   !==========================================================================!
   !==========================================================================!

end module debug_comp_mod


