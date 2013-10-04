
#include "debug_directives.f90" 

module debug_comp_mod                                                          
   use debug_common                                                            
   implicit none                                                               
   real, dimension(:), allocatable :: ar_x, ar_t, ar_y, ary_twindow

   contains                                                                    

   !==========================================================================!
   !==========================================================================!
                                                                               
   subroutine x_data( )                                                 
      implicit none                                                            
      integer :: i
         allocate( ar_x(dimx) )            
         allocate( cmp_data(20,dimx) )
         do i=1,dimx                                                     
            ar_x(i) = real(i)                                            
         enddo      
      return                                                                   
   end subroutine x_data
                                                                               
   !==========================================================================!
   !==========================================================================!

   subroutine al_data( )                                                 
      implicit none                                                            
         allocate( ar_Nvars( Nvars*dimx ) )
         allocate( ar_y(dimx) )                     
         allocate( ary_twindow(dimx) )                     
         allocate( ar_data(Nvars,dimy,dimx) )
      return                                                                   
   end subroutine al_data
                                                                                
   !==========================================================================!
   !==========================================================================!
                                                                               
   subroutine comp_data( i )                                                 
      implicit none      
      integer, intent(in) :: i                                                      
      integer :: j                                                             
         do j=1,Nvars
            call debug_dat( j, i+j ) 
         enddo
      return                                                                   
   end subroutine comp_data
                                                                              
   !==========================================================================!
   !==========================================================================!
                                                                               
   subroutine de_data( )                                                 
      implicit none                                                            
         deallocate( ar_Nvars )
         deallocate( ar_data )
         deallocate( ar_y)                     
         deallocate( ary_twindow)                     
      return                                                                   
   end subroutine de_data
                                                                               
   !==========================================================================!
   !==========================================================================!
     
      subroutine debug_dat( var_j, t )
         implicit none
         integer,intent(in) :: var_j, t
         integer ::i,j,k
            ary_twindow = real(0) 
            do i=1, dimy
               ar_y = real( ar_data(var_j,i,:) )    
               do k=1,dimx                                                     
                  if(isnan( ar_y(k) ) ) then
                     ar_y(k) =  -1.0
                  endif                                    
               enddo
               ary_twindow= ary_twindow + ar_y
            enddo
            ar_y = ary_twindow / real(dimy) 
            if(x_window>1) then
               call smooth_av( dimx, ar_x, ar_y, x_window ) 
            endif
            cmp_data(t,:) = ar_y
!if (t==1) print *, cmp_data(t,1:100)
         return   
      end subroutine debug_dat
      
      !=======================================================================!
      !=======================================================================!

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

      end subroutine smooth_av
   
      !==========================================================================!
      !==========================================================================!

      real function smooth(x_0,ax,af,n,var)
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
        end do 
      end if
      smooth=f
      weight=0.0
      x=0
      delta =0.0
      f=0.0
      deallocate (aweight)
      return
      end function smooth

   !==========================================================================!
   !==========================================================================!

end module debug_comp_mod


