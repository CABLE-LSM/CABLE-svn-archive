
module cable_alloc_setnan
   implicit none
   private
   public cable_safe_allocate 
   
   interface cable_safe_allocate
      module procedure setnan4, setnan8, setnanL, &
               isetnan4, isetnan8, setnan4_2d, setnan8_2d, &
               isetnan4_2d, isetnan8_2d, setnan4_3d, setnan8_3d
   end interface cable_safe_allocate
    
   !set via namelist
   LOGICAL  :: initnan = .TRUE.
  
   integer, parameter :: r4 = selected_real_kind(6)
   integer, parameter :: r8 = selected_real_kind(10)
   integer, parameter :: i4 = selected_int_kind(6)
   integer, parameter :: i8  = selected_int_kind(12)
   
   integer(kind=i4) inan4
   real(kind=r4) nan4
   equivalence (nan4,inan4)

   ! Signalling NaN value from Kahan's IEEE 754 description.
   data inan4 / Z'7F800001' /

   integer(kind=i4) :: inan8(2)
   real(kind=r8) :: nan8
   equivalence (nan8,inan8(1))
  
   ! For little endian machine, use this reversed order.
   data inan8 / Z'00000001', Z'7FF00000' /
   
contains

   subroutine setnan4(x, n, alloc_call )
      real(kind=r4), dimension(:), allocatable :: x
      integer :: n  
      integer :: i  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               x(i) = nan4
            enddo         
         endif

   end subroutine setnan4
   
   subroutine setnan8(x, n, alloc_call )
      real(kind=r8), dimension(:), allocatable :: x
      integer :: n  
      integer :: i  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               x(i) = nan8
            enddo         
         endif

   end subroutine setnan8

   subroutine isetnan4(x, n, alloc_call )
      integer(kind=i4), dimension(:), allocatable :: x
      integer :: n  
      integer :: i  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n) )
            endif
         endif
         
      
         if(initnan) then
            do i=1,n
               x(i) = inan4
            enddo         
         endif

   end subroutine isetnan4
   
   subroutine isetnan8(x, n, alloc_call )
      integer(kind=i8), dimension(:), allocatable :: x
      integer :: n  
      integer :: i  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n) )
            endif
         endif
         
         
         if(initnan) then
            do i=1,n
               x(i) = inan8(1)
            enddo         
         endif

   end subroutine isetnan8

   subroutine setnan4_2d(x, n, m, alloc_call )
      real(kind=r4), dimension(:,:), allocatable :: x
      integer :: n, m  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n,m) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               do j=1,m
                  x(i,j) = nan4
               enddo         
            enddo         
         endif

   end subroutine setnan4_2d
   
   subroutine setnan8_2d(x, n, m, alloc_call )
      real(kind=r8), dimension(:,:), allocatable :: x
      integer :: n, m  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n,m) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               do j=1,m
                  x(i,j) = nan8
               enddo         
            enddo         
         endif

   end subroutine setnan8_2d

   subroutine isetnan4_2d(x, n, m, alloc_call )
      integer(kind=i4), dimension(:,:), allocatable :: x
      integer :: n, m  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n,m) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               do j=1,m
                  x(i,j) = inan4
               enddo         
            enddo         
         endif

   end subroutine isetnan4_2d
   
   subroutine isetnan8_2d(x, n, m, alloc_call )
      integer(kind=i8), dimension(:,:), allocatable :: x
      integer :: n, m  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n,m) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               do j=1,m
                  x(i,j) = inan8(1)
               enddo         
            enddo         
         endif
 
   end subroutine isetnan8_2d

   subroutine setnan4_3d(x, n, m, p, alloc_call )
      real(kind=r4), dimension(:,:,:), allocatable :: x
      integer :: n, m, p  
      integer :: i, j, k  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n,m,p) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               do j=1,m
                  do k=1,p
                     x(i,j,k) = nan4
                  enddo         
               enddo         
            enddo         
         endif
 
   end subroutine setnan4_3d

   subroutine setnan8_3d(x, n, m, p, alloc_call )
      real(kind=r8), dimension(:,:,:), allocatable :: x
      integer :: n, m, p  
      integer :: i, j, k  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n,m,p) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               do j=1,m
                  do k=1,p
                     x(i,j,k) = nan8
                  enddo         
               enddo 
            enddo         
         endif
            
   end subroutine setnan8_3d

   subroutine setnanL(x, n, alloc_call )
      logical, dimension(:), allocatable :: x
      integer :: n  
      integer :: i  
      logical, optional :: alloc_call
      
         if( .NOT. allocated( x) ) then 
            if( present( alloc_call ) ) then 
               allocate( x(n) )
            endif
         endif
         
         if(initnan) then
            do i=1,n
               x(i) = .FALSE. 
            enddo         
         endif

   end subroutine setnanL
 

end module cable_alloc_setnan
