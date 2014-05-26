
module cable_BufVars_mod
   implicit none
   private
   public cable_BufVars
   public cable_setBufVars 
   public cable_testBufVars 
   
   interface cable_BufVars
      module procedure BufVars4, BufVars8, BufVarsL, &
               iBufVars4, iBufVars8, BufVars4_2d, BufVars8_2d, &
               iBufVars4_2d, iBufVars8_2d, BufVars4_3d, BufVars8_3d
   end interface cable_BufVars
    
   interface cable_setBufVars
      module procedure cable_setBufVars4, cable_setBufVars8, cable_setBufVarsL, &
               icable_setBufVars4, icable_setBufVars8, cable_setBufVars4_2d,     &
               cable_setBufVars8_2d, icable_setBufVars4_2d, icable_setBufVars8_2d, &
               cable_setBufVars4_3d, cable_setBufVars8_3d
   end interface cable_setBufVars
    
   interface cable_testBufVars
      module procedure cable_testBufVars4, cable_testBufVars8, cable_testBufVarsL, &
               icable_testBufVars4, icable_testBufVars8, cable_testBufVars4_2d,     &
               cable_testBufVars8_2d, icable_testBufVars4_2d, icable_testBufVars8_2d, &
               cable_testBufVars4_3d, cable_testBufVars8_3d
   end interface cable_testBufVars
    
   integer, parameter :: r4 = selected_real_kind(6)
   integer, parameter :: r8 = selected_real_kind(10)
   integer, parameter :: i4 = selected_int_kind(6)
   integer, parameter :: i8  = selected_int_kind(12)
   
   !!integer(kind=i4) inan4
   !!real(kind=r4) nan4
   !!equivalence (nan4,inan4)

   !!! Signalling NaN value from Kahan's IEEE 754 description.
   !!data inan4 / Z'7F800001' /

   !!integer(kind=i4) :: inan8(2)
   !!real(kind=r8) :: nan8
   !!equivalence (nan8,inan8(1))
  
   !!! For little endian machine, use this reversed order.
   !!data inan8 / Z'00000001', Z'7FF00000' /
   
contains

!------------------------------------------------------------------------------ 

   !elemental subroutine BufVars4(x, n)
   subroutine BufVars4(x, n, kend, alloc_call )
      real(kind=r4), dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
         
         if( .NOT. allocated(x) ) then  
            if( present( alloc_call ) ) then 
               allocate( x(n, kend) )
            endif
         endif
         
   end subroutine BufVars4
   
!------------------------------------------------------------------------------ 

   subroutine BufVars8(x, n, kend, alloc_call )
      real(kind=r8), dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
         
   end subroutine BufVars8

!------------------------------------------------------------------------------ 

   subroutine iBufVars4(x, n, kend, alloc_call )
      integer(kind=i4), dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
      
   end subroutine iBufVars4
   
!------------------------------------------------------------------------------ 

   subroutine iBufVars8(x, n, kend, alloc_call )
      integer(kind=i8), dimension(:,: ), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
         
   end subroutine iBufVars8

!------------------------------------------------------------------------------ 

   subroutine BufVars4_2d(x, n, m, kend, alloc_call )
      real(kind=r4), dimension(:,:,:), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine BufVars4_2d
   
!------------------------------------------------------------------------------ 

   subroutine BufVars8_2d(x, n, m,kend,  alloc_call )
      real(kind=r8), dimension(:,:,: ), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine BufVars8_2d

!------------------------------------------------------------------------------ 

   subroutine iBufVars4_2d(x, n, m, kend, alloc_call )
      integer(kind=i4), dimension(:,:,: ), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine iBufVars4_2d
  
!------------------------------------------------------------------------------ 
   subroutine iBufVars8_2d(x, n, m, kend, alloc_call )
      integer(kind=i8), dimension(:,:,:), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine iBufVars8_2d

!------------------------------------------------------------------------------ 

   subroutine BufVars4_3d(x, n, m, p, kend, alloc_call )
      real(kind=r4), dimension(:,:,:,:), allocatable :: x
      integer :: n, m, p  
      integer :: kend  
      integer :: i, j, k  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, p, kend) )
         endif
         
   end subroutine BufVars4_3d

!------------------------------------------------------------------------------ 

   subroutine BufVars8_3d(x, n, m, p, kend, alloc_call )
      real(kind=r8), dimension(:,:,:,:), allocatable :: x
      integer :: n, m, p  
      integer :: kend  
      integer :: i, j, k  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, p, kend) )
         endif
         
   end subroutine BufVars8_3d

!------------------------------------------------------------------------------ 

   subroutine BufVarsL(x, n, kend, alloc_call )
      logical, dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
         
   end subroutine BufVarsL
 
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 

   subroutine Cable_setBufVars4( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:) :: dup
      real(kind=r4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      !integer :: n  
      !integer :: i  
         
         !if( present( alloc_call ) ) then 
         !   allocate( x(n, kend) )
         !endif
         
         dup(:, ktau) = var(:) 
          
   end subroutine cable_setBufVars4
   
!------------------------------------------------------------------------------ 

   subroutine Cable_setBufVars8( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:) :: dup
      real(kind=r8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:, ktau) = var(:) 
          
   end subroutine cable_setBufVars8
 
!------------------------------------------------------------------------------ 
   subroutine iCable_setBufVars4( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:) :: dup
      integer(kind=i4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         dup(:, ktau) = var(:) 
          
   end subroutine icable_setBufVars4
   
!------------------------------------------------------------------------------ 

   subroutine iCable_setBufVars8( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:) :: dup
      integer(kind=i8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:, ktau) = var(:) 
          
   end subroutine icable_setBufVars8

!------------------------------------------------------------------------------ 

   subroutine Cable_setBufVars4_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:) :: dup
      real(kind=r4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine cable_setBufVars4_2d
   
!------------------------------------------------------------------------------ 

   subroutine Cable_setBufVars8_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:) :: dup
      real(kind=r8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine cable_setBufVars8_2d

!------------------------------------------------------------------------------ 

   subroutine iCable_setBufVars4_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:,:) :: dup
      integer(kind=i4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine icable_setBufVars4_2d
   
!------------------------------------------------------------------------------ 
   
   subroutine iCable_setBufVars8_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:,:) :: dup
      integer(kind=i8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine icable_setBufVars8_2d

!------------------------------------------------------------------------------ 

   subroutine Cable_setBufVars4_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:,:) :: dup
      real(kind=r4), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:,:, ktau) = var(:,:,:) 
          
   end subroutine cable_setBufVars4_3d

!------------------------------------------------------------------------------ 

   subroutine Cable_setBufVars8_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:,:) :: dup
      real(kind=r8), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:,:, ktau) = var(:,:,:) 
          
   end subroutine cable_setBufVars8_3d

!------------------------------------------------------------------------------ 

   subroutine cable_setBufVarsL( dup, var,  ktau, alloc_call ) 
      logical, dimension(:,:) :: dup
      logical, dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      
         dup(:, ktau) = var(:) 
         
   end subroutine cable_setBufVarsL
 
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 

   subroutine Cable_testBufVars4( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:) :: dup
      real(kind=r4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      !integer :: n  
      !integer :: i  
         
         var(:) = dup(:, ktau)
          
   end subroutine cable_testBufVars4
   
!------------------------------------------------------------------------------ 

   subroutine Cable_testBufVars8( dup, var, ktau, alloc_call ) 
      real(kind=r8), dimension(:,:) :: dup
      real(kind=r8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:) = dup(:, ktau)
          
   end subroutine cable_testBufVars8
 
!------------------------------------------------------------------------------ 
   subroutine iCable_testBufVars4( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:) :: dup
      integer(kind=i4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
          
         var(:) = dup(:, ktau)

   end subroutine icable_testBufVars4
   
!------------------------------------------------------------------------------ 

   subroutine iCable_testBufVars8( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:) :: dup
      integer(kind=i8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:) = dup(:, ktau)
          
   end subroutine icable_testBufVars8

!------------------------------------------------------------------------------ 

   subroutine Cable_testBufVars4_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:) :: dup
      real(kind=r4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine cable_testBufVars4_2d
   
!------------------------------------------------------------------------------ 

   subroutine Cable_testBufVars8_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:) :: dup
      real(kind=r8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine cable_testBufVars8_2d

!------------------------------------------------------------------------------ 

   subroutine iCable_testBufVars4_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:,:) :: dup
      integer(kind=i4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine icable_testBufVars4_2d
   
!------------------------------------------------------------------------------ 
   
   subroutine iCable_testBufVars8_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:,:) :: dup
      integer(kind=i8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine icable_testBufVars8_2d

!------------------------------------------------------------------------------ 

   subroutine Cable_testBufVars4_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:,:) :: dup
      real(kind=r4), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:,:) = dup(:,:,:, ktau)
          
   end subroutine cable_testBufVars4_3d

!------------------------------------------------------------------------------ 

   subroutine Cable_testBufVars8_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:,:) :: dup
      real(kind=r8), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:,:) = dup(:,:,:, ktau)
          
   end subroutine cable_testBufVars8_3d

!------------------------------------------------------------------------------ 

   subroutine cable_testBufVarsL( dup, var,  ktau, alloc_call ) 
      logical, dimension(:,:) :: dup
      logical, dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      
         var(:) = dup(:, ktau)
         
   end subroutine cable_testBufVarsL
 
!------------------------------------------------------------------------------ 





end module cable_BufVars_mod
