
module cable_dupVars_mod
   implicit none
   private
   public cable_dupVArs
   public cable_setDupVars 
   public cable_resetDupVars 
   
   interface cable_dupVars
      module procedure dupvars4, dupvars8, dupvarsL, &
               idupvars4, idupvars8, dupvars4_2d, dupvars8_2d, &
               idupvars4_2d, idupvars8_2d, dupvars4_3d, dupvars8_3d
   end interface cable_dupVArs
    
   interface cable_setdupvars
      module procedure cable_setdupvars4, cable_setdupvars8, cable_setdupvarsL, &
               icable_setdupvars4, icable_setdupvars8, cable_setdupvars4_2d,     &
               cable_setdupvars8_2d, icable_setdupvars4_2d, icable_setdupvars8_2d, &
               cable_setdupvars4_3d, cable_setdupvars8_3d
   end interface cable_setdupvars
    
   interface cable_resetdupvars
      module procedure cable_resetdupvars4, cable_resetdupvars8, cable_resetdupvarsL, &
               icable_resetdupvars4, icable_resetdupvars8, cable_resetdupvars4_2d,     &
               cable_resetdupvars8_2d, icable_resetdupvars4_2d, icable_resetdupvars8_2d, &
               cable_resetdupvars4_3d, cable_resetdupvars8_3d
   end interface cable_resetdupvars
    
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

   !elemental subroutine dupvars4(x, n)
   subroutine dupvars4(x, n, kend, alloc_call )
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
         
   end subroutine dupvars4
   
!------------------------------------------------------------------------------ 

   subroutine dupvars8(x, n, kend, alloc_call )
      real(kind=r8), dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
         
   end subroutine dupvars8

!------------------------------------------------------------------------------ 

   subroutine idupvars4(x, n, kend, alloc_call )
      integer(kind=i4), dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
      
   end subroutine idupvars4
   
!------------------------------------------------------------------------------ 

   subroutine idupvars8(x, n, kend, alloc_call )
      integer(kind=i8), dimension(:,: ), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
         
   end subroutine idupvars8

!------------------------------------------------------------------------------ 

   subroutine dupvars4_2d(x, n, m, kend, alloc_call )
      real(kind=r4), dimension(:,:,:), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine dupvars4_2d
   
!------------------------------------------------------------------------------ 

   subroutine dupvars8_2d(x, n, m,kend,  alloc_call )
      real(kind=r8), dimension(:,:,: ), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine dupvars8_2d

!------------------------------------------------------------------------------ 

   subroutine idupvars4_2d(x, n, m, kend, alloc_call )
      integer(kind=i4), dimension(:,:,: ), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine idupvars4_2d
  
!------------------------------------------------------------------------------ 
   subroutine idupvars8_2d(x, n, m, kend, alloc_call )
      integer(kind=i8), dimension(:,:,:), allocatable :: x
      integer :: n, m  
      integer :: kend  
      integer :: i, j  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, kend) )
         endif
         
   end subroutine idupvars8_2d

!------------------------------------------------------------------------------ 

   subroutine dupvars4_3d(x, n, m, p, kend, alloc_call )
      real(kind=r4), dimension(:,:,:,:), allocatable :: x
      integer :: n, m, p  
      integer :: kend  
      integer :: i, j, k  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, p, kend) )
         endif
         
   end subroutine dupvars4_3d

!------------------------------------------------------------------------------ 

   subroutine dupvars8_3d(x, n, m, p, kend, alloc_call )
      real(kind=r8), dimension(:,:,:,:), allocatable :: x
      integer :: n, m, p  
      integer :: kend  
      integer :: i, j, k  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, m, p, kend) )
         endif
         
   end subroutine dupvars8_3d

!------------------------------------------------------------------------------ 

   subroutine dupvarsL(x, n, kend, alloc_call )
      logical, dimension(:,:), allocatable :: x
      integer :: n  
      integer :: kend  
      integer :: i  
      logical, optional :: alloc_call
      
         if( present( alloc_call ) ) then 
            allocate( x(n, kend) )
         endif
         
   end subroutine dupvarsL
 
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 

   subroutine Cable_setDupVars4( dup, var,  ktau, alloc_call ) 
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
          
   end subroutine cable_setDupVars4
   
!------------------------------------------------------------------------------ 

   subroutine Cable_setDupVars8( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:) :: dup
      real(kind=r8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:, ktau) = var(:) 
          
   end subroutine cable_setDupVars8
 
!------------------------------------------------------------------------------ 
   subroutine iCable_setDupVars4( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:) :: dup
      integer(kind=i4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         dup(:, ktau) = var(:) 
          
   end subroutine icable_setDupVars4
   
!------------------------------------------------------------------------------ 

   subroutine iCable_setDupVars8( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:) :: dup
      integer(kind=i8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:, ktau) = var(:) 
          
   end subroutine icable_setDupVars8

!------------------------------------------------------------------------------ 

   subroutine Cable_setDupVars4_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:) :: dup
      real(kind=r4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine cable_setDupVars4_2d
   
!------------------------------------------------------------------------------ 

   subroutine Cable_setDupVars8_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:) :: dup
      real(kind=r8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine cable_setDupVars8_2d

!------------------------------------------------------------------------------ 

   subroutine iCable_setDupVars4_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:,:) :: dup
      integer(kind=i4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine icable_setDupVars4_2d
   
!------------------------------------------------------------------------------ 
   
   subroutine iCable_setDupVars8_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:,:) :: dup
      integer(kind=i8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:, ktau) = var(:,:) 
          
   end subroutine icable_setDupVars8_2d

!------------------------------------------------------------------------------ 

   subroutine Cable_setDupVars4_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:,:) :: dup
      real(kind=r4), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:,:, ktau) = var(:,:,:) 
          
   end subroutine cable_setDupVars4_3d

!------------------------------------------------------------------------------ 

   subroutine Cable_setDupVars8_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:,:) :: dup
      real(kind=r8), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         dup(:,:,:, ktau) = var(:,:,:) 
          
   end subroutine cable_setDupVars8_3d

!------------------------------------------------------------------------------ 

   subroutine cable_setdupvarsL( dup, var,  ktau, alloc_call ) 
      logical, dimension(:,:) :: dup
      logical, dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      
         dup(:, ktau) = var(:) 
         
   end subroutine cable_setdupvarsL
 
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 

   subroutine cable_resetDupVars4( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:) :: dup
      real(kind=r4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      !integer :: n  
      !integer :: i  
         
         !if( present( alloc_call ) ) then 
         !   allocate( x(n, kend) )
         !endif
         
         var(:) = dup(:, ktau)
          
   end subroutine cable_resetDupVars4
   
!------------------------------------------------------------------------------ 

   subroutine cable_resetDupVars8( dup, var, ktau, alloc_call ) 
      real(kind=r8), dimension(:,:) :: dup
      real(kind=r8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:) = dup(:, ktau)
          
   end subroutine cable_resetDupVars8
 
!------------------------------------------------------------------------------ 
   subroutine iCable_resetDupVars4( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:) :: dup
      integer(kind=i4), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
          
         var(:) = dup(:, ktau)

   end subroutine iCable_resetDupVars4
   
!------------------------------------------------------------------------------ 

   subroutine iCable_resetDupVars8( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:) :: dup
      integer(kind=i8), dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:) = dup(:, ktau)
          
   end subroutine iCable_resetDupVars8

!------------------------------------------------------------------------------ 

   subroutine cable_resetDupVars4_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:) :: dup
      real(kind=r4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine cable_resetDupVars4_2d
   
!------------------------------------------------------------------------------ 

   subroutine cable_resetDupVars8_2d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:) :: dup
      real(kind=r8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine cable_resetDupVars8_2d

!------------------------------------------------------------------------------ 

   subroutine iCable_resetDupVars4_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i4), dimension(:,:,:) :: dup
      integer(kind=i4), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine iCable_resetDupVars4_2d
   
!------------------------------------------------------------------------------ 
   
   subroutine iCable_resetDupVars8_2d( dup, var,  ktau, alloc_call ) 
      integer(kind=i8), dimension(:,:,:) :: dup
      integer(kind=i8), dimension(:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:) = dup(:,:, ktau)
          
   end subroutine iCable_resetDupVars8_2d

!------------------------------------------------------------------------------ 

   subroutine cable_resetDupVars4_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r4), dimension(:,:,:,:) :: dup
      real(kind=r4), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:,:) = dup(:,:,:, ktau)
          
   end subroutine cable_resetDupVars4_3d

!------------------------------------------------------------------------------ 

   subroutine cable_resetDupVars8_3d( dup, var,  ktau, alloc_call ) 
      real(kind=r8), dimension(:,:,:,:) :: dup
      real(kind=r8), dimension(:,:,:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
         
         var(:,:,:) = dup(:,:,:, ktau)
          
   end subroutine cable_resetDupVars8_3d

!------------------------------------------------------------------------------ 

   subroutine cable_resetdupvarsL( dup, var,  ktau, alloc_call ) 
      logical, dimension(:,:) :: dup
      logical, dimension(:) :: var 
      integer :: ktau  
      logical, optional :: alloc_call
      
         var(:) = dup(:, ktau)
         
   end subroutine cable_resetdupvarsL
 
!------------------------------------------------------------------------------ 





end module cable_dupvars_mod
