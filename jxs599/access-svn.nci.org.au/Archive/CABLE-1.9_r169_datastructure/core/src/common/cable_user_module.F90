!================================================================================
!=== This module is intended to enable users to easily test functionality without 
!=== having to follow the technical layout of the CABLE code. This module 
!=== therefore provides a means for the user to:
!=== 1. specify/call/use alternative subroutines
!=== 2. pass variables
!=== ----------------------------------------------------------------
!=== 1. in some cases it might be clearer to include the new (version of a)
!===    subroutine and "USE" this module where appropriate
!=== 2. the user might want access to a variable which exists at some higher level 
!=== of the code, but is not available in that particular function. for the sake 
!=== of clarity and ease, rather that passing this variable through argument lists 
!=== this module makes it possible to access any variable, anywhere. NB. this 
!=== is NOT automatically packed into CABLE format as is done through the 
!=== interface. To use this feature:
!===     * declare the variable name in the header of this module (as below) as 
!===     you wish to refer to it. (=var)
!===     * "use cable_user_module" in the subroutine from which the variable 
!===     originates.
!===     * "call cable_user_var(dim1 [, dim2, dim3], var)" in that subroutine.
!===     dim1 (and optionally dim2, dim3) refer to the shape of the variable
!===     at present only "REAL" variables are catered for of up to 3 dimensions 
!===     * then in that subroutine, where the existing variable you wish passed 
!===     is called for e.g. "x", include the line "var = x" 
!===     * finally "use cable_user_module" in the subroutine where you want to 
!===     access the variable (var)
!=== **NB** 
!=== this makes a copy of "x" and changing "var" does not feed back into "x" 
!================================================================================

module cable_user_module
   implicit none
   public  
  
   !--- variable declarations which can be seen wherever this module is "USED"
   !--- for example 
   !real, allocatable, dimension(:,:) :: tstar_tile_eg
   
   !--- in this example we wish to access the "real" variable 
   !---     tstar_tile(land_pts,ntile)
   !--- which exists in sf_exch.F90. we will refer to it elsewher in the code as
   !---     tstar_tile_eg
   !--- thus in sf_exch() we add the lines. 
   !--- .........
   !---     use cable_user_module 
   !--- .........
   !--- .........
   !---     call cable_user_var(land_pts,ntiles,tstar_tile_eg) 
   !---     tstar_tile_eg = tstar_tile
   !--- then to use access this var in any subroutine simply
   !---     use cable_user_module 
   !--- and do whatever you like with tstar_tile_eg 

   interface cable_user_var
      module procedure alloc_user_var1,alloc_user_var2, alloc_user_var3
   end interface

    
   contains

   subroutine alloc_user_var1(dim1,var)
      use cable_common_module
      implicit none 
      integer :: dim1   
      real, allocatable, dimension(:) :: var 
         if (.NOT. allocated(var))  then
            allocate(var(dim1))
            var=0.
         endif   
      return
   end subroutine alloc_user_var1

   subroutine alloc_user_var2(dim1,dim2,var)
      use cable_common_module
      implicit none 
      integer :: dim1,dim2,i,j   
      real, allocatable, dimension(:,:) :: var 
         if (.NOT. allocated(var)) then 
            allocate(var(dim1,dim2))
            var=0.
         endif   
      return
   end subroutine alloc_user_var2

   subroutine alloc_user_var3(dim1,dim2,var)
      use cable_common_module
      implicit none 
      integer :: dim1,dim2,dim3   
      real, allocatable, dimension(:,:,:) :: var 
         if (.NOT. allocated(var)) then 
            allocate(var(dim1,dim2,dim3))
            var=0.
         endif   
      return
   end subroutine alloc_user_var3

end module cable_user_module
