
program debug 
   use debug_read_mod   !only called funcs are public in this mod 
   use debug_write_mod  !only called funcs are public in this modi
   use IFPORT
   implicit none

   character(len=30) :: file_var, varname

   real*8, dimension(:,:), allocatable ::var 
   integer, dimension(:), allocatable ::  timestep
   
   integer :: nlat, nlon, ntile, ntime, nvar, nmp
   
   real, dimension(:,:,:,:), allocatable ::newvar 

   !--- vars read by subr. read_args -also  from file input.dat
   !--- Nvars = # vars contained in binary file output by host program 
   !--- dimx = typically # landpoints over which the var is specified at each timestep 
   !--- dimy = # timesteps
   integer :: dimx, dimy 
     
   integer :: i,j,k,l,m

   real :: lon_dx
      
      !single site
      nlat =1
      nlon =1 

      !======================================================================!
      !=== read perl script interp. (input.dat) of command line args      ===!
      !--- which determine behaviour of program. which file to process,   ===!
      !=== plot/write text file, how to smooth the data                   ===! 
      !======================================================================!
      call read_args(file_var)

         print *, "Arguments file has been read"
      !======================================================================!
      !=== 
      !--- 
      !=== 
      !======================================================================!
         call read_file( trim(file_var), varname, dimx, dimy)
         allocate( var(dimy,dimx) )
         var = 0.
         call read_dat_file( trim(file_var), var, dimx, dimy)

print *, 'jhan', shape( var)         
         print *, "Variable file has been read"

         
         stop
!print *, var

!         nvar = dimx
!         ntime = dimy
!         allocate( timestep(ntime) )
!         do i=1,ntime
!            timestep(i) = i
!         enddo        
!
!         allocate ( newvar(nlon, nlat, ntile, ntime) ) 
!         
!         newvar = -2000000000. 
!
!         lat_j = 0
!         lon_k = 0
!
!         do i=1,nmp
!            do j=1,nlat 
!               if(lat_index(i) == lat(j)) then
!                  lat_j(i) = j
!                  exit   
!               endif           
!            enddo
!            do k=1,nlon 
!               if(lon_index(i) == lon(k)) then
!                  lon_k(i) = k
!                  exit   
!               endif           
!           enddo
!            if(lat_j(i) ==0) print *, i, "lat failed to tag"
!            if(lon_k(i) ==0) print *, i, "lon failed to tag"
!         enddo
!
!         do i=1,nmp
!            !newvar( lon_k(i), lat_j(i), int( tile_index(i) ), : )  = real( var(:,i) )         
!            newvar( lon_k(i), lat_j(i), int( tile_index(i) ), : )  = real( tile_frac(i) * var(:,i) )         
!         enddo
! 
!         lon_dx = 360./(nlon)
!
!         !lon(1) = 0.
!         !do i=2,nlon
!         !   lon(i) = lon(i-1) + lon_dx
!         !enddo
!      
!       
!         print *, nlat, nlon, ntile, ntime, nvar!, real(lat), &
!                !real(lon)
!        !stop 
!         call write_ncdf_template( trim(trim(varname)//".nc"), nlat, nlon, ntile, ntime, nvar, real(lat), &
!                real(lon), tile, timestep, newvar, varname )
!          
   stop
end program debug 
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!




!=============================================================================!
!=== subr. to read command line args interpreted by perl script.comm. line ===!
!=== see top description of program for further explanation                ===! 
!=============================================================================!
   
subroutine read_args(file_var)
   implicit none
   character(len=30), intent(out) :: file_var
   integer, parameter :: gok=0
   integer :: gopenstatus
   
   open(unit=1,file='input.dat', status="unknown",action="read", iostat=gopenstatus )
      if(gopenstatus==gok) then
         read(1,*), file_var
      else
         stop 'input.dat NOT found to read'
      endif
   close(1)
   
   return
end subroutine read_args
   
!=======================================================================!
!=======================================================================!






