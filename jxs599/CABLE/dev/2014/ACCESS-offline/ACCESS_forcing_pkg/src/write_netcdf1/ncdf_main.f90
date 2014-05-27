!are we writing the template
#define TEMPLATE 0 

!--- base of filename read by subr. read_args from file input.dat
!--- which is created by script cable_diag.pl after reading files created
!--- host program (i.e. UM, Mk3L, offline CABLE)
!===================================================================
!===DIMENSIONS
!===================================================================

!---template only : LATITUDE (analagous to x-axis)
!                   LONGITUDE (analagous to y-axis)
!                   TILE_DEF (analagous to z-axis)
!                   TIME_DEF (analagous to time-axis)
#define LATITUDE trim(dir_catted)//"latitude"
#define LONGITUDE trim(dir_catted)//"longitude"
#define NTILE_DEF 17 
#define TILE_DEF (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17/) 

!---template only per consistent time length: else read and set 
#define TIME_DEF "time"  


!===================================================================
!===INDEXES 
!===================================================================

!---template only : LAT_INDEX - latitude of mp point
!                   LON_INDEX - longitude of mp point
!                   TILE_INDEX - tile of mp point
#define LAT_INDEX trim(dir_catted)//"lat_index"
#define LON_INDEX trim(dir_catted)//"lon_index"
#define TILE_INDEX trim(dir_catted)//"tile_index"




program debug 
   use debug_read_mod   !only called funcs are public in this mod 
   use debug_write_mod  !only called funcs are public in this modi

   implicit none

   character(len=2000) :: dir_catted 
   character(len=30) :: file_var, varname, newfile

   real*8, dimension(:,:), allocatable ::var 
   real*8, dimension(:), allocatable :: lat, lon
   real*8, dimension(:), allocatable :: lat_index, lon_index, tile_index 
   
   integer, dimension(:), allocatable ::  timestep, tile
   
   integer :: nlat, nlon, ntile, ntime, nvar, nmp
   
   real, dimension(:,:,:,:), allocatable ::newvar 

   !--- vars read by subr. read_args -also  from file input.dat
   !--- Nvars = # vars contained in binary file output by host program 
   !--- dimx = typically # landpoints over which the var is specified at each timestep 
   !--- dimy = # timesteps
   integer :: dimx, dimy 
     
   integer :: i,j,k,l,m, logu
   integer, dimension(:), allocatable :: lon_k, lat_j

   real :: lon_dx
   character(len=1) :: dummy 
   character(len=*), parameter :: logfile = '/home/599/jxs599/ncdf.log'

     ! get CLI args, basename of file & how many nodes
      IF( IARGC() > 0 ) THEN
         CALL GETARG(1, dir_catted)
      ENDIF
      
      logu=444 
      open(unit=logu,file=logfile, status="unknown",action="write" )
      write(logu,*) "Fortran executable: write netcdf"
      write(logu,*) "called with args: " 
      write(logu,*) "diectory of mapping data: ",trim( dir_catted )
      write(logu,*) "executing.... " 



      !======================================================================!
      !=== read perl script interp. (input.dat) of command line args      ===!
      !--- which determine behaviour of program. which file to process,   ===!
      !=== plot/write text file, how to smooth the data                   ===! 
      !======================================================================!
      call read_args(file_var, newfile)

      !======================================================================!
      !=== read info about the spec. binary data which was created by the ===!
      !--- host so we know how many vars are contained within, how many   ===!
      !=== points there are at each timestep, how many timesteps.         ===! 
      !======================================================================!
         ! latitude
         call read_file( trim(LATITUDE), varname, dimx, dimy)
         allocate( lat(dimx) )
         call read_dat_file( trim(LATITUDE), lat, dimx, dimy)
         nlat = dimx
         
         ! longitude
         call read_file( trim(LONGITUDE), varname, dimx, dimy)
         allocate( lon(dimx) )
         call read_dat_file( trim(LONGITUDE), lon, dimx, dimy)
         nlon = dimx

         ! latitude index
         call read_file( trim(LAT_INDEX), varname, dimx, dimy)
         allocate( lat_index(dimx) )
         call read_dat_file( trim(LAT_INDEX), lat_index, dimx, dimy)
         nmp = dimx
         allocate( lon_k(nmp), lat_j(nmp) )
            
         ! longitude index
         call read_file( trim(LON_INDEX), varname, dimx, dimy)
         allocate( lon_index(dimx) )
         call read_dat_file( trim(LON_INDEX), lon_index, dimx, dimy)

         ! tile index
         call read_file( trim(TILE_INDEX), varname, dimx, dimy)
         allocate( tile_index(dimx) )
         call read_dat_file( trim(TILE_INDEX), tile_index, dimx, dimy)
         ntile = NTILE_DEF 
         allocate( tile(ntile) )
         tile = TILE_DEF
         
         ! variable
         !#@print *, ""
         !#@print *, ""
         !#@print *, ""
         !#@print *, "tile_frac ",  trim(trim(dir_catted)//trim(file_var))
         !#@print *, ""
         !#@print *, ""
         !#@print *, ""
         !#@stop
         call read_file( trim(trim(dir_catted)//trim(file_var)), varname, dimx, dimy)
         allocate( var(dimy,dimx) )
         call read_dat_file( trim(trim(dir_catted)//trim(file_var)), var, dimx, dimy)
         nvar = dimx
         ntime = dimy
         allocate( timestep(ntime) )
         do i=1,ntime
            timestep(i) = i
         enddo        
         
         ! allocate mem for variable to netcdf
         allocate ( newvar(nlon, nlat, ntile, ntime) ) 
         
         ! initialize, also serves as flag value 
         newvar = 0./0. 
         lat_j = 0
         lon_k = 0
        
         ! index mp patches with corresponding lat/long
         ! lat_index is mp long gives index  into lat array
         do i=1,nmp
            do j=1,nlat 
               if(lat_index(i) == lat(j)) then
                  lat_j(i) = j
                  exit   
               endif           
            enddo
            do k=1,nlon 
               !print *, "lon_index(i), lon(k)", lon_index(i), lon(k)
               if(lon_index(i) == lon(k)) then
                  lon_k(i) = k
                  exit   
               endif           
            enddo
            if(lat_j(i) ==0) print *, i, "lat failed to tag"
            if(lon_k(i) ==0) print *, i, "lon failed to tag"
         enddo
!stop
         !do i=1,nmp
         !   print *, "i ", i
         !   print *, "lat_j ", lat_j(i)
         !   print *, "lon_k ", lon_k(i)
         !   !read(*,*) dummy
         !enddo
   !print *, "lat_j", lat_j
   !print *, "lon_k", lon_k

         do i=1,nmp
            newvar( lon_k(i), lat_j(i), int( tile_index(i) ), : )  = real( var(:,i) )
         enddo
   !print *, "newvar ", newvar
         ! as recorded longitude goes from 0 to 180 and back to zero 
         lon_dx = 360./(nlon)
         !print *,nlon, lon_dx
         lon(1) = 0.
         do i=2,nlon
            lon(i) = lon(i-1) + lon_dx
         enddo
         !print *,lon
       
         call write_ncdf_template( newfile, nlat, nlon, ntile, ntime, nvar,    &
                                   real(lat),                                  & 
                                   real(lon), tile, timestep, newvar )
          
      write(logu,*) "Fortran executable: write netcdf"
      write(logu,*) "Finished." 
      close(logu)
   stop
end program debug 
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!




!=============================================================================!
!=== subr. to read command line args interpreted by perl script.comm. line ===!
!=== see top description of program for further explanation                ===! 
!=============================================================================!
   
subroutine read_args(file_var, newfile)
   implicit none
   character(len=30), intent(out) :: file_var, newfile
   integer, parameter :: gok=0
   integer :: gopenstatus
   
   open(unit=1,file='input.dat', status="unknown",action="read", iostat=gopenstatus )
      if(gopenstatus==gok) then
         read(1,*), file_var
         read(1,*), newfile
      else
         stop 'input.dat NOT found to read'
      endif
   close(1)
   
   return
end subroutine read_args
   
!=======================================================================!
!=======================================================================!






