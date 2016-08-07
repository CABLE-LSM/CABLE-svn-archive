module cable_routing


  use cable_common_module, only : cable_user,ktau_gl
  use cable_def_types_mod, only : r_2, ms, mp,mlat,mlon,mland,&
                                  soil_snow_type,soil_parameter_type,        &
                                  veg_parameter_type, canopy_type, met_type

  use cable_IO_vars_module, only : landpt, patch, max_vegpatches, metGrid,  &
                                  land_x,land_y,latitude,longitude,         &
                                  lat_all, lon_all,                         &
                                  logn,output,xdimsize,ydimsize,check,mask
  use cable_rrm_nc_names
  use netcdf

  implicit none


  public
  !**************************************************************************!
  !  Temporary to avoid having to link to the CABLE mods while testing       !
  !**************************************************************************!  
  
  !parameters
  integer,   parameter :: n_river_tracers = 1   !number of tracers...liq  -->carbon,nitrogen and ice??.  not in use as of now
  integer,   parameter :: max_n_ovrlap = 16!324   !assume at most 1x1 to 0.0625x0.0625 -->16x16  do 18x18=324 to be safe
  real(r_2), parameter :: re = 6371000.0              !radius of the earth (km)
  real(r_2), parameter :: deg2rad = 3.14159/180.0     !constant converts degrees to radians
  real(r_2), parameter :: eps=1e-5                    !tolerance parameter.  UNUSED NOW 

  real(r_2), parameter :: river_theta_param = 0.5
  real(r_2), parameter :: source_area_channel_cutoff = 50000.0

  real(r_2), parameter :: fNaN = -1e36
  integer, parameter   :: iNaN = -999999

  real(r_2)  :: del_river = 600.0
  
  type river_grid_type    !will move to cable types eventually    
  
    integer                            :: nlat
    integer                            :: nlon
    integer                            :: npts
    integer                            :: nbasins
    integer                            :: nrr_cells   !number of active cells
    
    !move below into a separate data container?? nah...
    real(r_2), allocatable, dimension(:) :: lat   !1dimensional array stored as lat=0,lon={1,nlon}; lat=1,lon={1,nlon}, etc.  compressed doesn't contain all lat/lon
    real(r_2), allocatable, dimension(:) :: lon
    real(r_2), allocatable, dimension(:) :: length
    real(r_2), allocatable, dimension(:) :: slope
    real(r_2), allocatable, dimension(:) :: elev
    real(r_2), allocatable, dimension(:) :: area
    real(r_2), allocatable, dimension(:) :: source_area     !upstream source area draining into grid cell
    integer,   allocatable, dimension(:) :: land_mask       !1=land,0=ocean,2=land bordering ocean
    integer,   allocatable, dimension(:) :: dwnstrm_index   !index of cell flow goes towards
    integer,   allocatable, dimension(:) :: ocean_outlet    !index of the ending grid cell
    integer,   allocatable, dimension(:) :: upstrm_number   !number of grid cells upstream from the current
    integer,   allocatable, dimension(:) :: active_cell     !1=cell is active, 0=cell not active
    integer,   allocatable, dimension(:) :: direction       !direction of the river flow
    integer,   allocatable, dimension(:) :: is_main_channel !main river channel = 1, not main channel = 0
    integer,   allocatable, dimension(:) :: orig_ind        !the index of the original global array that has all lat/lon points

    
  end type river_grid_type
  
  
  type river_flow_type     !one global.  make 2D and have more than one tracer??????
  
    real(r_2), allocatable, dimension(:) :: mass        !river water mass  (total.  not per m2 like the lsm)
    real(r_2), allocatable, dimension(:) :: vel         !river water speed
    real(r_2), allocatable, dimension(:) :: hgt         !river water height 
    real(r_2), allocatable, dimension(:) :: srf_runoff_flux !surface runoff flux (+ to river) from lsm
    real(r_2), allocatable, dimension(:) :: Fin
    real(r_2), allocatable, dimension(:) :: Fout
    real(r_2), allocatable, dimension(:) :: Fin_n  !flow into cell at new timestep
    real(r_2), allocatable, dimension(:) :: Mass_n !mass in cell at new timestep
    real(r_2), allocatable, dimension(:) :: ocean_outflow  !keep track of outflow to ocean or off of grid
    real(r_2), allocatable, dimension(:) :: theta           !wave speed
    real(r_2)                            ::  total_outflow 
  end type river_flow_type
  
  !below will go into cable_routing_main_routine.  !global on myrank =0, local on myrank=1->nprocs
  type(river_grid_type), save :: river_grid 
  type(river_flow_type), save :: river_var
  

  interface create
    module procedure alloc_river_flow ! interface for a module
    module procedure alloc_river_grid ! procedure is implicit
  end interface

  interface destroy
    module procedure dealloc_river_flow ! interface for a module
    module procedure dealloc_river_grid ! procedure is implicit
  end interface
  
  interface read_nc
    module procedure read_nc_flt_1d
    module procedure read_nc_int_1d
    module procedure read_nc_flt_2d
    module procedure read_nc_int_2d
  end interface
 

contains


 subroutine run_river_route_model(soil,ssnow,river_filename,is_global,dels_river,dels)
     implicit none
     type(soil_snow_type), intent(inout) :: ssnow
     type(soil_parameter_type), intent(in) :: soil
     character(len=99), intent(in) :: river_filename
     logical, intent(in)       :: is_global
     real, intent(in)          :: dels_river,dels

     integer :: npts,nsub_steps,i

     !npts = size(soil%hksat(:,1),dim=1)
     npts = mland  !no tiling for now

     if (ktau_gl .le. 1) then  !initialize and read data
        call create(river_var,npts)
        call create(river_grid,npts)
        call init_river_vars(river_var)
        call get_river_route_data(river_grid,river_filename)
        call find_downstream_index(river_grid,is_global)
        call find_main_river_channels(river_grid)
     else

       river_var%mass(:) = ssnow%river_mass(:)*river_grid%area(:)
     end if

     del_river = real(dels_river,r_2)

     call map_lsm_flux_to_river(river_var,river_grid,ssnow,dels)

     !substep time for stability
     nsub_steps = int(dels/dels_river)

     river_var%ocean_outflow(:) = 0._r_2  !ocean for global or else at edge of grid 

     do i=1,nsub_steps
        call step_river_srf_subsrf_routing_kinematic(river_var,river_grid)
     end do

     call store_river_into_ssnow(river_var,river_grid,ssnow)

  end subroutine run_river_route_model


  subroutine map_lsm_flux_to_river(river,river_grid,ssnow,dels)
     implicit none
     type(river_flow_type), intent(inout) :: river
     type(river_grid_type), intent(in)   :: river_grid
     type(soil_snow_type), intent(in)    :: ssnow
     real, intent(in)                    :: dels

     river%srf_runoff_flux(:) = river_grid%area(:)*(ssnow%rnof1+ssnow%rnof2)  ! in kg/s

  end subroutine map_lsm_flux_to_river

  subroutine store_river_into_ssnow(river,river_grid,ssnow)
     type(river_flow_type), intent(inout) :: river
     type(river_grid_type), intent(in) :: river_grid
     type(soil_snow_type), intent(inout)  :: ssnow

     ssnow%river_mass(:) = river%mass(:)/river_grid%area(:)

  end subroutine store_river_into_ssnow

  subroutine init_river_vars(river)
     implicit none
     type(river_flow_type), intent(inout) :: river

     river%Fin(:)   = 0._r_2
     river%Fout(:)  = 0._r_2
     river%vel(:)   = 0._r_2
     river%mass(:)  = 100._r_2

  end subroutine init_river_vars


!----------------------------------------------------------------------------!

  subroutine get_river_route_data(grid_var,filename)
  !i am assuming the river route grid matches that of the model, 
  !not really valid if doing regional from a global dataset
  !
  
  !reads while file.  need to determine if covered by lsm grid later
    implicit none

    type(river_grid_type), intent(inout)    :: grid_var    
    character(len=99), intent(in)          :: filename
    
    integer :: ncid_river
    integer :: rr_lat_dim_id,rr_lon_dim_id,rr_lat_var_id,rr_lon_var_id
    integer :: nlat_rr_file, nlon_rr_file, npts_rr_file
    
    
    real(r_2), dimension(:)  , allocatable :: lat_rvr
    real(r_2), dimension(:)  , allocatable :: lon_rvr
    
    integer :: i,j,k     !integers to loop through lsm data
    
    integer, dimension(2) :: start_inds
    integer, dimension(2) :: end_inds
    integer :: nc_check
    integer, dimension(:), allocatable :: lon_inds_rvr,lat_inds_rvr

    allocate(lat_inds_rvr(mland))
    allocate(lon_inds_rvr(mland))
    
    call check_nc( nf90_open(trim(filename), nf90_nowrite, ncid_river) )

    !get dim ids
    write(*,*) trim(rr_lat_dim_name)
    write(*,*) trim(rr_lon_dim_name)

    call check_nc( nf90_inq_dimid(ncid_river, rr_lat_dim_name, rr_lat_dim_id) )
    call check_nc( nf90_inq_dimid(ncid_river, rr_lon_dim_name, rr_lon_dim_id) )
    
    !get dim lengths
    call check_nc( nf90_inquire_dimension(ncid_river,rr_lat_dim_id,len=nlat_rr_file) )
    call check_nc( nf90_inquire_dimension(ncid_river,rr_lon_dim_id,len=nlon_rr_file) )
    
    allocate(lat_rvr(nlat_rr_file))
    allocate(lon_rvr(nlon_rr_file))
    
    !read the latitude and longitudes
    call check_nc( nf90_inq_varid(ncid_river,rr_lat_var_name,rr_lat_var_id) )
    call check_nc( nf90_inq_varid(ncid_river,rr_lon_var_name,rr_lon_var_id) )
    
    !
    call check_nc( nf90_get_var(ncid_river,rr_lat_var_id,lat_rvr(:)) )
    call check_nc( nf90_get_var(ncid_river,rr_lon_var_id,lon_rvr(:)) )
    
    start_inds = (/1           , 1          /)
    end_inds   = (/nlon_rr_file,nlat_rr_file/)
    
    npts_rr_file = nlon_rr_file * nlat_rr_file

    !find indieces that match the lsm grid
    !note this assumes that they share same resolution!
    call find_matching_rvr_cells(lon_rvr,lat_rvr,lon_inds_rvr,lat_inds_rvr)
    
    !allocate variable for the river grid
    call alloc_river_grid(grid_var,mland)  !allocate for the lsm grid, not native river data grid
    
    grid_var%npts = mland  !includes all lat/lon even ocean
    grid_var%nlat = mlat
    grid_var%nlon = mlon
    !fill lat/lon in grid_var variable
    !do j=1,mlat
    !  do i=1,mlon
    !    k = (j-1)*nlon_rr_file + i
    !    grid_var%lat(k) = lat_rvr(j)
    !    grid_var%lon(k) = lon_rvr(i)
    !    grid_var%orig_ind(k) = k
    !  end do
    !end do
    grid_var%lat(:) = latitude(:)
    grid_var%lon(:) = longitude(:)
    
    !read the data
    !call read_nc(ncid_river,mask_name    ,start_inds,end_inds,grid_var%land_mask(:),lon_inds_rvr,lat_inds_rvr)
    call read_nc(ncid_river,rdir_name    ,start_inds,end_inds,grid_var%direction(:),lon_inds_rvr,lat_inds_rvr)
    call read_nc(ncid_river,length_name  ,start_inds,end_inds,grid_var%length(:),lon_inds_rvr,lat_inds_rvr)
    !call read_nc(ncid_river,slope_name,start_inds,end_inds,grid_var%slope(:),lon_inds_rvr,lat_inds_rvr)
    !call read_nc(ncid_river,elev_name,start_inds,end_inds,grid_var%elev(:),lon_inds_rvr,lat_inds_rvr)
    call read_nc(ncid_river,src_area_name,start_inds,end_inds,grid_var%source_area(:),lon_inds_rvr,lat_inds_rvr)
    
    nc_check = nf90_close(ncid_river)

    deallocate(lat_inds_rvr)
    deallocate(lon_inds_rvr)
     
  contains
    subroutine find_matching_rvr_cells(lon_rvr,lat_rvr,lon_inds_rvr,lat_inds_rvr)

      real, dimension(:), intent(in) :: lon_rvr
      real, dimension(:), intent(in) :: lat_rvr
      integer, dimension(:), intent(out) :: lon_inds_rvr                                       
      integer, dimension(:), intent(out) :: lat_inds_rvr                                       

      integer :: nlon_rvr,nlat_rvr
      real :: distance,newLength
      integer :: k,i,j

      nlon_rvr = size(lon_rvr,dim=1) 
      nlat_rvr = size(lat_rvr,dim=1) 

      ! and longitude(:) has already been converted to -180 to 180 for CABLE.
      lon_inds_rvr(:) = -999 
      lat_inds_rvr(:) = -999 
      do k = 1, mland
        distance = 300.0 ! initialise, units are degrees
        do j = 1, nlat_rvr 
        do i = 1, nlon_rvr
            newLength = SQRT((lon_rvr(i) - longitude(k))**2      &
                           + (lat_rvr(j) -  latitude(k))**2)
            if (newLength < distance) then
              distance = newLength
              lon_inds_rvr(k) = i 
              lat_inds_rvr(k) = j 
            end if
        end do
        end do
        if (lon_inds_rvr(k) < -900 .or. lat_inds_rvr(k) < -900) then
          write(*,*) 'Land point ', k, ' cannot find the nearest grid!'
          write(*,*) 'lon, lat = ', longitude(k), latitude(k)
          stop
        end if
     end do

    end subroutine find_matching_rvr_cells
    

  end subroutine get_river_route_data

!-----------------------------------------------------------------------------
!             !current data set doesn't use these numbers as it is: 1,2,4,8,16,32,64,128,256
!                       32 64  128    1.0  lat(3)
!                       16     1      0.0  lat(2)
!                        8  4  2     -1.0  lat(1)

  integer function dirc2latindex(in_dirc)   !assuming lat goes from south to north....
    implicit none
    integer, intent(in) :: in_dirc

    if (in_dirc .eq. 2 .or. in_dirc .eq. 4 .or. in_dirc .eq. 8) then
       dirc2latindex = -1
    elseif (in_dirc .eq. 32 .or. in_dirc .eq. 64 .or. in_dirc .eq. 128) then
       dirc2latindex = 1
    else
       dirc2latindex = 0
    end if

  end function dirc2latindex
  
!----------------------------------------------------------------------------!
!                      179 180   181
!                     l(1) l(2) l(3)
!                       32 64  128    1.0  lat(3)
!                       16     1      0.0  lat(2)
!                        8  4  2     -1.0  lat(1)

  integer function dirc2lonindex(in_dirc)   !lon goes from west to east
    implicit none
    integer, intent(in) :: in_dirc

    dirc2lonindex = 0
    if (in_dirc .eq. 1 .or. in_dirc .eq. 2 .or. in_dirc .eq. 128) then
      dirc2lonindex  = 1
    elseif (in_dirc .eq. 8 .or. in_dirc .eq. 16 .or. in_dirc .eq. 32) then
      dirc2lonindex = -1
    else
      dirc2lonindex = 0
    end if

  end function dirc2lonindex

!----------------------------------------------------------------------------!
!  call find_downstream_index(grid_var%dwnstrm_index,river_dirc)  !this only works if it is on the global grid

  subroutine find_downstream_index(grid_var,is_global)
  
    implicit none
    
    type(river_grid_type), intent(inout) :: grid_var
    logical,               intent(in   ) :: is_global
    
    integer :: i,j,k,ii,jj,kk,is
    integer :: lon_wrap
    logical :: off_grid,keep_search

    !make sure we called this prior to reordering or removing grid points
!    if (grid_var%npts .ne. grid_var%nlat*grid_var%nlon) &
!         stop "Must find the downstream index using the entire global grid not a subsection"

!    if (is_global) then
!      lon_wrap = grid_var%nlon
!    else
!      lon_wrap = 1
!    end if
    
    grid_var%dwnstrm_index(:) = 0
    
    do k=1,mland
        off_grid = .false.
      
        if (grid_var%direction(k) /= -9999) then
        
          ii = land_x(k) + dirc2latindex(grid_var%direction(k))
          jj = land_y(k) + dirc2lonindex(grid_var%direction(k))
        
          if (is_global) then
            if (ii .lt. 1 )            ii = ii + grid_var%nlon  
            if (ii .gt. grid_var%nlon) ii = ii - grid_var%nlon

            if (jj .lt. 1) jj = 1
            if (jj .gt. grid_var%nlat) jj = grid_var%nlat

          else
             if ((ii .lt. 1) .or. (ii .gt. grid_var%nlon) .or. &
                 (jj .lt. 1) .or. (jj .gt. grid_var%nlat)) then
                off_grid = .true.
             end if 
          end if 

          if (jj .lt. 1 .or. jj .gt. grid_var%nlat .or. ii .lt. 1 .or. ii .gt. grid_var%nlon) then
            write(*,*) ii
            write(*,*) jj
            stop "needs checks yo"
          endif

          !need to search the land_x and land_y arrays to find vaue of ii and jj
          kk = -1
          is = 1
          keep_search = .true.
          do while (keep_search)
             if ((land_x(is) .eq. ii) .and. (land_y(is) .eq. jj)) then
               kk = is
               keep_search = .false.
             else
               is = is + 1
             end if
             if (is .gt. mland) keep_search = .false.
          end do

          
          if (off_grid) kk = -1

          grid_var%dwnstrm_index(k) = kk

        endif
    enddo
    

  end subroutine find_downstream_index
  
    
!----------------------------------------------------------------------------! 
  subroutine find_main_river_channels(grid_var)
    implicit none

    type(river_grid_type), intent(inout) :: grid_var

    integer :: k,kk,j
    integer :: ntot
    logical :: keep_looping

    ntot = size(grid_var%source_area(:),dim=1)


    do k=1,ntot
      if (grid_var%source_area(k) .ge. source_area_channel_cutoff) then
        grid_var%is_main_channel = 1
      else
        grid_var%is_main_channel = 0
      end if
    end do

  end subroutine find_main_river_channels

!-----------------------------------------------------------------------------   
  subroutine step_river_routing(river,grid_var,dels)
     implicit none
     
    type(river_flow_type),  intent(inout) :: river   !contains mass,flow variables
    type(river_grid_type),  intent(in)    :: grid_var
    real,                   intent(in)    :: dels
     
    integer :: kk_begind, kk_endind
    integer :: i,j,ii,kk
       
      kk_begind = 1
      kk_endind = grid_var%nrr_cells
    

      river%Fin(kk_begind:kk_endind) = 0._r_2   !zero out input fluxes

      do kk = kk_begind, kk_endind

        river%Fout(kk) = river%mass(kk)*river%vel(kk)/grid_var%length(kk)

         if ((grid_var%dwnstrm_index(kk)) .gt. 1 .and. (grid_var%dwnstrm_index(kk) .le. grid_var%nrr_cells)) then
            river%Fin(grid_var%dwnstrm_index(kk)) = river%Fin(grid_var%dwnstrm_index(kk)) + river%Fout(kk)
         else
            river%total_outflow = river%total_outflow  + river%Fout(kk)
         end if

      end do

      river%mass(kk_begind:kk_endind) =  river%mass(kk_begind:kk_endind) - &
                                        (river%Fout(kk_begind:kk_endind) + river%Fin(kk_begind:kk_endind))*dels
                                        
    end subroutine step_river_routing
    
!----------------------------------------------------------------------------!

  subroutine step_river_srf_subsrf_routing_kinematic(river,grid_var)
    implicit none
     
    type(river_flow_type),          intent(inout) :: river   !contains mass,flow variables
    type(river_grid_type),          intent(in)    :: grid_var
     
    integer :: kk_begind, kk_endind
    integer :: i,j,k,ii,kk
    
    
      kk_begind = 1                !loop over all points for this mpi task.
      kk_endind = size(river%mass(:),dim=1)
           
      river%Fin_n(:) = 0._r_2
      river%Mass_n(:) = 0._r_2
    
      do kk = kk_begind, kk_endind

        river%theta(kk) = river_theta_param*del_river / grid_var%length(kk)
      
        k = grid_var%dwnstrm_index(kk)

        if (k .gt. 0) then
        
           !  main channel flow.  use manning and a rectangular channel eventually.  kinematic place holder for now
           !compute new river mass at current cell
          river%Mass_n(kk) = (1.-river%theta(kk)) * river%mass(kk) + river%Fin(kk) +&
                             river%srf_runoff_flux(kk) * del_river
          !compute flux to the downstream cell
          river%Fin_n(k)  = river%Fin_n(k) + river%theta(kk) * river%Mass_n(kk)          

        else  !outlet is off of the grid
           river%Mass_n(kk) = (1.-river%theta(kk)) * river%mass(kk) + river%Fin(kk) +&
                             river%srf_runoff_flux(kk) * del_river
           !keep track out total outflow from all cells
           river%total_outflow = river%total_outflow + river%theta(kk) * river%Mass_n(kk)
           !ocean outflow for this cell in kg/s
           river%ocean_outflow(kk) = river%ocean_outflow(kk) + river%theta(kk) * river%Mass_n(kk) / del_river
        end if

      end do
      
      do kk = kk_begind, kk_endind
      
        river%Fin(kk)  = river%Fin_n(kk)  !store new values
        
        river%vel(kk) = river%theta(kk) / del_river * river%mass(kk)   !should be mass flux not vel
        river%mass(kk) = river%Mass_n(kk)
      
     end do                        

     
   end subroutine step_river_srf_subsrf_routing_kinematic
    
!----------------------------------------------------------------------------!

!**********Below are the Routines For Allocating/Deallocating****************!

  subroutine alloc_river_flow(var,npts)
    implicit none
    type(river_flow_type),intent(inout)  :: var
    integer, intent(in) :: npts

    allocate(var%mass(npts))
    var%mass(:) = fNaN
    
    allocate(var%vel(npts))
    var%vel(:) = fNaN
    
    allocate(var%hgt(npts))
    var%hgt(:) = fNaN    
    
    allocate(var%srf_runoff_flux(npts))
    var%srf_runoff_flux(:) = fNaN
    
    allocate(var%Fin(npts))
    var%Fin(:) = fNaN
    
    allocate(var%Fout(npts))
    var%Fout(:) = fNaN       

    allocate(var%Mass_n(npts))
    var%Mass_n(:) = fNaN

    allocate(var%Fin_n(npts))
    var%Fin_n(:) = fNaN

    allocate(var%ocean_outflow(npts))
    var%ocean_outflow(:) = fNaN

    var%total_outflow = 0.

    allocate(var%theta(npts))
    var%theta(:) = fNaN

  end subroutine alloc_river_flow
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_flow(var)
    implicit none
    type(river_flow_type),intent(inout) :: var

    deallocate(var%mass)
    deallocate(var%vel)
    deallocate(var%hgt)
    deallocate(var%srf_runoff_flux)
    deallocate(var%Fin)
    deallocate(var%Fout)
    deallocate(var%Mass_n)
    deallocate(var%Fin_n)
    deallocate(var%ocean_outflow)
    deallocate(var%theta)

  end subroutine dealloc_river_flow
  
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!

  subroutine alloc_river_grid(var,npts)
    implicit none
    type(river_grid_type), intent(inout) :: var
    integer, intent(in)                  :: npts

    allocate(var%lat(npts))
    var%lat(:) = fNaN
    allocate(var%lon(npts))
    var%lon(:) = fNaN
    allocate(var%length(npts))
    var%length(:) = fNaN
    allocate(var%slope(npts))
    var%slope(:) = fNaN   
    allocate(var%dwnstrm_index(npts))
    var%dwnstrm_index(:) = iNaN  
    allocate(var%upstrm_number(npts))
    var%upstrm_number(:) = iNan
    allocate(var%ocean_outlet(npts))
    var%ocean_outlet(:) = iNaN
    allocate(var%elev(npts))
    var%elev(:) = fNaN
    allocate(var%land_mask(npts))
    var%land_mask(:) = iNaN
    allocate(var%active_cell(npts))
    var%active_cell(:) = iNaN
    allocate(var%direction(npts))
    var%direction(:) = iNaN
    allocate(var%area(npts))
    var%area(:) = fNaN
    allocate(var%is_main_channel(npts))
    var%is_main_channel(:) = iNaN
    allocate(var%source_area(npts))
    var%source_area(:) = fNaN
    allocate(var%orig_ind(npts))
    var%orig_ind(:) = iNaN
    !since we know the total number of points set it
    var%npts = npts


  end subroutine alloc_river_grid
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_grid(var)
    implicit none
    type(river_grid_type), intent(inout) :: var

    write(*,*) 'starting to dealloc grid var'
    deallocate(var%lat)
    deallocate(var%lon)
    deallocate(var%slope)
    deallocate(var%length)
    deallocate(var%dwnstrm_index)
    deallocate(var%elev)
    deallocate(var%land_mask)
    deallocate(var%active_cell)
    deallocate(var%direction)
    deallocate(var%area)
    deallocate(var%is_main_channel)
    deallocate(var%source_area)
    deallocate(var%orig_ind)
!    deallocate(var%topo_ind)
!    deallocate(var%basin_ind)
    write(*,*) 'dealloced gird var now dealloc maps'
    write(*,*) 'dealloced grid'
  end subroutine dealloc_river_grid

!----------------------------------------------------------------------------!


!**********Below are the Routines For Reading Netcdf Files*******************!

  subroutine check_nc(nc_status)
    implicit none
    integer, intent(in) :: nc_status

    if (nc_status .ne. NF90_NOERR) then
      write(*,*)  trim(nf90_strerror(nc_status))
      stop "Stopped in River Routing"
    end if

  end subroutine check_nc   

!----------------------------------------------------------------------------!    
    
  subroutine read_nc_int_1d(ncid,var_name,start_inds,end_inds,var_data,lon_inds_rvr,lat_inds_rvr,mask)
    use netcdf
    implicit none
    
    integer,               intent(in)    :: ncid
    character(len=*)   , intent(in)    :: var_name
    integer, dimension(:), intent(in)    :: start_inds,end_inds
    integer, dimension(:), intent(inout) :: var_data
    integer, dimension(:), intent(in)    :: lon_inds_rvr,lat_inds_rvr
    integer, dimension(:), optional    :: mask
    
    integer :: nc_check,nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    integer, dimension(:,:), allocatable :: int_data
    integer                              :: npts
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1
    npts    = nlon_rr*nlat_rr
    allocate(int_data(nlon_rr,nlat_rr))
    
    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,int_data(:,:),start_inds,end_inds) )

    do k=1,mland
       var_data(k) = int_data(lon_inds_rvr(k),lat_inds_rvr(k))
    end do

    if (present(mask)) then
      where(mask .ge. 1) var_data(:) = -9999
    end if
    
    deallocate(int_data)
     
   end subroutine read_nc_int_1d

!----------------------------------------------------------------------------!   
 
  subroutine read_nc_flt_1d(ncid,var_name,start_inds,end_inds,var_data,lon_inds_rvr,lat_inds_rvr,mask)
    use netcdf
    implicit none
    
    integer,                 intent(in)    :: ncid
    character(len=*)       , intent(in)    :: var_name
    integer  , dimension(:), intent(in)    :: start_inds,end_inds
    real(r_2), dimension(:), intent(inout) :: var_data
    integer, dimension(:), intent(in)    :: lon_inds_rvr,lat_inds_rvr
    integer, dimension(:), optional    :: mask    
    
    integer :: nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    real(r_2), dimension(:,:), allocatable :: flt_data
    integer                                :: npts
    
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1
    npts    = nlon_rr*nlat_rr
    allocate(flt_data(nlon_rr,nlat_rr))
    
    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,flt_data(:,:),start_inds,end_inds) )

    do k=1,mland
       var_data(k) = flt_data(lon_inds_rvr(k),lat_inds_rvr(k))
    end do

    if (present(mask)) then
      where(mask .ge. 1) var_data(:) = -9999
    end if

    deallocate(flt_data)

     
  end subroutine read_nc_flt_1d

!----------------------------------------------------------------------------!   
  
  subroutine read_nc_flt_2d(ncid,var_name,start_inds,end_inds,var_data,mask)
    use netcdf
    implicit none
    
    integer,                 intent(in)    :: ncid
    character(len=*)       , intent(in)    :: var_name
    integer  , dimension(:), intent(in)    :: start_inds,end_inds
    real(r_2), dimension(:,:), intent(inout) :: var_data
    integer, dimension(:,:), optional    :: mask    
    
    integer :: nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    real(r_2), dimension(:,:), allocatable :: flt_data
    integer  , dimension(:)  , allocatable :: mask_1d
    integer                                :: npts
    
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1

    if ((size(var_data,dim=1) .ne. nlat_rr) .or. (size(var_data,dim=2) .ne. nlon_rr)) then
       write(*,*) 'wrong size data passed to read_nc_flt_2d'
       stop
    end if

    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,var_data(:,:),start_inds,end_inds) )

    if (present(mask)) then
      where(mask .ge. 1) var_data(:,:) = -9999.
    end if

     
  end subroutine read_nc_flt_2d

!----------------------------------------------------------------------------!   
  
  subroutine read_nc_int_2d(ncid,var_name,start_inds,end_inds,var_data,mask)
    use netcdf
    implicit none
    
    integer,                 intent(in)    :: ncid
    character(len=*)       , intent(in)    :: var_name
    integer  , dimension(:), intent(in)    :: start_inds,end_inds
    integer, dimension(:,:), intent(inout) :: var_data
    integer, dimension(:,:), optional    :: mask    
    
    integer :: nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    real(r_2), dimension(:,:), allocatable :: flt_data
    integer  , dimension(:)  , allocatable :: mask_1d
    integer                                :: npts
    
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1

    if ((size(var_data,dim=1) .ne. nlat_rr) .or. (size(var_data,dim=2) .ne. nlon_rr)) then
       write(*,*) 'wrong size data passed to read_nc_flt_2d'
       stop
    end if

    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,var_data(:,:),start_inds,end_inds) )

    if (present(mask)) then
      where(mask .ge. 1) var_data(:,:) = -9999.
    end if

     
  end subroutine read_nc_int_2d
  
!----------------------------------------------------------------------------!   
  
!----------------------------------------------------------------------------!   
!****Below are the Routines For MPI settings/gettings/partitioning***********!
!----------------------------------------------------------------------------!  
  

end module cable_routing
