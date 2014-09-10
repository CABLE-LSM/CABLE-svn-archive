module cable_routing

  use netcdf
  !use cable_types
  !use cable_common_module, only : cable_user
  !use cable_def_types_mod, only : r_2, ms, mp,mlat,mlon

  !use cable_IO_vars_module, only : landpt, patch, max_vegpatches, metGrid,  &
  !                                land_x,land_y,latitude,longitude,         &
  !                                lat_all, lon_all,                         &
  !                                logn,output,xdimsize,ydimsize,check,mask


  implicit none
  
  !**************************************************************************!
  !  Temporary to avoid having to link to the CABLE mods while testing       !
  !**************************************************************************!
  integer, parameter :: r_2  = SELECTED_REAL_KIND(12, 50)
  integer   :: mlat = 108
  integer   :: mlon = 243
  integer   :: mp = 243*108
  real(r_2) :: dels = 1800.0
  
  real(r_2), dimension(:), allocatable, save :: latitude
  real(r_2), dimension(:), allocatable, save :: longitude

  real(r_2), dimension(:,:), allocatable, save :: lat_all
  real(r_2), dimension(:,:), allocatable, save :: lon_all   
  !**************************************************************************!
  !  Temporary to avoid having to link to the CABLE mods while testing       !
  !**************************************************************************!  
  
  !parameters
  integer,   parameter :: n_river_tracers = 1   !number of tracers...liq  -->carbon,nitrogen and ice??.  not in use as of now
  integer,   parameter :: max_n_ovrlap = 324   !assume at most 1x1 to 0.0625x0.0625 -->16x16  do 18x18=324 to be safe
  real(r_2), parameter :: re = 6371000.0              !radius of the earth (km)
  real(r_2), parameter :: deg2rad = 3.14159/180.0     !constant converts degrees to radians
  real(r_2), parameter :: eps=1e-5                    !tolerance parameter.  UNUSED NOW 

  real(r_2), parameter :: fNaN = -1e36
  integer, parameter   :: iNaN = -999999
  
  
  type map_grid_type
    
    integer, pointer, dimension(:,:)   :: ind_lgr    !index of the Land cell  Given the River cell
    real(r_2), pointer, dimension(:,:) :: weight_lgr   
    integer, pointer, dimension(:)     :: n_ovrlap_lgr   !number of rivers cells that over lap the land cell  
    
    contains
    procedure :: create  => alloc_river_map_grid
    procedure :: destroy => dealloc_river_map_grid 
        
  end type map_grid_type

  
  type river_grid_type    !will move to cable types eventually    
  
    integer                            :: nlat
    integer                            :: nlon
    integer                            :: npts
    integer                            :: nbasins
    integer                            :: nrr_cells
    
    !move below into a separate data container??
    real(r_2), pointer, dimension(:) :: lat   !1dimensional array stored as lat=0,lon={1,nlon}; lat=1,lon={1,nlon}, etc
    real(r_2), pointer, dimension(:) :: lon
    real(r_2), pointer, dimension(:) :: length
    real(r_2), pointer, dimension(:) :: slope
    real(r_2), pointer, dimension(:) :: elev
    real(r_2), pointer, dimension(:) :: area
    
    integer,   pointer, dimension(:) :: land_mask      !1=land,0=ocean,2=land bordering ocean
    integer,   pointer, dimension(:) :: dwnstrm_index  !index of cell flow goes towards
    integer,   pointer, dimension(:) :: ocean_outlet   !index of the ending grid cell
    integer,   pointer, dimension(:) :: upstrm_number  !number of grid cells upstream from the current
    integer,   pointer, dimension(:) :: active_cell    !1=cell is active, 0=cell not active
    integer,   pointer, dimension(:) :: direction      !direction of the river flow

    type(map_grid_type) :: maps
   
    contains
    procedure :: create             => alloc_river_grid                  !river_grid%create(npts)
    procedure :: destroy            => dealloc_river_grid                !river_gric%destroy()
    
    procedure :: mapit              => determine_hilo_res_mapping        !river_grid%mapit(lsm_latitude,lsm_longitude)
    procedure :: get_data           => get_river_route_data              !river_grid%get_data(filename)
    procedure :: find_dwnstrm       => find_downstream_index             !river_grid%find_dwnstrm()
    procedure :: find_outlets       => associate_ocean_outlet_points     !river_grid%find_outlets()
    procedure :: reorder            => reorder_grid_by_basin             !river_grid%arrange(basins)   !one for grid and one for river_flow?
    procedure :: collapse           => remove_inactive_land_basins       !river_grid%collapse(basins)
       
    
  end type river_grid_type
  
  
  type river_flow_type     !one global.  make 2D and have more than one tracer??????
  
    real(r_2), pointer, dimension(:) :: mass_init   !river water mass at start of timestep (per m2)
    real(r_2), pointer, dimension(:) :: mass        !river water mass  (per m2)
    real(r_2), pointer, dimension(:) :: vel         !river water speed
    real(r_2), pointer, dimension(:) :: hgt         !river water height 
    real(r_2), pointer, dimension(:) :: runoff_flux !total flux (+ to river) from lsm
    real(r_2), pointer, dimension(:) :: Fin
    real(r_2), pointer, dimension(:) :: Fout
    
    contains
    procedure :: create     => alloc_river_flow
    procedure :: destroy    => dealloc_river_flow   
    
    procedure :: step_time  => step_river_routing
    procedure :: balance    => compute_global_mass_balance   !not added yet
    
  end type river_flow_type
  
  type basin_type
  
    integer                        :: begind         !basin index start in reordered global 1d array
    integer                        :: endind         !basin end index in global reordered array
    integer                        :: n_basin_cells  !number of river cells in the basin
    integer, pointer, dimension(:) :: river_points   !the indices for the basin in the unordered global 1d array
    
    contains
    procedure :: create     => alloc_basin
    procedure :: destroy    => dealloc_basin
    
  end type basin_type
  
  type(river_grid_type), TARGET, SAVE :: river_grid
  type(river_flow_type), TARGET, SAVE :: river
  type(basin_type), pointer, save, dimension(:) :: basins
   
  
  !outline
  !if timestep==0 --> 
  !  call river_route_init(river_grid,filename)
  !    call get_river_route_data(filename,river_grid)
  !                                          => grid%{nlat,nlon,npts,lat,lon,direction,length,slope,elevation,land_mask}
  !       !only alloc and fills river data for region covered by the lsm grid.  still need check to ensure that only complete
  !       !basins are included?  this won't work if lsm region ends in the middle of a basin.....this check needs to occur
  !       !after reading whole file and calculating the basins
  !
  !       call alloc_river_route_vars()  !-->call within
  !
  !    call find_downstream_index(river_grid)
  !                            => river_grid%{dwnstrm_index,river_dirc}
  !
  !    call associate_ocean_outlet_points(river_grid)
  !                                       river_grid%{ocean_outlet,upstrm_number,nbasins,nrr_cells,land_mask}
  !
  !    call reorder_grid_by_basin(river,river_grid,basin)
  !                                             => basin%{n_basin_cells,river_points}
  !
  !    call determine_hilo_res_mappings(river_grid,latitude,longitude)
  !                                  => river_grid%{lat,lon,maps%weight_river_to_land,river_grid%maps%ind_river_to_land}
  !
  !    call river_grid%collapse => remove_inactive_basins
  !
  !
  !    call partiion_basins_to_pes(river_grids,basin)
  !                             => river_grids%{nrr_cells,nbasins
  !                             => basin%{n_basin_cells,
  !
  !
  !call map_qrunoff_lsm_to_river(Qsrf,Qsubsrf,river,river_grid)
  !                                        => river%{wat_mass,wat_vol,wat_hgt,wat_vel}
  !                                        => river_grid%maps{ind_river_to_land,weight_river_to_land,n_ovrlap}
  !call step_river(river,river_grid)
  !             => river%{wat_mass,wat_vol,wat_hgt,wat_vel}
  !             => river_grid%{distance,mpi%bg,mpi%ed}
  !   call river_mass_balance()
  !
  !call save_river_output()
  !call map_qriver_river_to_lsm()
  !
  !if timestep==final --> call river_route_end()
  !   call dealloc_river_route_vars()

    !ideas.....    
    !determine runoff input in river cells from lsm:
    !Q_river_runoff(:) = 0._r_2
    !do kk=1,river_grid%npts
    !  do i=1,n_ovrlap_lgr(kk)
    !    Q_river_runoff(kk) = Q_river_runoff(kk)  + Q_runoff(ind_lgr(kk,i))*weight_lgr(kk,i)
    !  end do
    !end do
    
    !map the river water to each grid cell
    !wb_river(:) = 0._r_2
    !do kk=1,river_grid%npts
    !  do i=1,n_ovrlap_lgr(kk)
    !    wb_river_land(ind_lgr(kk,i)) = wb_river_land(ind_lgr(kk,i)) + river%wb(kk)*weight_lgr(kk,i)   !check this
    !  end do
    !end do    


  

contains

!----------------------------------------------------------------------------!
  function compute_global_mass_balance(river_var,grid_var,basin_var) result(mass_error)
    implicit none
    class(river_type)     , intent(inout) :: river_var
    class(river_grid_type), intent(in)    :: grid_var
    class(basin_type)     , intent(in)    :: basin_var
    real(r_2)                             :: mass_error
    
    !local variables
    real(r_2) :: total_mass, total_lsm_flux,total_outflow, init_total_mass
    integer :: i,j,k
    
    init_total_mass = 0._r_2
    total_mass      = 0._r_2
    total_lsm_flux  = 0._r_2
    mass_error      = 0._r_2
    
    total_mass      = sum(river_var%mass(:))
    init_total_mass = sum(river_var%mass_init(:))
    total_lsm_flux  = sum(river_var%runoff_flux(:)) * dels
    !need to compute outflow
    
    total_outflow = 0._r_2
    
    do i=1,grid_var%nbasins
    
      bg = basin_var(i)%begind
      ed = basin_var(i)%endind  
      j = maxloc(grid_var%upstrm_index(bg:ed),dim=1) + bg - 1 !maxloc returns relative to indices
      
      total_outflow = total_outflow + river_var%Fout(j)*grid_var%area(j)*dels
    end do
    
    mass_error = total_mass - init_total_mass + total_lsm_flux - total_outflow
    
  end function compute_global_mass_balance
    


!----------------------------------------------------------------------------!

  subroutine get_river_route_data(river_grid,filename)
  !reads while file.  need to determine if covered by lsm grid later
    use netcdf
    implicit none

    class(river_grid_type), intent(inout)   :: river_grid    
    character(len=250), intent(in)          :: filename
    
    integer :: ncid_river
    integer :: rr_lat_dim_id,rr_lon_dim_id,rr_lat_var_id,rr_lon_var_id
    integer :: nlat_rr_file, nlon_rr_file, npts_rr_file
    integer :: nlat_rr,nlon_rr,npts_rr
    
    
    character(len=*), parameter :: mask_name   = "land_mask"
    character(len=*), parameter :: length_name = "length"
    character(len=*), parameter :: slope_name  = "slope"
    character(len=*), parameter :: elev_name   = "elevation"
    character(len=*), parameter :: rdir_name   = "direction"
    character(len=*), parameter :: rr_lat_dim_name = "lat"
    character(len=*), parameter :: rr_lon_dim_name = "lon"
    character(len=*), parameter :: rr_lat_var_name = "latitude"
    character(len=*), parameter :: rr_lon_var_name = "longitude"
    
    real(r_2), dimension(:)  , allocatable :: lat_data
    real(r_2), dimension(:)  , allocatable :: lon_data
    integer  , dimension(:,:), allocatable :: tmp_mask_data
    
    integer :: i,j,k     !integers to loop through lsm data
    integer :: ri,rj,rk  !integers to loop through river
    integer :: nc_check
    
    integer :: lon_start,lon_end,lat_start,lat_end
    integer, dimension(2) :: start_inds
    integer, dimension(2) :: end_inds
    
    real :: dlat_lsm,dlon_lsm
    
    call check_nc( nf90_open(trim(filename), nf90_nowrite, ncid_river) )

    !get dim ids
    call check_nc( nf90_inq_dimid(ncid_river, rr_lat_dim_name, rr_lat_dim_id) )
    call check_nc( nf90_inq_dimid(ncid_river, rr_lon_dim_name, rr_lon_dim_id) )
    
    !get dim lengths
    call check_nc( nf90_inquire_dimension(ncid_river,rr_lat_dim_id,len=nlat_rr_file) )
    call check_nc( nf90_inquire_dimension(ncid_river,rr_lon_dim_id,len=nlon_rr_file) )
    
    allocate(lat_data(nlat_rr_file))
    allocate(lon_data(nlon_rr_file))
    
    !read the latitude and longitudes
    call check_nc( nf90_inq_varid(ncid_river,rr_lat_var_name,rr_lat_var_id) )
    call check_nc( nf90_inq_varid(ncid_river,rr_lon_var_name,rr_lon_var_id) )
    
    !
    call check_nc( nf90_get_var(ncid_river,rr_lat_var_id,lat_data(:)) )
    call check_nc( nf90_get_var(ncid_river,rr_lon_var_id,lon_data(:)) )
    
    start_inds = (/1           , 1          /)
    end_inds   = (/nlon_rr_file,nlat_rr_file/)
    
    !allocate variable for the river grid
    call river_grid%create(npts_rr_file)
    
    river_grid%npts = npts_rr_file  !includes all lat/lon even ocean
    river_grid%nlat = nlat_rr_file  
    river_grid%nlon = nlon_rr_file
    !fill lat/lon in river_grid variable
    k=0
    do j=1,nlat_rr_file
      do i=1,nlon_rr_file
        k = k + 1
        river_grid%lat(k) = lat_data(j)
        river_grid%lon(k) = lon_data(i)
      end do
    end do
    
    !read the data
    call read_nc_int(ncid_river,mask_name,start_inds,end_inds,river_grid%land_mask(:))
    call read_nc_int(ncid_river,rdir_name,start_inds,end_inds,river_grid%direction(:))
    
    call read_nc_flt(ncid_river,length_name,start_inds,end_inds,river_grid%length(:))
    call read_nc_flt(ncid_river,slope_name ,start_inds,end_inds,river_grid%slope(:))
    call read_nc_flt(ncid_river,elev_name  ,start_inds,end_inds,river_grid%elev(:))
    
    nc_check = nf90_close(ncid_river)
    
    deallocate(tmp_mask_data)
        
  end subroutine get_river_route_data
    
!-----------------------------------------------------------------------------
!             !current data set doesn't use these numbers as it is: 1,2,4,8,16,32,64,128,256
!                       32 64  128
!                       16     1
!                        8  4  2

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

  integer function dirc2lonindex(in_dirc)   !lon goes from west to east
    implicit none
    integer, intent(in) :: in_dirc

    dirc2lonindex = 0
    if (in_dirc .eq. 1 .or. in_dirc .eq. 2 .or. in_dirc .eq. 8) then
      dirc2lonindex  = 1
    elseif (in_dirc .eq. 8 .or. in_dirc .eq. 16 .or. in_dirc .eq. 32) then
      dirc2lonindex = -1
    else
      dirc2lonindex = 0
    end if

  end function dirc2lonindex

!----------------------------------------------------------------------------!
!  call find_downstream_index(river_grid%dwnstrm_index,river_dirc)  !this only works if it is on the global grid

  subroutine find_downstream_index(river_grid)
    implicit none
    
    class(river_grid_type), intent(inout) :: river_grid
    
    integer :: i,j,k,ii,jj,kk
    
    river_grid%dwnstrm_index(:) = 0
    
    do j=1,river_grid%nlat
      do i=1,river_grid%nlon
      
        k = i + (j-1)*river_grid%nlon
        
        if (river_grid%direction(k) /= -9999) then
        
          ii = i + dirc2latindex(river_grid%direction(k))
          jj = j + dirc2lonindex(river_grid%direction(k))
          
          if (ii .lt. 1     )          ii = ii + river_grid%nlon  
          if (ii .gt. river_grid%nlon) ii = ii - river_grid%nlon
          if (jj .lt. 1 .or. jj .gt. river_grid%nlat .or. ii .lt. 1 .or. ii .gt. river_grid%nlon) then
            stop
          endif
          
          kk = ii + (jj-1)*river_grid%nlat
          
          river_grid%dwnstrm_index(k) = kk

        endif
      enddo
    enddo

  end subroutine find_downstream_index
  
!----------------------------------------------------------------------------! 

  subroutine associate_ocean_outlet_points(river_grid)
    implicit none
    
    class(river_grid_type), intent(inout) :: river_grid
  
    !local variables
    integer :: i,j,k
  
    river_grid%ocean_outlet(:)  = -1
    river_grid%upstrm_number(:) = 0
    
    do i=1,river_grid%npts
      j = i
      
      if (river_grid%land_mask(i) .eq. 1) then   !it is a land cell
      
        k = 0
        do while ((river_grid%land_mask(j) .eq. 1) .and. k < river_grid%npts)
           j = river_grid%dwnstrm_index(j)               !index of the most downstream point.  ie end of the line
           k = k + 1
        end do
         
        if (river_grid%land_mask(j) .eq. 2) then  !ended at an ocean point
          river_grid%ocean_outlet(i) = j
          river_grid%upstrm_number(j) = river_grid%upstrm_number(j) + 1
        elseif (river_grid%land_mask(j) .eq. 1) then          !ended at a land point
          river_grid%ocean_outlet(i) = j
          river_grid%upstrm_number(j) = river_grid%upstrm_number(j) + 1
        else
          stop
        end if
         
      else
        river_grid%ocean_outlet(i) = j                                 !if it is ocean it is its own outlet
        river_grid%upstrm_number(j) = river_grid%upstrm_number(j) + 1  !and has only itself as an upstream cell
      end if
    end do
         
    !count the total number of basins.
    river_grid%nrr_cells = 0
    river_grid%nbasins   = 0
    do i=1,river_grid%npts
      if (river_grid%upstrm_number(i) .gt. 0) then   !upstrm > 1 only for
         river_grid%nbasins = river_grid%nbasins + 1
         river_grid%nrr_cells = river_grid%nrr_cells + 1
      end if
    end do
    
  end subroutine associate_ocean_outlet_points    
  
!----------------------------------------------------------------------------!

  subroutine reorder_grid_by_basin(grid_var,basins)
    implicit none
    class(river_grid_type),                   intent(inout) :: grid_var    
    class(basin_type), pointer, dimension(:), intent(inout) :: basins
    !loop: basins
    !loop : find all for basin i
    !  create vector of indices for basin i --> output: basin_indices(nbasins,npts).
    !             basins won't have contiguous memory...hmmmm

    integer :: cnt, i,ii,j,jj,k, kk, total_nbasins, total_land_cells, ncells
    integer, allocatable, dimension(:) :: tmp_indices
    
    type(river_grid_type) :: ord_grid_var   !grid variable ordered so basins are continuous
    
    
    call ord_grid_var%create(grid_var%npts)  !this is the total number of possible river points
    
    ord_grid_var%npts = grid_var%npts
    ord_grid_var%nlon = grid_var%nlon
    ord_grid_var%nlat = grid_var%nlat        
    ord_grid_var%nrr_cells = grid_var%nrr_cells
    ord_grid_var%nbasins   = grid_var%nbasins    
    
    allocate(tmp_indices(grid_var%npts))

    total_nbasins = grid_var%nbasins
    !find the total number of basins with multiple upstream river cells
    j = 0
    do i=1,total_nbasins
      if (grid_var%upstrm_number(i) .gt. 2) then
        j = j + 1
      end if
    end do
    
    grid_var%nbasins = j
    allocate(basins(grid_var%nbasins))    
    
    total_land_cells = 0
    j = 0
    do i=1,total_nbasins
       tmp_indices(:) = 0
       cnt = 0
       !is this a basin with > 2 river cells?  i.e. not just an ocean cell?
       if (grid_var%upstrm_number(i) .gt. 2) then
         j=j+1
       
         do kk=1,grid_var%npts
           if (grid_var%ocean_outlet(kk) .eq. i)  then     !check for land point here?
             cnt = cnt + 1
             tmp_indices(cnt) = kk
           end if
         end do
        
         total_land_cells = total_land_cells + cnt
         ncells = cnt
         basins(j)%n_basin_cells = cnt
         call basins(j)%create(cnt)                      !or can use alloc here as below
         !allocate(basins(i)%river_points(cnt))
         basins(j)%river_points(:) = tmp_indices(1:cnt)  !i like this solution.  simply pass basin indices to loop over. 
                                                        !will need to put these in contiguous array to pass back to master.                                          
       end if
    end do
    
    !the ocean cells have been eliminated through the above process.  npts should also change
    
        ! I need to reorder the basin cells so they are contiguous.
        !then each loop can go from i=pe_start,pe_end; wat_mass(i) = dt*(Fin-Fout) + wat_mass(i)
    cnt=1
    do i=1,grid_var%nbasins
      basins(i)%begind = cnt
      basins(i)%endind = cnt + basins(i)%n_basin_cells - 1
      do kk=1,basins(i)%n_basin_cells
        k = basins(i)%river_points(kk)
        
        ord_grid_var%dwnstrm_index(cnt) = grid_var%dwnstrm_index(k)   !set ordered by basin to the non ordered values
        ord_grid_var%ocean_outlet(cnt)  = grid_var%ocean_outlet(k)
        ord_grid_var%upstrm_number(cnt) = grid_var%upstrm_number(k)
        ord_grid_var%maps%ind_lgr(cnt,:) = grid_var%maps%ind_lgr(k,:)
        ord_grid_var%maps%weight_lgr(cnt,:) = grid_var%maps%weight_lgr(k,:)
        ord_grid_var%maps%n_ovrlap_lgr(cnt) = grid_var%maps%n_ovrlap_lgr(k)
        ord_grid_var%lat(cnt)             = grid_var%lat(k)
        ord_grid_var%lon(cnt)             = grid_var%lon(k)
        ord_grid_var%slope(cnt)           = grid_var%slope(k)
        ord_grid_var%length(cnt)          = grid_var%length(k)
        ord_grid_var%land_mask(cnt)       = grid_var%land_mask(k)
        ord_grid_var%active_cell(cnt)     = grid_var%active_cell(k)
        
        !flow variables yet to be defined.  no need to remap
        ! wat_mass, wat_vol, wat_hgt, wat_length
        cnt = cnt + 1
      end do
    end do 
                                                       
    deallocate(tmp_indices)
    
    !destroy grid var.  make it noew with fewer points (doesn't include the ocean now)
    call grid_var%destroy()
    call grid_var%create(total_land_cells)
    
    grid_var%npts = total_land_cells
    grid_var%nlat = ord_grid_var%nlat
    grid_var%nlon = ord_grid_var%nlon
    
    !then set to grid_var which is now continuous in terms of basins.  map same as before reordered
    !grid_var = ord_grid_var  !have not orderloaded = operator.  this won't work
    grid_var%dwnstrm_index(:) = ord_grid_var%dwnstrm_index(1:total_land_cells)   !set ordered by basin to the non ordered values
    grid_var%ocean_outlet(:)  = ord_grid_var%ocean_outlet(1:total_land_cells)
    grid_var%upstrm_number(:) = ord_grid_var%upstrm_number(1:total_land_cells)
    grid_var%maps%ind_lgr(:,:) = ord_grid_var%maps%ind_lgr(1:total_land_cells,:)
    grid_var%maps%weight_lgr(:,:) = ord_grid_var%maps%weight_lgr(1:total_land_cells,:)
    grid_var%maps%n_ovrlap_lgr(:) = ord_grid_var%maps%n_ovrlap_lgr(1:total_land_cells)
    grid_var%lat(:)             = ord_grid_var%lat(1:total_land_cells)
    grid_var%lon(:)             = ord_grid_var%lon(1:total_land_cells)
    grid_var%slope(:)           = ord_grid_var%slope(1:total_land_cells)
    grid_var%length(:)          = ord_grid_var%length(1:total_land_cells)
    grid_var%land_mask(:)       = ord_grid_var%land_mask(1:total_land_cells)    
    grid_var%active_cell(:)     = ord_grid_var%active_cell(1:total_land_cells)
    
    
    call ord_grid_var%destroy()  !clean up

  end subroutine reorder_grid_by_basin
  
!----------------------------------------------------------------------------!

  subroutine remove_inactive_land_basins(grid_var,basins)
    implicit none
    
    class(river_grid_type),                   intent(inout) :: grid_var
    class(basin_type), pointer, dimension(:), intent(inout) :: basins  
    
    type(basin_type), allocatable, dimension(:) :: cmp_basins
    type(river_grid_type)                       :: cmp_grid_var
    
    integer :: i,j,k,ii,jj,kk,cnt
    integer :: n_active_cells
    integer :: total_active_cells
    
    integer, allocatable, dimension(:) :: active_basin
    
    !find the total number of active cells.  use this to create new grid and basins
    total_active_cells = sum(grid_var%active_cell(:))
    
    
    allocate(active_basin(grid_var%nbasins))
    active_basin(:) = 0    
    
    call cmp_grid_var%create(total_active_cells)
    cmp_grid_var%nrr_cells = total_active_cells
    cmp_grid_var%npts      = total_active_cells
    cmp_grid_var%nlat      = grid_var%nlat
    cmp_grid_var%nlon      = grid_var%nlon
    
    !  determine if each basin is active.  only active if all cells covered by land model.  make cut off 90%???
    do i=1,grid_var%nbasins
      
      n_active_cells = sum(grid_var%active_cell(basins(i)%begind:basins(i)%endind))  !count the number of active cells in the basin
      
      if (n_active_cells .ge. int(0.75*basins(i)%n_basin_cells)) then   !compute basin if we have forcing for > 3/4 of it
        active_basin(i) = 1
      else
        active_basin(i) = 0
      end if
    end do
    
    cmp_grid_var%nbasins = sum(active_basin(:))   
    allocate(cmp_basins(cmp_grid_var%nbasins))
    
    if (cmp_grid_var%nbasins .lt. grid_var%nbasins) then   !there are some basins that aren't active
    
      cnt=1
      do i=1,grid_var%nbasins
    
        if (active_basin(i) .eq. 1) then
      
          cmp_basins(i)%begind = cnt
          cmp_basins(i)%endind = cnt + basins(i)%n_basin_cells - 1
          cmp_basins(i)%n_basin_cells = basins(i)%n_basin_cells
          !use temporaries to make code shorter
          j  = cmp_basins(i)%begind 
          jj = cmp_basins(i)%endind
        
          k  = basins(i)%begind
          kk = basins(i)%endind
      
          !overloading = operator would make this much cleaner
          cmp_grid_var%dwnstrm_index(j:jj) = grid_var%dwnstrm_index(k:kk)
          cmp_grid_var%ocean_outlet(j:jj)  = grid_var%ocean_outlet(k:kk)
          cmp_grid_var%upstrm_number(j:jj)  = grid_var%upstrm_number(k:kk)
          cmp_grid_var%maps%ind_lgr(j:jj,:) = grid_var%maps%ind_lgr(k:kk,:)
          cmp_grid_var%maps%weight_lgr(j:jj,:) = grid_var%maps%weight_lgr(k:kk,:)
          cmp_grid_var%maps%n_ovrlap_lgr(j:jj) = grid_var%maps%n_ovrlap_lgr(k:kk)

          cmp_grid_var%lat(j:jj)             = grid_var%lat(k:kk)
          cmp_grid_var%lon(j:jj)             = grid_var%lon(k:kk)
          cmp_grid_var%slope(j:jj)           = grid_var%slope(k:kk)
          cmp_grid_var%length(j:jj)          = grid_var%length(k:kk)
          cmp_grid_var%land_mask(j:jj)       = grid_var%land_mask(k:kk)      
        
          cnt = cmp_basins(i)%endind + 1
        
        end if
      
      end do
    
      !remove original grid_var variable.  reallocate new one with only active routing cells
      call grid_var%destroy()
    
      call grid_var%create(total_active_cells) 
    
      !copy over all data from the temporary compact river grid (cmp_river_grid)
      grid_var%npts = total_active_cells
      grid_var%nlat = cmp_grid_var%nlat  !nlat for global grid.  not all used
      grid_var%nlon = cmp_grid_var%nlon
      grid_var%nbasins = cmp_grid_var%nbasins
      grid_var%nrr_cells = cmp_grid_var%nrr_cells
    
      grid_var%active_cell(:) = 1          !all river cells are now active
      grid_var%lat(:) = cmp_grid_var%lat(:)
      grid_var%lon(:) = cmp_grid_var%lon(:)
      grid_var%slope(:) = cmp_grid_var%slope(:)
      grid_var%length(:) = cmp_grid_var%length(:)
      grid_var%land_mask(:) = cmp_grid_var%land_mask(:)

      grid_var%dwnstrm_index(:) = cmp_grid_var%dwnstrm_index(:)
      grid_var%ocean_outlet(:)  = cmp_grid_var%ocean_outlet(:)
      grid_var%upstrm_number(:)  = cmp_grid_var%upstrm_number(:)
      grid_var%maps%ind_lgr(:,:) = cmp_grid_var%maps%ind_lgr(:,:)
      grid_var%maps%weight_lgr(:,:) = cmp_grid_var%maps%weight_lgr(:,:)
      grid_var%maps%n_ovrlap_lgr(:) = cmp_grid_var%maps%n_ovrlap_lgr(:)
    
      !do the same for the basin variable
      do i=1,size(basins(:))
        call basins(i)%destroy()
      end do
      deallocate(basins)
    
      allocate(basins(grid_var%nbasins))

      do i=1,grid_var%nbasins
        basins(i)%n_basin_cells = cmp_basins(i)%n_basin_cells
        call basins(i)%create(basins(i)%n_basin_cells)
        basins(i)%begind = cmp_basins(i)%begind
        basins(i)%endind = cmp_basins(i)%endind
        do k=1,basins(i)%n_basin_cells
          basins(i)%river_points(k) = k + cmp_basins(i)%begind - 1
        end do
      
      end do
      
    end if  !some basins are not active.
    
    
    do i=1,size(cmp_basins(:))
      call cmp_basins(i)%destroy()
    end do
    
    deallocate(cmp_basins)
    deallocate(active_basin)
    
    

    
  end subroutine remove_inactive_land_basins
    
!----------------------------------------------------------------------------! 
  !Determine the wieghts and mapping from 1D global land grid to the 1D river grid
  !river_grid%{lat,lon} are assumed higher resolution than lat_out lon_out?
  
  subroutine determine_hilo_res_mapping(river_grid,lat_lo,lon_lo)
  
    implicit none
    
    class(river_grid_type), intent(inout) :: river_grid
    real(r_2), dimension(:),     intent(in)  :: lat_lo,lon_lo  !the lo resolution grid

    !local variables
    integer   :: i,j,ii,jj,k,kk,k_tmp,kk_tmp  !integer counters for the loops

    real(r_2) :: dlat_lo,dlon_lo   !grid cell size (degrees) of lo res grid
    real(r_2) :: dlat_hi,dlon_hi !grid cell size (degrees) of hi res grid

    real(r_2) :: dlone,dlonw,dx,dlats,dlatn,dy   !cell size (degrees) and m (dx) in hi-lo  cells that overlap
    
    real(r_2) :: dy_lo,dy_hi,area_lo,area_hi  !size of the north-south side of the lo resolution grid (lo), area of the lo res grid
    
    real(r_2) :: Eedge_hi, Eedge_lo   !eastern  edge (lon) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Wedge_hi, Wedge_lo   !western  edge (lon) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Sedge_hi, Sedge_lo   !southern edge (lat) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Nedge_hi, Nedge_lo   !northern edge (lat) of high resolution (hi) and lo res (lo) grid cells

    !are lat and lon global 1D vectors or only vectors with lat lon values?
    if (size(river_grid%lat) .gt. river_grid%nlat .and. size(river_grid%lon) .gt. river_grid%nlon) then
       dlat_hi  = river_grid%lat(river_grid%nlon+1)-river_grid%lat(1)
       dlon_hi  = river_grid%lon(2) - river_grid%lon(1)
    else
       dlat_hi  = river_grid%lat(2) - river_grid%lat(1)
       dlon_hi  = river_grid%lon(2) - river_grid%lon(1)
    end if

    if (size(lat_lo) .gt. mlat .and. size(lon_lo) .gt. mlon) then
       dlat_lo  = lat_all(1,2)-lat_all(1,1)  !array is compacted,  must use lat_all lon_all brought in from io_vars module
       dlon_lo  = lon_all(2,1) - lon_all(1,1)
    else
       dlat_lo  = lat_lo(2) - lat_lo(1)
       dlon_lo  = lon_lo(2) - lon_lo(1)
    end if

    !init overlap mappings and areas
    river_grid%maps%n_ovrlap_lgr(:) = 0
    river_grid%maps%weight_lgr(:,:) = 0._r_2
    river_grid%maps%ind_lgr(:,:)    = -1
    
    !initialize all cells to inactive
    river_grid%active_cell(:)                = 0

    !map each input lat/lon to output grid and calc the fraction of cell within (1 cell can map to upto four output)
    do k=1,mp   !loop over all the land grid points
    
      Sedge_lo = lat_lo(k)  - dlat_lo/2._r_2
      Nedge_lo = lat_lo(k)  + dlat_lo/2._r_2    
      
      Wedge_lo = lon_lo(k)  - dlon_lo/2._r_2
      Eedge_lo = lon_lo(k)  + dlon_lo/2._r_2
    
      dy_lo = sin(deg2rad*Nedge_lo) - sin(deg2rad*Sedge_lo)
      area_lo = dy_lo * dlon_lo * re * re  !adjust for fraction of grid cells??
      
      do kk=1,river_grid%npts

        Sedge_hi = river_grid%lat(kk) - dlat_hi/2._r_2  !avoid some calcs dooing it here
        Nedge_hi = river_grid%lat(kk) + dlat_hi/2._r_2
        
        dy_hi = sin(deg2rad*Nedge_hi) - sin(deg2rad*Sedge_hi)
        area_hi = dy_hi * dlon_hi * re * re  !adjust for fraction of grid cells??

        Wedge_hi = river_grid%lon(kk) - dlon_hi/2._r_2
        Eedge_hi = river_grid%lon(kk) + dlon_hi/2._r_2

        if ((Wedge_hi .le. Eedge_lo) .and. (Eedge_hi .ge. Wedge_lo) .and. &
            (Sedge_hi .le. Nedge_lo) .and. (Nedge_hi .ge. Sedge_lo)) then
              
          river_grid%maps%n_ovrlap_lgr(kk) = river_grid%maps%n_ovrlap_lgr(kk) + 1  !number of overlapping cells
          
          dlone = min(Eedge_lo,Eedge_hi)*deg2rad !determine area of input cell within the output grid cell
          dlonw = max(Wedge_lo,Wedge_hi)*deg2rad 
          dx = max(0.0,(dlone-dlonw))

          dlatn = min(Nedge_lo,Nedge_hi)*deg2rad 
          dlats = max(Sedge_lo,Sedge_hi)*deg2rad 
          dy = max(0.0,(sin(dlatn)-sin(dlats)))

          river_grid%maps%ind_lgr(kk,river_grid%maps%n_ovrlap_lgr(kk))    = k                   !lsm point for given river point
          river_grid%maps%weight_lgr(kk,river_grid%maps%n_ovrlap_lgr(kk)) = dx*dy / area_lo     !fraction of lsm cell k occupied by river cell
              
          river_grid%active_cell(kk) = 1             !mark this river cell as active.

        end if  !test lon overlap
        
      end do !kk loop over input lat
      
    end do  !loop over mp land points
    
  end subroutine determine_hilo_res_mapping 
  
!----------------------------------------------------------------------------!

   subroutine step_river_routing(river,river_grid,basins,basins_pe_start,basins_pe_end)
     implicit none
     
     class(river_flow_type),          intent(inout) :: river   !contains mass,flow variables
     class(river_grid_type),          intent(in)    :: river_grid
     class(basin_type), dimension(:), intent(in)    :: basins          !contains info on each basin    
     integer,                         intent(in)    :: basins_pe_start,basins_pe_end
     
!      real(r_2), pointer, dimension(:) :: river_water
!      real(r_2), pointer, dimension(:) :: river_vel
!      real(r_2), pointer, dimension(:) :: river_hgt
!      real(r_2), pointer, dimension(:) :: river_distance
!      real(r_2), pointer, dimension(:) :: Fout
!      real(r_2), pointer, dimension(:) :: Fin     
     
!      integer :: npoints_local_basin
!      integer :: i,j,k,ii,jj,kk
     
!      allocate(Fin(river_grid%npts))
!      allocate(Fout(river_grid%npts))
     
!      river_water => river%mass(:)                            !pointers to full arrays.  each pe has full copy of array?
!      river_vel   => river%vel(:)
!      river_hgt   => river%hgt(:)
!      river_distance => river_grid%length(:)
!      Fin         => river%Fin(:)
!      Fout        => river%Fout(:)
     
!      do i=basins_pe_start,basins_pe_end                       !loop over basins assinged to this pe
!        npoints_local_basin = basins(i)%n_basin_cells
!        do k=1,npoints_local_basin
!          kk = basins(i)%river_points(k)
!          Fin(kk) = 0._r_2
!        end do
!
!      do i=basins_pe_start,basins_pe_end
!        j  = basins(i)%begind
!        jj = basins(i)%endind
!        river%Fin(j:jj) = 0._r_2
!        do k=j,jj
!          river%Fout(k) = river%mass(k)*river%vel(k)/river_grid%length(k)
!          river%Fin(river_grid%dwnstrm_index(k)) = river%Fin(river_grid%dwnstrm_index(k)) + Fout(k)
!        end do
!
!        river%mass(j:jj) = river%mass(j:jj) - (river%Fout(j:jj) + river%Fin(j:jj))*dels
!
!        call basins(i)%get_outflow(river)  !write a subroutine to compute total basin outflow
!
!
!      end do  !loop over this pe basins
       
!        do k=1,npoints_local_basin
!          kk = basins(i)%river_points(k)
!          Fout(kk) = river_water(kk)*river_vel(kk)/river_disctance(kk)/river_distance(kk)*dels  !made up
!          Fin(grid_var%dwnstrm_index(kk)) = Fin(grid_var%dwnstrm_index(kk)) + Fout(kk)
!        end do
       
!        do k=1,npoints_local_basin
!          kk = basins(i)%river_points(k)
!          river_water(kk) = river_water(kk) + (Fin(kk)-Fout(kk))*dels
!          if (grid_var%is_outlet .eq. 1) then
!            total_outflow(i) = Fout(kk)*dels
!          end if
!        end do
     
!      end do
     
    end subroutine step_river_routing
    
!----------------------------------------------------------------------------!

!**********Below are the Routines For Allocating/Deallocating****************!

  subroutine alloc_river_flow(var,npts)
    implicit none
    class(river_flow_type),intent(inout)  :: var
    integer, intent(in) :: npts

    allocate(var%mass_init(npts))
    var%vol(:) = fNaN
    
    allocate(var%mass(npts))
    var%mass(:) = fNaN
    
    allocate(var%vel(npts))
    var%vel(:) = fNaN
    
    allocate(var%hgt(npts))
    var%hgt(:) = fNaN    
    
    allocate(var%runoff_flux(npts))
    var%runoff_flux(:) = fNaN
    
    allocate(var%Fin(npts))
    var%Fin(:) = fNaN      
    
    allocate(var%Fout(npts))
    var%Fout(:) = fNaN       

  end subroutine alloc_river_flow
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_flow(var)
    implicit none
    class(river_flow_type),intent(inout) :: var

    deallocate(var%mass_init)
    deallocate(var%mass)
    deallocate(var%vel)
    deallocate(var%hgt)
    deallocate(var%runoff_flux)
    deallocate(var%Fin)
    deallocate(var%Fout)

  end subroutine dealloc_river_flow
  
!----------------------------------------------------------------------------!

  subroutine alloc_river_map_grid(var,npts)
    implicit none
    class(map_grid_type),intent(inout)  :: var
    integer, intent(in)                :: npts

    allocate(var%ind_lgr(npts,max_n_ovrlap))  !outut land cell number form a given river cell number
    var%ind_lgr(:,:) = iNaN
    
    allocate(var%weight_lgr(npts,max_n_ovrlap))
    var%weight_lgr(:,:) = fNaN
    
    allocate(var%n_ovrlap_lgr(npts))
    var%n_ovrlap_lgr(:) = iNaN    
    

  end subroutine alloc_river_map_grid
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_map_grid(var)
    implicit none
    class(map_grid_type), intent(inout) :: var
    
    deallocate(var%ind_lgr)    
    deallocate(var%weight_lgr)    
    deallocate(var%n_ovrlap_lgr)    
    
  end subroutine dealloc_river_map_grid
  
!----------------------------------------------------------------------------!

  subroutine alloc_river_grid(var,npts)
    implicit none
    class(river_grid_type), intent(inout) :: var
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

    call var%maps%create(npts)         


  end subroutine alloc_river_grid
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_grid(var)
    implicit none
    class(river_grid_type), intent(inout) :: var

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
    
    call var%maps%destroy()

  end subroutine dealloc_river_grid

!----------------------------------------------------------------------------!

  subroutine alloc_basin(var,npts)
    implicit none
    class(basin_type), intent(inout) :: var
    integer,           intent(in)    :: npts
    
    allocate(var%river_points(npts))
    var%river_points(:) = iNaN
    
  end subroutine alloc_basin
  
!----------------------------------------------------------------------------!

  subroutine dealloc_basin(var)
    implicit none
    class(basin_type), intent(inout) :: var
    
    deallocate(var%river_points)
    
  end subroutine dealloc_basin  

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
    
  subroutine read_nc_int(ncid,var_name,start_inds,end_inds,var_data,mask)
    use netcdf
    implicit none
    
    integer,               intent(in)    :: ncid
    character(len=*)   , intent(in)    :: var_name
    integer, dimension(:), intent(in)    :: start_inds,end_inds
    integer, dimension(:), intent(inout) :: var_data
    integer, dimension(:,:), optional    :: mask
    
    integer :: nc_check,nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    integer, dimension(:,:), allocatable :: int_data
    
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1
    
    allocate(int_data(nlon_rr,nlat_rr))
    
    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,int_data(:,:),start_inds,end_inds) )
    
    if (present(mask)) then
    
      k=1
      do j=1,nlat_rr
        do i=1,nlon_rr
          if (mask(i,j) .ge. 1) then
            var_data(k) = int_data(i,j)
            k=k+1
          end if
        end do
      end do
      
    else   !mask data not present.  get all data
    
      k=1
      do j=1,nlat_rr
        do i=1,nlon_rr
          var_data(k) = int_data(i,j)
          k=k+1
        end do
      end do
    
    end if
     
     deallocate(int_data)
     
   end subroutine read_nc_int

!----------------------------------------------------------------------------!   
 
  subroutine read_nc_flt(ncid,var_name,start_inds,end_inds,var_data,mask)
    use netcdf
    implicit none
    
    integer,                 intent(in)    :: ncid
    character(len=*)       , intent(in)    :: var_name
    integer  , dimension(:), intent(in)    :: start_inds,end_inds
    real(r_2), dimension(:), intent(inout) :: var_data
    integer, dimension(:,:), optional    :: mask    
    
    integer :: nc_check,nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    real(r_2), dimension(:,:), allocatable :: flt_data
    
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1
    
    allocate(flt_data(nlon_rr,nlat_rr))
    
    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,flt_data(:,:),start_inds,end_inds) )

    if (present(mask)) then
    
      k=1
      do j=1,nlat_rr
        do i=1,nlon_rr
          if (mask(i,j) .ge. 1) then
            var_data(k) = flt_data(i,j)
            k=k+1
          end if
        end do
      end do
      
    else   !mask data not present.  get all data
    
      k=1
      do j=1,nlat_rr
        do i=1,nlon_rr
          var_data(k) = flt_data(i,j)
          k=k+1
        end do
      end do
    
    end if

    deallocate(flt_data)
     
  end subroutine read_nc_flt 



end module cable_routing
