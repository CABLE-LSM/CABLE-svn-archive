!compile using preproc for using mpi.
! compile with -Dname=def which is the same as
!  #define name def
!
!  this module needs to have -DRRM=1  (RRM => River Routing Module)
!
!  so if using mpi (River Routing Modue _ MPI )
!  -DRRM_MPI=1
!  not using mpi don't include this option
!
!
!  -Dwith_cable  link to the cable modules
!  else define local parameters that are needed from cable

module cable_routing

#ifdef RRM

  use netcdf
  use cable_rrm_nc_names  
  
#ifdef with_cable
  use cable_types
  use cable_common_module, only : cable_user
  use cable_def_types_mod, only : r_2, ms, mp,mlat,mlon

  use cable_IO_vars_module, only : landpt, patch, max_vegpatches, metGrid,  &
                                  land_x,land_y,latitude,longitude,         &
                                  lat_all, lon_all,                         &
                                  logn,output,xdimsize,ydimsize,check,mask
#endif

  implicit none

#ifdef RRM_MPI
  include 'mpif.h'
#endif


  public
  !**************************************************************************!
  !  Temporary to avoid having to link to the CABLE mods while testing       !
  !**************************************************************************!
#ifndef with_cable 
  integer, parameter :: r_2  = SELECTED_REAL_KIND(8)
  integer, parameter :: mlat = 35!
  integer, parameter :: mlon = 50!360
  integer, parameter :: mp = 50*35!360*150
  real(r_2) :: dels = 1800.0
  
  real(r_2), dimension(:), allocatable, save :: latitude
  real(r_2), dimension(:), allocatable, save :: longitude

  real(r_2), dimension(:,:), allocatable, save :: lat_all
  real(r_2), dimension(:,:), allocatable, save :: lon_all   
  real(r_2), dimension(:), allocatable, save :: pft_frac_lo

  !MPI
  integer :: comm,ierr
#endif  
  !**************************************************************************!
  !  Temporary to avoid having to link to the CABLE mods while testing       !
  !**************************************************************************!  
  
  !parameters
  integer,   parameter :: n_river_tracers = 1   !number of tracers...liq  -->carbon,nitrogen and ice??.  not in use as of now
  integer,   parameter :: max_n_ovrlap = 16!324   !assume at most 1x1 to 0.0625x0.0625 -->16x16  do 18x18=324 to be safe
  real(r_2), parameter :: re = 6371000.0              !radius of the earth (km)
  real(r_2), parameter :: deg2rad = 3.14159/180.0     !constant converts degrees to radians
  real(r_2), parameter :: eps=1e-5                    !tolerance parameter.  UNUSED NOW 

  real(r_2), parameter :: river_theta = 0.5
  real(r_2), parameter :: subsurf_theta = 0.4
  real(r_2), parameter :: source_area_channel_cutoff = 50000.0
  real(r_2), parameter :: del_river = 600.0

  real(r_2), parameter :: fNaN = -1e36
  integer, parameter   :: iNaN = -999999
  integer, parameter   :: min_pts_per_basin = 25  !min number of cells to call group a basin
  
  
  type map_grid_type
    
    integer  , allocatable, dimension(:,:)   :: ind_lgr        !index of the Land cell  Given the River cell
    real(r_2), allocatable, dimension(:,:)   :: weight_lgr     !fraction of the river cell covered by land cell relative to total river cell area
                                                         !Land Going to River
    real(r_2), allocatable, dimension(:,:)   :: weight_rgl     !fraction of the river cell covered by land cell relative to total land cell area
                                                         !River Going to Land
    integer  , allocatable, dimension(:)     :: n_ovrlap_lgr   !number of land cells that over lap the river cell  

  end type map_grid_type

  
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
    !real(r_2), allocatable, dimension(:) :: topo_ind        !topographic index -- ln(area/tan(slope))
    !integer,   allocatable, dimension(:) :: basin_ind       !basin index to group basins -- number large to small
    integer,   allocatable, dimension(:) :: land_mask       !1=land,0=ocean,2=land bordering ocean
    integer,   allocatable, dimension(:) :: dwnstrm_index   !index of cell flow goes towards
    integer,   allocatable, dimension(:) :: ocean_outlet    !index of the ending grid cell
    integer,   allocatable, dimension(:) :: upstrm_number   !number of grid cells upstream from the current
    integer,   allocatable, dimension(:) :: active_cell     !1=cell is active, 0=cell not active
    integer,   allocatable, dimension(:) :: direction       !direction of the river flow
    integer,   allocatable, dimension(:) :: is_main_channel !main river channel = 1, not main channel = 0
    integer,   allocatable, dimension(:) :: orig_ind        !the index of the original global array that has all lat/lon points

    type(map_grid_type) :: maps
       
    
  end type river_grid_type
  
  
  type river_flow_type     !one global.  make 2D and have more than one tracer??????
  
    real(r_2), allocatable, dimension(:) :: mass_init   !river water mass at start of timestep (total. not per m2 like the lsm)
    real(r_2), allocatable, dimension(:) :: mass        !river water mass  (total.  not per m2 like the lsm)
    real(r_2), allocatable, dimension(:) :: mass_sub    !river water mass  (total.  not per m2 like the lsm)
    real(r_2), allocatable, dimension(:) :: vel         !river water speed
    real(r_2), allocatable, dimension(:) :: hgt         !river water height 
    real(r_2), allocatable, dimension(:) :: srf_runoff_flux !surface runoff flux (+ to river) from lsm
    real(r_2), allocatable, dimension(:) :: sub_runoff_flux !subsurface runoff flux (+ to river) from lsm    
    real(r_2), allocatable, dimension(:) :: Fin
    real(r_2), allocatable, dimension(:) :: Fin_sub
    real(r_2), allocatable, dimension(:) :: Fout
    
  end type river_flow_type
  
  type basin_type
  
    integer, allocatable, dimension(:)   :: begind            !basin index start in reordered global 1d array
    integer, allocatable, dimension(:)   :: endind            !basin end index in global reordered array
    integer, allocatable, dimension(:)   :: n_basin_cells     !number of river cells in the basin
    integer, allocatable, dimension(:,:) :: river_points  !the indices for the basin in the unordered global 1d array
    
  end type basin_type

  !below will go into cable_routing_main_routine.  !global on myrank =0, local on myrank=1->nprocs
  type(river_grid_type), pointer, save :: global_river_grid      , local_river_grid
  type(river_flow_type), pointer, save :: global_river           , local_river
  type(basin_type)     , pointer, save :: global_basins, local_basins
  
  
  !outline
  !if timestep==0 --> 
  ! if myrank = 0:
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
  !                                         !currently ocean cells are single cell basins
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
  !    Send patition info to all workers
  !
  !end if myrank = 0
  ! if myrank /= 0 
  !   get data from 0
  !   allocate local river,river_grid,basin variables on each non-0 myrank
  !
  !start time stepping
  !
  !  call river%get_runoff(Qrunoff_lsm,river_grid)
  !  if myrank is 0.  after mpimaster has gathered output from the workers....
  !call map_qrunoff_lsm_to_river(Qsrf,Qsubsrf,river,river_grid)  !treat surface and subsurface the same
  !
  !  !send from myrank=0 to the rest.
  !
  !call step_river(river,river_grid)                             !!!!split into main river channel segments and overland/subsurface routing flow
  !             => river%{wat_mass,wat_vol,wat_hgt,wat_vel}
  !             => river_grid%{distance,mpi%bg,mpi%ed}
  !   call river_mass_balance()
  !
  !  send back to myrank = 0
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

  public map_lsm_runoff_to_river,get_river_route_data,find_downstream_index,associate_ocean_outlet_points
  public reorder_grid_by_basin,remove_inactive_land_basins,find_main_river_channels,check_nc


  interface create
    module procedure alloc_river_flow ! interface for a module
    module procedure alloc_river_grid ! procedure is implicit
    module procedure alloc_basin
    module procedure alloc_river_map_grid
  end interface

  interface destroy
    module procedure dealloc_river_flow ! interface for a module
    module procedure dealloc_river_grid ! procedure is implicit
    module procedure dealloc_basin
    module procedure dealloc_river_map_grid
  end interface
  
  interface read_nc
    module procedure read_nc_flt
    module procedure read_nc_int
  end interface
 

contains

!----------------------------------------------------------------------------!

  subroutine map_lsm_runoff_to_river(river_var,grid_var,Qsrf_runoff_lsm,Qsub_runoff_lsm)
    implicit none
    type(river_flow_type), target,intent(inout) :: river_var
    real(r_2), dimension(:), intent(in)   :: Qsrf_runoff_lsm, Qsub_runoff_lsm
    type(river_grid_type), target,intent(in)    :: grid_var

    !local pointers to the derived variables
    real(r_2), pointer, dimension(:)   :: srf_runoff_flux
    real(r_2), pointer, dimension(:)   :: sub_runoff_flux
    real(r_2), pointer, dimension(:)   :: area
    real(r_2), pointer, dimension(:,:) :: weight
    integer,   pointer, dimension(:,:) :: ind
    integer,   pointer, dimension(:)   :: n_ovrlap

    integer :: kk,i,ii

    !local pointers to the derived variables.  shorten names 
    srf_runoff_flux => river_var%srf_runoff_flux(:)
    sub_runoff_flux => river_var%sub_runoff_flux(:)
    weight          => grid_var%maps%weight_lgr(:,:)
    area            => grid_var%area(:)
    ind             => grid_var%maps%ind_lgr(:,:)
    n_ovrlap        => grid_var%maps%n_ovrlap_lgr(:)

    !determine runoff input in river cells from lsm:
    srf_runoff_flux(:) = 0._r_2
    sub_runoff_flux(:) = 0._r_2    
    do kk=1,grid_var%npts

      do i=1,n_ovrlap(kk)

        ii = ind(kk,i)
        srf_runoff_flux(kk) = srf_runoff_flux(kk) +  Qsrf_runoff_lsm(ii)*weight(kk,i)*area(kk)  !convert from mm/m2 => mm

        sub_runoff_flux(kk) = sub_runoff_flux(kk) + Qsub_runoff_lsm(ii)*weight(kk,i)*area(kk)  !convert from mm/m2 => mm
      end do

    end do
    
    nullify(srf_runoff_flux)
    nullify(sub_runoff_flux)
    nullify(weight)
    nullify(ind)
    nullify(n_ovrlap)

  end subroutine map_lsm_runoff_to_river

!----------------------------------------------------------------------------!

  subroutine map_river_mass_to_lsm(river_var,grid_var,river_mass_lsm,subsrf_mass_lsm)
    implicit none
    type(river_flow_type), target,intent(in)      :: river_var
    type(river_grid_type), target,intent(in)      :: grid_var
    real(r_2),  intent(inout)               :: river_mass_lsm(:)
    real(r_2), intent(inout)                :: subsrf_mass_lsm(:)
    
    !local pointers to the derived variables
    real(r_2), pointer, dimension(:)   :: area
    real(r_2), pointer, dimension(:,:) :: weight
    integer,   pointer, dimension(:,:) :: ind
    integer,   pointer, dimension(:)   :: n_ovrlap    

    integer :: kk,i,ii
    
    !local pointers to the derived variables.  shorten names 
    weight          => grid_var%maps%weight_rgl(:,:)     
    area            => grid_var%area(:)
    ind             => grid_var%maps%ind_lgr(:,:)
    n_ovrlap        => grid_var%maps%n_ovrlap_lgr(:)    

    !determine runoff input in river cells from lsm:
    river_mass_lsm(:) = 0.
    !subsrf_mass_lsm(:) = 0.

    do kk=1,grid_var%npts

      do i=1,n_ovrlap(kk)

        ii = ind(kk,i) 
        river_mass_lsm(ii)   = river_mass_lsm(ii)  + river_var%mass(kk)    * weight(kk,i)/area(kk)  !mm => mm/m2
        subsrf_mass_lsm(ii)  = subsrf_mass_lsm(ii) + river_var%mass_sub(kk)* weight(kk,i)/area(kk)  !mm => mm/m2
        
      end do

    end do
    
    nullify(weight)
    nullify(area)
    nullify(ind)
    nullify(n_ovrlap)


  end subroutine map_river_mass_to_lsm  


!----------------------------------------------------------------------------!

  subroutine create_river_grid_copy(grid_var,new_grid_var,ncells)
  
    type(river_grid_type), intent(inout),pointer    :: grid_var
    type(river_grid_type), intent(inout),pointer :: new_grid_var    
    integer, optional,      intent(in)    :: ncells
    !local variables
    integer :: ntot

    if (present(ncells)) then !not used anymore
      ntot = ncells
    else
      ntot = grid_var%npts
    end if

    new_grid_var => grid_var

!not done

  end subroutine create_river_grid_copy

!----------------------------------------------------------------------------!

  subroutine copy_river_grid_vector_values(grid_var,new_grid_var,n1,n2,n3,n4)
  
    !usage: new_var(n1:n2) = old_var(n3:n4).  if single value use n1=n2 and n3=n4
    !if n3 and n4 are missing then old_var(1:old_var%npts)
    !if all are missing then it is new_var(1:npts)=old_var(1:npts) where npts = old_var%npts
    implicit none
    type(river_grid_type), intent(inout) :: grid_var
    type(river_grid_type), intent(inout) :: new_grid_var
    
    integer, optional :: n1
    integer, optional :: n2
    integer, optional :: n3
    integer, optional :: n4
    integer :: m1,m2,m3,m4
    
    m4 = new_grid_var%npts
    if (present(n4)) m4 = n4
    m3 = 1
    if (present(n3)) m3 = n3
    m2 = grid_var%npts
    if (present(n2)) m2 = n2
    m1 = 1
    if (present(n1)) m1 = n1

    new_grid_var%dwnstrm_index(m1:m2)   = grid_var%dwnstrm_index(m3:m4)   !copy all the values
    new_grid_var%ocean_outlet(m1:m2)    = grid_var%ocean_outlet(m3:m4)
    new_grid_var%upstrm_number(m1:m2)   = grid_var%upstrm_number(m3:m4)
    new_grid_var%lat(m1:m2)             = grid_var%lat(m3:m4)
    new_grid_var%lon(m1:m2)             = grid_var%lon(m3:m4)
    new_grid_var%slope(m1:m2)           = grid_var%slope(m3:m4)
    new_grid_var%length(m1:m2)          = grid_var%length(m3:m4)
    new_grid_var%land_mask(m1:m2)       = grid_var%land_mask(m3:m4)    
    new_grid_var%active_cell(m1:m2)     = grid_var%active_cell(m3:m4)
    new_grid_var%is_main_channel(m1:m2) = grid_var%is_main_channel(m3:m4) 
    new_grid_var%orig_ind(m1:m2)        = grid_var%orig_ind(m3:m4)
    !don't forget the maps derived type!
    new_grid_var%maps%ind_lgr(m1:m2,:)    = grid_var%maps%ind_lgr(m3:m4,:)
    new_grid_var%maps%weight_lgr(m1:m2,:) = grid_var%maps%weight_lgr(m3:m4,:)
    new_grid_var%maps%weight_rgl(m1:m2,:) = grid_var%maps%weight_rgl(m3:m4,:)
    new_grid_var%maps%n_ovrlap_lgr(m1:m2) = grid_var%maps%n_ovrlap_lgr(m3:m4)


  end subroutine copy_river_grid_vector_values


!----------------------------------------------------------------------------!
  function compute_global_mass_balance_error(river_var,grid_var,basin_var) result(mass_error)
  
    implicit none
    type(river_flow_type),  target,intent(in)  :: river_var
    type(river_grid_type),  target,intent(in)  :: grid_var
    type(basin_type),              intent(in)  :: basin_var
    real(r_2) :: mass_error
    !local variables
    !pointers to derived type arrays
    real(r_2), pointer, dimension(:)  :: mass_srf
    real(r_2), pointer, dimension(:)  :: mass_sub
    real(r_2), pointer, dimension(:)  :: mass_init
    real(r_2), pointer, dimension(:)  :: mass_sub_init
    real(r_2), pointer, dimension(:)  :: srf_runoff_flux
    real(r_2), pointer, dimension(:)  :: sub_runoff_flux
    integer  , pointer, dimension(:)  :: upstrm_number
    real(r_2), pointer, dimension(:)  :: area
    real(r_2), pointer, dimension(:)  :: Fout
    
    real(r_2) :: total_mass, total_lsm_flux,total_outflow, init_total_mass
    integer :: i,j,bg,ed
    
    
    mass_srf        => river_var%mass(:)
    mass_sub        => river_var%mass_sub(:)
    mass_init       => river_var%mass_init(:)
    srf_runoff_flux => river_var%srf_runoff_flux(:)
    sub_runoff_flux => river_var%sub_runoff_flux(:)
    Fout            => river_var%Fout(:)
    
    upstrm_number   => grid_var%upstrm_number(:)
    
    
    init_total_mass = 0._r_2
    total_mass      = 0._r_2
    total_lsm_flux  = 0._r_2
    mass_error      = 0._r_2
    
    total_mass      = sum(mass_srf(:)) + sum(mass_sub(:))
    init_total_mass = sum(mass_init(:))
    total_lsm_flux  = (sum(srf_runoff_flux(:)) + sum(sub_runoff_flux(:))) * dels
    !need to compute outflow
    
    total_outflow = 0._r_2
    
    do i=1,grid_var%nbasins
    
      bg = basin_var%begind(i)
      ed = basin_var%endind(i) 
      
      j = maxloc(upstrm_number(bg:ed),dim=1) + bg - 1 !maxloc returns relative to indices.  cell with most upstream values is the basin outlet
      
      total_outflow = total_outflow + Fout(j)*area(j)*dels
    end do
    
    mass_error = total_mass - init_total_mass + total_lsm_flux - total_outflow
    
    nullify(Fout)
    nullify(sub_runoff_flux)
    nullify(srf_runoff_flux)
    nullify(mass_init)
    nullify(mass_srf)
    nullify(mass_sub)
    
    
  end function compute_global_mass_balance_error
    

!----------------------------------------------------------------------------!

  subroutine get_river_route_data(grid_var,filename)
  
  !reads while file.  need to determine if covered by lsm grid later
    implicit none

    type(river_grid_type), intent(inout),pointer  :: grid_var    
    character(len=250), intent(in)                :: filename
    
    integer :: ncid_river
    integer :: rr_lat_dim_id,rr_lon_dim_id,rr_lat_var_id,rr_lon_var_id
    integer :: nlat_rr_file, nlon_rr_file, npts_rr_file
    
    
    real(r_2), dimension(:)  , allocatable :: lat_data
    real(r_2), dimension(:)  , allocatable :: lon_data
    
    integer :: i,j,k     !integers to loop through lsm data
    
    integer, dimension(2) :: start_inds
    integer, dimension(2) :: end_inds
    integer :: nc_check
    
    call check_nc( nf90_open(trim(filename), nf90_nowrite, ncid_river) )

    !get dim ids
    write(*,*) trim(rr_lat_dim_name)
    write(*,*) trim(rr_lon_dim_name)

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
    
    npts_rr_file = nlon_rr_file * nlat_rr_file
    
    !allocate variable for the river grid
    call alloc_river_grid(grid_var,npts_rr_file)
    
    grid_var%npts = npts_rr_file  !includes all lat/lon even ocean
    grid_var%nlat = nlat_rr_file  
    grid_var%nlon = nlon_rr_file
    !fill lat/lon in grid_var variable
    do j=1,nlat_rr_file
      do i=1,nlon_rr_file
        k = (j-1)*nlon_rr_file + i
        grid_var%lat(k) = lat_data(j)
        grid_var%lon(k) = lon_data(i)
        grid_var%orig_ind(k) = k
      end do
    end do

    
    !read the data
    call read_nc(ncid_river,mask_name    ,start_inds,end_inds,grid_var%land_mask(:))
    call read_nc(ncid_river,rdir_name    ,start_inds,end_inds,grid_var%direction(:))
    call read_nc(ncid_river,length_name  ,start_inds,end_inds,grid_var%length(:))
    call read_nc(ncid_river,slope_name   ,start_inds,end_inds,grid_var%slope(:))
    call read_nc(ncid_river,elev_name    ,start_inds,end_inds,grid_var%elev(:))
    call read_nc(ncid_river,src_area_name,start_inds,end_inds,grid_var%source_area(:))
    !call read_nc(ncid_river,topo_ind_name,start_inds,end_inds,grid_var%topo_ind(:))
    !call read_nc(ncid_river,basin_ind_name,start_inds,end_inds,grid_var%basin_ind(:))
    
    nc_check = nf90_close(ncid_river)
    
        
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
    
    integer :: i,j,k,ii,jj,kk
    integer :: lon_wrap

    !make sure we called this prior to reordering or removing grid points
    if (grid_var%npts .ne. grid_var%nlat*grid_var%nlon) &
         stop "Must find the downstream index using the entire global grid not a subsection"

    if (is_global) then
      lon_wrap = grid_var%nlon
    else
      lon_wrap = 1
    end if
    
    grid_var%dwnstrm_index(:) = 0
    
    do j=1,grid_var%nlat
      do i=1,grid_var%nlon
      
        k = i + (j-1)*grid_var%nlon
        
        if (grid_var%direction(k) /= -9999) then
        
          ii = i + dirc2latindex(grid_var%direction(k))
          jj = j + dirc2lonindex(grid_var%direction(k))
         
         if (ii .lt. 1 )            ii = ii + lon_wrap!grid_var%nlon  
         if (ii .gt. grid_var%nlon) ii = ii - lon_wrap!grid_var%nlon

         if (jj .lt. 1) jj = 1
         if (jj .gt. grid_var%nlat) jj = grid_var%nlat
          

          if (jj .lt. 1 .or. jj .gt. grid_var%nlat .or. ii .lt. 1 .or. ii .gt. grid_var%nlon) then
            write(*,*) ii
            write(*,*) jj
            stop "needs checks yo"
          endif
          
          kk = ii + (jj-1)*grid_var%nlon
          
          grid_var%dwnstrm_index(k) = kk

        endif
      enddo
    enddo
    

  end subroutine find_downstream_index
  
!----------------------------------------------------------------------------! 

  subroutine associate_ocean_outlet_points(grid_var)
  
    implicit none
    
    type(river_grid_type), intent(inout) :: grid_var
  
    !local variables
    integer :: i,j,k,j_prev
    logical :: continue_search
  
    grid_var%ocean_outlet(:)  = -1
    grid_var%upstrm_number(:) = 0   !upstrm_number only non-zero for outlet cells
  
!$OMP PARALLEL DO PRIVATE(i,j,k,j_prev) 
    do i=1,grid_var%npts
      j = i
      if (grid_var%land_mask(i) .eq. 1) then   !it is a land cell
      
        k = 0
        j_prev = j
        do while (k < grid_var%npts .and. grid_var%land_mask(j) .ne. 0)
           j_prev = j                                  
           j = grid_var%dwnstrm_index(j)      !index of the next downstream point.  j_prev now the most recent upstream cell
           k = k + 1
        end do
         
        if (grid_var%land_mask(j) .eq. 0) then  !ended at an ocean point.  previous land pt is the outlet cell
          grid_var%ocean_outlet(i) = j_prev
!$OMP ATOMIC 
          grid_var%upstrm_number(j_prev) = grid_var%upstrm_number(j_prev) + 1
        elseif (grid_var%land_mask(j) .eq. 1) then          !ended at a land point.  this is the outlet cell
          grid_var%ocean_outlet(i) = j
!$OMP ATOMIC
          grid_var%upstrm_number(j) = grid_var%upstrm_number(j) + 1
        else
          stop
        end if
         
      else
        grid_var%ocean_outlet(i) = j                               !if it is ocean it is its own outlet
        grid_var%upstrm_number(j) = 0                              !and has only itself as an upstream cell
      end if
    end do
!$OMP END PARALLEL DO
         
    !count the total number of basins.   

    write(*,*) maxval(grid_var%upstrm_number(:))

    grid_var%nrr_cells = 0
    grid_var%nbasins   = 0
    do i=1,grid_var%npts
      if (grid_var%upstrm_number(i) .gt. 0) then   !count ocean cells that are now single cell basins.  dealt with in reorder_grid_by_basin
         grid_var%nbasins = grid_var%nbasins + 1
         grid_var%nrr_cells = grid_var%nrr_cells + grid_var%upstrm_number(i)
      end if
    end do

    write(*,*) 'ended associate_ocean_outlet_pts'
    
  end subroutine associate_ocean_outlet_points    

!----------------------------------------------------------------------------!
!  not going to use this
! going to go it on 2d grid just like 2d gw to make it easier
! both have same haloes

  subroutine reorder_grid_by_basin(grid_var,basins)
  
    implicit none
    type(river_grid_type), pointer, intent(inout) :: grid_var    
    type(basin_type)     , pointer, intent(inout) :: basins

    integer :: cnt, i,ii,j,jj,k,kk, total_nbasins, total_land_cells, ncells, partial_nbasins
    integer, allocatable, dimension(:)   :: tmp_indices
    integer, allocatable, dimension(:)   :: basin_map  !maps the basin number to the outlet cell number in the global array
    integer, allocatable, dimension(:)   :: basin_num_points
!    integer, allocatable, dimension(:,:) :: basin_points
    
    type(river_grid_type),pointer :: ord_grid_var   !grid variable ordered so basins are continuous
   
    write(*,*) 'in reorder'    
 
    allocate(ord_grid_var)
    call create(ord_grid_var, grid_var%npts)  !this is the total number of possible river points
   
    write(*,*) 'created ord_grid_var' 
    ord_grid_var%npts = grid_var%npts
    ord_grid_var%nlon = grid_var%nlon
    ord_grid_var%nlat = grid_var%nlat        
    ord_grid_var%nrr_cells = grid_var%nrr_cells
    ord_grid_var%nbasins   = grid_var%nbasins    
    
!    allocate(tmp_indices(grid_var%npts))   !tmp_indices large enough to hold all points rather than reallocating specific size

    !basin numbering is based on global array index.  not coninuous.  num basins < possible index vals
    total_nbasins = maxval(grid_var%ocean_outlet(:))!grid_var%nbasins  
    !total_nbasins  = maxval(grid_var%basin_ind(:))

    allocate(basin_num_points(total_nbasins))
    basin_num_points(:) = 0
    write(*,*) 'find basin_num pts'
    do k=1,grid_var%npts
      j = grid_var%ocean_outlet(k)
      !j = grid_var%basin_ind(k)
      basin_num_points(j) = basin_num_points(j) + 1
    end do 

    write(*,*) 'count basins'
    partial_nbasins = 0
    do i=1,total_nbasins
      if (basin_num_points(i) .gt. 2) then
        partial_nbasins = partial_nbasins + 1
      end if
    end do  

!    allocate(basin_points(total_nbasins,maxval(basin_num_points)))
!
!    write(*,*) total_nbasins,grid_var%npts

!    basin_num_points(:) = 0
!    basin_points(:,:)   = 0
!    do k=1,grid_var%npts
!      if (grid_var%upstrm_number(k) .gt. 0) then
!        j = grid_var%ocean_outlet(k)
!        basin_num_points(j) = basin_num_points(j) + 1
!        basin_points(j,basin_num_points(j)) = k
!      end if 
!    end do 
    write(*,*) 'partial nbasins - ',partial_nbasins

    ord_grid_var%nbasins = partial_nbasins
    if (associated(basins)) write(*,*) 'basins is already allocated'
    allocate(basins)

    write(*,*) 'toal_nbasins-',total_nbasins

    j = 0
    do i=1,total_nbasins  !count number of basins we include
      if (basin_num_points(i) .ge. min_pts_per_basin) then
        j = j + 1
      end if  
    end do 

    !allocate with largest basin to hole all.
    !wastes space but better than allocating thousands of derived type basins
    call alloc_basin(basins,j,maxval(basin_num_points))  

    !write(*,*) 'num points per basin:'  
    total_land_cells = 0
    j = 0
    do i=1,total_nbasins      !basin_num_points contains basins we don't include.  loop over total basins not partial
      cnt = basin_num_points(i)
      if (cnt .gt. min_pts_per_basin) then
        j = j + 1
        !write(*,*) 'org basin-',i,'new basin-',j,' #-',cnt
        !call alloc_basin(basins(j),cnt)
        basins%n_basin_cells(j)   = cnt

        kk=0
        do k=1,grid_var%npts
          jj = grid_var%ocean_outlet(k)
          !jj = grid_var%basin_ind(k)
          if (jj .eq. i) then
            kk = kk + 1
            basins%river_points(kk,j) = k
          end if
        end do

        total_land_cells = total_land_cells + cnt
      end if 
    end do 
    !sanity check that j == ord_grid_var%nbasins (ie partial_nbasins)

    !the ocean cells have been eliminated through the above process.  npts should also change
    
        ! I need to reorder the basin cells so they are contiguous.
        !then each loop can go from i=pe_start,pe_end if each holds global array
    write(*,*) 'start reordering'
    cnt=1
    do i=1,ord_grid_var%nbasins
    
      basins%begind(i) = cnt
      basins%endind(i) = cnt + basins%n_basin_cells(i) - 1
      
      do kk=1,basins%n_basin_cells(i)
      
        k = basins%river_points(kk,i)
        call copy_river_grid_vector_values(grid_var,ord_grid_var,cnt,cnt,k,k)  !copies only a single point for all vectors from grid_var to ord_grid_var
        !Note: flow variables yet to be defined.  no need to remap
        cnt = cnt + 1
        
      end do
      
    end do 
    write(*,*) 'reordered'
                                                       
    deallocate(basin_num_points)

    !destroy grid var.  make it new with fewer points (doesn't include the ocean now)
    deallocate(grid_var)
    grid_var => ord_grid_var

  end subroutine reorder_grid_by_basin
  
!----------------------------------------------------------------------------!

  subroutine remove_inactive_land_basins(grid_var,basins)
  
    implicit none
    
    type(river_grid_type),pointer, intent(inout) :: grid_var
    type(basin_type),     pointer, intent(inout) :: basins  
    
    type(basin_type), pointer     :: cmp_basins
    type(basin_type), pointer     :: tmp_basins
    type(river_grid_type),pointer :: cmp_grid_var
    
    integer :: i,j,k,ii,kk,cnt,js,je,ks,ke
    integer :: n_active_cells
    integer :: total_active_cells
    
    integer, allocatable, dimension(:) :: active_basin,basin_ind_map
    
    allocate(active_basin(grid_var%nbasins))
    active_basin(:) = 0    
   
   
    !  determine if each basin is active.  only active if all cells covered by land model.  make cut off 90%???
    total_active_cells = 0
    do i=1,grid_var%nbasins
      
      n_active_cells = sum(grid_var%active_cell(basins%begind(i):basins%endind(i)))  !count the number of active cells in the basin
      
      if (n_active_cells .ge. int(0.75*basins%n_basin_cells(i))) then   !compute basin if we have forcing for > 3/4 of it
        active_basin(i) = 1
        total_active_cells = total_active_cells + n_active_cells
      else
        active_basin(i) = 0
      end if
    end do

    allocate(cmp_grid_var) 
    call create(cmp_grid_var,total_active_cells)
    write(*,*) 'created cmp_grid_var'
    cmp_grid_var%nrr_cells = total_active_cells
    cmp_grid_var%npts      = total_active_cells
    cmp_grid_var%nlat      = grid_var%nlat
    cmp_grid_var%nlon      = grid_var%nlon
    cmp_grid_var%nbasins = sum(active_basin(:))   

    allocate(cmp_basins)
    call alloc_basin(cmp_basins,cmp_grid_var%nbasins,maxval(basins%n_basin_cells(:)))
    allocate(basin_ind_map(cmp_grid_var%nbasins))

    k=0
    do i=1,grid_var%nbasins
       if (active_basin(i) .eq. 1) then
         k = k + 1
         basin_ind_map(k) = i
       end if
    end do

    write(*,*) 'created cmp_basins with ',cmp_grid_var%nbasins,' basins out of ',grid_var%nbasins,' total basins'
    write(*,*) 'cmp_grid_var -> npts=',cmp_grid_var%npts
    write(*,*) 'grid_var -> npts=',grid_var%npts
    write(*,*) 'number of active cells by basin --'
    do i=1,size(active_basin)
      if (active_basin(i) .eq. 1) then
        write(*,*) 'basin number ',i,' with ',sum(grid_var%active_cell(basins%begind(i):basins%endind(i))),' active cells'
      end if
    end do
    
    if (cmp_grid_var%nbasins .lt. grid_var%nbasins) then   !there are some basins that aren't active
    
      cnt = 1    !keep track of starting index of current basins in the global vector of river cells
      do k=1,cmp_grid_var%nbasins
        i = basin_ind_map(k)
        cmp_basins%begind(k) = cnt
        cmp_basins%endind(k) = cnt + basins%n_basin_cells(i) - 1
        cmp_basins%n_basin_cells(k) = basins%n_basin_cells(i)
        !use temporaries to make code shorter
        js = cmp_basins%begind(k)   !j start
        je = cmp_basins%endind(k)    !j end
        
        ks = basins%begind(i)
        ke = basins%endind(i)

        !overloading = operator would make this a little cleaner
        call copy_river_grid_vector_values(grid_var,cmp_grid_var,js,je,ks,ke) 
        
        cnt = cmp_basins%endind(k) + 1
      
      end do

      write(*,*) 'finished copying to cmp_grid_var'
    
      !remove original grid_var variable.  reallocate new one with only active routing cells
      !call destroy(grid_var)
      write(*,*) 'destroy grid_var'
      !deallocate(grid_var)
      grid_var => cmp_grid_var

      !call create_river_grid_copy(cmp_grid_var,grid_var,total_active_cells) 
      write(*,*) 'called create copy'

      grid_var%nrr_cells = total_active_cells
      grid_var%active_cell(:) = 1  !all cells now active.
    
      !do the same for the basin variable
      write(*,*) 'destroy basins'

      write(*,*) 'point basins to cmp_basins'
      basins => cmp_basins
      write(*,*) 'done the poiting yo'


    else

      grid_var => cmp_grid_var


    end if  !some basins are not active.
    
    if (allocated(active_basin)) deallocate(active_basin)

    write(*,*) 'leaving remove_inactive_land_basins'
    
  end subroutine remove_inactive_land_basins

!----------------------------------------------------------------------------! 
  subroutine find_main_river_channels(grid_var)
    implicit none

    type(river_grid_type), intent(inout) :: grid_var

    integer :: k,kk,j
    integer :: ntot
    logical :: keep_looping

    ntot = size(grid_var%source_area(:))


    do k=1,ntot
      if (grid_var%source_area(k) .ge. source_area_channel_cutoff) then
        grid_var%is_main_channel = 1
      else
        grid_var%is_main_channel = 0
      end if
    end do

    !check that main channels are continuous
    do k=1,ntot
      write(*,*) real(k)/real(ntot)
      keep_looping = .true.
      if (grid_var%is_main_channel(k) .eq. 1) then
        j  = k
        kk = 1
        do while (keep_looping)

          j = grid_var%dwnstrm_index(j)
          kk = kk + 1
          if (kk .ge. ntot) then
            keep_looping = .false.
          end if

          if (grid_var%is_main_channel(j) .ne. 1) then
            write(*,*) 'a main channel become a non main channel.  something is wrong with main channel or down stream index'
            stop
          end if

        end do

      end if
    end do 

  end subroutine find_main_river_channels

!----------------------------------------------------------------------------! 
  !Determine the wieghts and mapping from 1D global land grid to the 1D river grid
  !river_grid%{lat,lon} are assumed higher resolution than lat_out lon_out?
  !assumes river grid is lat/lon and so is lsm grid
  
  subroutine determine_hilo_res_mapping(grid_var,lat_lo,lon_lo,pft_frac_lo)
  
    implicit none
    
    type(river_grid_type),  intent(inout) :: grid_var
    real(r_2), dimension(:), intent(in)    :: lat_lo
    real(r_2), dimension(:), intent(in)    :: lon_lo  !the lo resolution grid
    real(r_2), dimension(:), intent(in)    :: pft_frac_lo   !for tiled it is the fraction fo grid cell occupied by the pft

    !local variables
    integer   :: i,j,ii,k,kk   !integer counters for the loops

    real(r_2) :: dlat_lo,dlon_lo   !grid cell size (degrees) of lo res grid
    real(r_2) :: dlat_hi,dlon_hi !grid cell size (degrees) of hi res grid

    real(r_2) :: dlone,dlonw,dx,dlats,dlatn,dy   !cell size (degrees) and m (dx) in hi-lo  cells that overlap
    
    real(r_2) :: dy_lo,dy_hi,area_lo,area_hi  !size of the north-south side of the lo resolution grid (lo), area of the lo res grid
    
    real(r_2) :: Eedge_hi, Eedge_lo   !eastern  edge (lon) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Wedge_hi, Wedge_lo   !western  edge (lon) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Sedge_hi, Sedge_lo   !southern edge (lat) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Nedge_hi, Nedge_lo   !northern edge (lat) of high resolution (hi) and lo res (lo) grid cells

    !are lat and lon global 1D vectors or only vectors with lat lon values?
    if (size(grid_var%lat) .gt. grid_var%nlat .and. size(grid_var%lon) .gt. grid_var%nlon) then
       dlat_hi  = abs(grid_var%lat(grid_var%nlon+1)-grid_var%lat(1))
       dlon_hi  = abs(grid_var%lon(2) - grid_var%lon(1))
    else
       dlat_hi  = abs(grid_var%lat(2) - grid_var%lat(1))
       dlon_hi  = abs(grid_var%lon(2) - grid_var%lon(1))
    end if

    if (size(lat_lo) .gt. mlat .and. size(lon_lo) .gt. mlon) then
       dlat_lo  = abs(lat_lo(mlon+1)-lat_lo(1) )  !array is compacted,  must use lat_all lon_all brought in from io_vars module
       dlon_lo  = abs(lon_lo(2) - lon_lo(1) )
    else
       dlat_lo  = abs(lat_lo(2) - lat_lo(1))
       dlon_lo  = abs(lon_lo(2) - lon_lo(1))
    end if

    !init overlap mappings and areas
    grid_var%maps%n_ovrlap_lgr(:) = 0
    grid_var%maps%weight_lgr(:,:) = 0._r_2
    grid_var%maps%weight_rgl(:,:) = 0._r_2
    grid_var%maps%ind_lgr(:,:)    = -1
    
    !initialize all cells to inactive
    grid_var%active_cell(:)    = 0

    !map each input lat/lon to output grid and calc the fraction of cell within (1 cell can map to upto four output)
    do k=1,mp   !loop over all the land grid points
    
      Sedge_lo = lat_lo(k)  - dlat_lo/2._r_2
      Nedge_lo = lat_lo(k)  + dlat_lo/2._r_2    
      
      Wedge_lo = lon_lo(k)  - dlon_lo/2._r_2
      Eedge_lo = lon_lo(k)  + dlon_lo/2._r_2
    
      dy_lo = sin(deg2rad*Nedge_lo) - sin(deg2rad*Sedge_lo)
      area_lo = dy_lo * dlon_lo * re * re * pft_frac_lo(k) !adjust for fraction of grid cells?? need to for sure
      
      do kk=1,grid_var%npts

        Sedge_hi = grid_var%lat(kk) - dlat_hi/2._r_2  !avoid some calcs dooing it here
        Nedge_hi = grid_var%lat(kk) + dlat_hi/2._r_2
        
        dy_hi = sin(deg2rad*Nedge_hi) - sin(deg2rad*Sedge_hi)
        area_hi = dy_hi * dlon_hi * re * re  !adjust for fraction of grid cells??

        Wedge_hi = grid_var%lon(kk) - dlon_hi/2._r_2
        Eedge_hi = grid_var%lon(kk) + dlon_hi/2._r_2

        if ((Wedge_hi .le. Eedge_lo) .and. (Eedge_hi .ge. Wedge_lo) .and. &
            (Sedge_hi .le. Nedge_lo) .and. (Nedge_hi .ge. Sedge_lo)) then
              
          grid_var%maps%n_ovrlap_lgr(kk) = grid_var%maps%n_ovrlap_lgr(kk) + 1  !number of overlapping cells
          grid_var%active_cell(kk) = 1             !mark this river cell as active.  only active if covered by the lsm
          
          dlone = min(Eedge_lo,Eedge_hi)*deg2rad !determine area of input cell within the output grid cell
          dlonw = max(Wedge_lo,Wedge_hi)*deg2rad 
          dx = max(0.0,(dlone-dlonw))

          dlatn = min(Nedge_lo,Nedge_hi)*deg2rad 
          dlats = max(Sedge_lo,Sedge_hi)*deg2rad 
          dy = max(0.0,(sin(dlatn)-sin(dlats)))

          grid_var%maps%ind_lgr(kk,grid_var%maps%n_ovrlap_lgr(kk))    = k                   !lsm point for given river point
          grid_var%maps%weight_lgr(kk,grid_var%maps%n_ovrlap_lgr(kk)) = dx*dy / area_lo     !fraction of river cell k overlapping by lsm cell k.  need pft frac somehwere
          grid_var%maps%weight_rgl(kk,grid_var%maps%n_ovrlap_lgr(kk)) = dx*dy / area_hi

        end if  !test lon overlap
        
      end do !kk loop over input lat
      
    end do  !loop over mp land points
    
  end subroutine determine_hilo_res_mapping
!----------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------   
  subroutine step_river_routing(river,grid_var)
     implicit none
     
    type(river_flow_type),          intent(inout) :: river   !contains mass,flow variables
    type(river_grid_type),          intent(in)    :: grid_var
     
    integer :: kk_begind, kk_endind
    integer :: i,j,ii,kk
       
    !do i=basins_pe_start,basins_pe_end!    loop over a subsection of all of the basins

    !  kk_begind  = basins(i)%begind
    !  kk_endind  = basins(i)%endind
      kk_begind = 1
      kk_endind = grid_var%nrr_cells
    

      river%Fin(kk_begind:kk_endind) = 0._r_2   !zero out input fluxes

      do kk = kk_begind, kk_endind

        river%Fout(kk) = river%mass(kk)*river%vel(kk)/grid_var%length(kk)

        river%Fin(grid_var%dwnstrm_index(kk)) = river%Fin(grid_var%dwnstrm_index(kk)) + river%Fout(kk)

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
    
    real(r_2), allocatable, dimension(:)  :: river_fin_n, river_mass_n  !overland / river components
    real(r_2), allocatable, dimension(:)  :: subsurf_fin_n, subsurf_mass_n  !subsurface store, fluxes
    
    !do i=basins_pe_start,basins_pe_end!    loop over a subsection of all of the basins

    !  kk_begind  = basins(i)%begind
    !  kk_endind  = basins(i)%endind
    
      kk_begind = 1                !loop over all points for this mpi task.
      kk_endind = grid_var%nrr_cells
           
      allocate(river_fin_n(kk_endind))
      allocate(river_mass_n(kk_endind))   
      allocate(subsurf_mass_n(kk_endind))
      allocate(subsurf_fin_n(kk_endind))

      river_fin_n(:) = 0._r_2
      river_mass_n(:) = 0._r_2
      subsurf_mass_n(:) = 0._r_2
      subsurf_fin_n(:) = 0._r_2

      do kk = kk_begind, kk_endind
      
        k = grid_var%dwnstrm_index(kk)
        
        if (grid_var%is_main_channel(kk) .ne. 0 .and. grid_var%is_main_channel(k) .ne. 0) then  !current and downstream are main channels
          !main channel flow.  this will use manning and a rectangular channel eventually.  kinematic place holder.
          river_mass_n(kk) = (1.-river_theta) * river%mass(kk) + river%Fin(kk) +&
                             river%srf_runoff_flux(kk) + river%sub_runoff_flux(kk)
          river_fin_n(k)  = river_fin_n(k) + river_theta * river_mass_n(kk)          
          !subsurf
          subsurf_mass_n(kk) = 0.  !no subsurface for the main channel
          subsurf_fin_n(kk) = 0.
        elseif (grid_var%is_main_channel(k) .ne. 0) then    !not a main channel.  overland/subsurface routing for current, downstream is a main channel
          !overland flow
          river_mass_n(kk) = (1.-river_theta) * river%mass(kk) + river%Fin(kk) + river%srf_runoff_flux(kk)
          river_fin_n(k)  = river_fin_n(k) + river_theta * river_mass_n(kk)
          !subsurface flow.  leaves subsurface into main river channel
          subsurf_mass_n(kk) = (1. - subsurf_theta) * river%mass_sub(kk) + &
                                river%Fin_sub(kk) + river%sub_runoff_flux(kk)
          river_fin_n(k)   = river_fin_n(k) + subsurf_theta * subsurf_mass_n(kk)
          
        else   !both are not main channels
          !overland flow
          river_mass_n(kk) = (1.-river_theta) * river%mass(kk) + river%Fin(kk) + river%srf_runoff_flux(kk)
          river_fin_n(k)  = river_fin_n(k) + river_theta * river_mass_n(kk)
          !subsurface flow.  leaves subsurface into main river channel
          subsurf_mass_n(kk) = (1. - subsurf_theta) * river%mass_sub(kk) + river%Fin_sub(kk) + river%sub_runoff_flux(kk)
          subsurf_fin_n(k)   = subsurf_fin_n(k) + subsurf_theta * subsurf_mass_n(kk)
          
        end if
        
      end do
      
      do kk = kk_begind, kk_endind
      
        river%Fin(kk)  = river_fin_n(kk)  !store new values
        river%Fin_sub(kk) = subsurf_fin_n(kk)
        
        river%vel(kk) = river_theta / del_river * river%mass(kk)   !should be mass flux not vel
        river%mass(kk) = river_mass_n(kk)
        river%mass_sub(kk) = subsurf_mass_n(kk)
      
      end do                        

!        call basins(i)%get_outflow(river)  !write a subroutine to compute total basin outflow

    !end do  !loop over this pe basins
    
      deallocate(river_fin_n)
      deallocate(river_mass_n)
      deallocate(subsurf_fin_n)
      deallocate(subsurf_mass_n)
     
    end subroutine step_river_srf_subsrf_routing_kinematic
    
!----------------------------------------------------------------------------!

!**********Below are the Routines For Allocating/Deallocating****************!

  subroutine alloc_river_flow(var,npts)
    implicit none
    type(river_flow_type),intent(inout)  :: var
    integer, intent(in) :: npts

    allocate(var%mass_init(npts))
    var%mass_init(:) = fNaN
    
    allocate(var%mass(npts))
    var%mass(:) = fNaN
    
    allocate(var%mass_sub(npts))
    var%mass_sub(:) = fNaN    
    
    allocate(var%vel(npts))
    var%vel(:) = fNaN
    
    allocate(var%hgt(npts))
    var%hgt(:) = fNaN    
    
    allocate(var%srf_runoff_flux(npts))
    var%srf_runoff_flux(:) = fNaN
    
    allocate(var%sub_runoff_flux(npts))
    var%sub_runoff_flux(:) = fNaN    
    
    allocate(var%Fin(npts))
    var%Fin(:) = fNaN
    
    allocate(var%Fin_sub(npts))
    var%Fin_sub(:) = fNaN       
    
    allocate(var%Fout(npts))
    var%Fout(:) = fNaN       

  end subroutine alloc_river_flow
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_flow(var)
    implicit none
    type(river_flow_type),intent(inout) :: var

    deallocate(var%mass_init)
    deallocate(var%mass)
    deallocate(var%mass_sub)
    deallocate(var%vel)
    deallocate(var%hgt)
    deallocate(var%srf_runoff_flux)
    deallocate(var%sub_runoff_flux)    
    deallocate(var%Fin)
    deallocate(var%Fin_sub)
    deallocate(var%Fout)

  end subroutine dealloc_river_flow
  
!----------------------------------------------------------------------------!

  subroutine alloc_river_map_grid(var,npts)
    implicit none
    type(map_grid_type),intent(inout)  :: var
    integer, intent(in)                :: npts

    allocate(var%ind_lgr(npts,max_n_ovrlap))  !outut land cell number form a given river cell number
    var%ind_lgr(:,:) = iNaN
    
    allocate(var%weight_lgr(npts,max_n_ovrlap))
    var%weight_lgr(:,:) = fNaN
    
    allocate(var%weight_rgl(npts,max_n_ovrlap))
    var%weight_rgl(:,:) = fNaN    
    
    allocate(var%n_ovrlap_lgr(npts))
    var%n_ovrlap_lgr(:) = iNaN    
    

  end subroutine alloc_river_map_grid
  
!----------------------------------------------------------------------------!

  subroutine dealloc_river_map_grid(var)
    implicit none
    type(map_grid_type), intent(inout) :: var
    write(*,*) 'deallocating maps' 
    deallocate(var%ind_lgr)    
    deallocate(var%weight_lgr)    
    deallocate(var%weight_rgl)     
    deallocate(var%n_ovrlap_lgr)    
    write(*,*) 'deallocated maps' 
  end subroutine dealloc_river_map_grid
  
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
!    allocate(var%topo_in(npts))
!    var%topo_ind(:) = fNaN
!    allocate(var%basin_ind(npts))
!    var%basin_ind(:) = iNaN

    !since we know the total number of points set it
    var%npts = npts

    call create(var%maps,npts)         


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
    call destroy(var%maps)
    write(*,*) 'dealloced grid'
  end subroutine dealloc_river_grid

!----------------------------------------------------------------------------!

  subroutine alloc_basin(var,nbasins,npts)
    implicit none
    type(basin_type), intent(inout) :: var
    integer,          intent(in)    :: nbasins 
    integer,          intent(in)    :: npts
    
    allocate(var%begind(nbasins))
    var%begind(:) = iNaN
    allocate(var%endind(nbasins))
    var%endind(:) = iNaN
    allocate(var%n_basin_cells(nbasins))
    var%n_basin_cells(:) = iNaN 
    allocate(var%river_points(npts,nbasins))
    var%river_points(:,:) = iNaN
    
  end subroutine alloc_basin
  
!----------------------------------------------------------------------------!

  subroutine dealloc_basin(var)
    implicit none
    type(basin_type), intent(inout) :: var
    
    deallocate( var%begind )
    deallocate( var%endind )
    deallocate( var%n_basin_cells )
    deallocate( var%river_points )
    
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
    integer, dimension(:)  , allocatable :: mask_1d
    integer                              :: npts
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1
    npts    = nlon_rr*nlat_rr
    allocate(int_data(nlon_rr,nlat_rr))
    
    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,int_data(:,:),start_inds,end_inds) )


    var_data = reshape(int_data,(/npts/))

    if (present(mask)) then
      allocate(mask_1d(npts))
      mask_1d = reshape(mask,(/npts/))
      where(mask_1d .ge. 1) var_data(:) = -9999
    end if
    
    deallocate(int_data)
    if (allocated(mask_1d)) deallocate(mask_1d)
     
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
    
    integer :: nlat_rr,nlon_rr,var_id
    integer :: i,j,k
    
    real(r_2), dimension(:,:), allocatable :: flt_data
    integer  , dimension(:)  , allocatable :: mask_1d
    integer                                :: npts
    
    nlat_rr = end_inds(2) - start_inds(2) + 1
    nlon_rr = end_inds(1) - start_inds(1) + 1
    npts    = nlon_rr*nlat_rr
    allocate(flt_data(nlon_rr,nlat_rr))
    
    call check_nc( nf90_inq_varid(ncid,var_name,var_id) )
    call check_nc( nf90_get_var(ncid,var_id,flt_data(:,:),start_inds,end_inds) )


    var_data = reshape(flt_data,(/npts/))

    if (present(mask)) then
      allocate(mask_1d(npts))
      mask_1d = reshape(mask,(/npts/))
      where(mask_1d .ge. 1) var_data(:) = -9999
    end if

    deallocate(flt_data)
    if (allocated(mask_1d)) deallocate(mask_1d)

     
  end subroutine read_nc_flt 
  
!----------------------------------------------------------------------------!   
!****Below are the Routines For MPI settings/gettings/partitioning***********!
!----------------------------------------------------------------------------!  

#ifdef RRM_MPI

  subroutine calculate_basins_per_pe(grid_var,basins,my_basin_start,my_basin_end,&
                                     global_index_start,global_index_end,np,mrnk)
  !all procs call this routine.
  !all calc it but only save the basin start and end inds for myrank
  !or I could run only on myrank = 1, bcast to other ranks
    implicit none
    type(river_grid_type),                       intent(inout) :: grid_var
    type(basin_type), dimension(:), allocatable, intent(inout) :: basins
    integer, intent(out)                  :: my_basin_start     !number of starting basin to do 
    integer, intent(out)                  :: my_basin_end       !basin number of ending basins
    integer, intent(out)                  :: global_index_start !starting in dex in global river array
    integer, intent(out)                  :: global_index_end   !endind index in global river array
    integer, optional,intent(in)          :: np     !tmp name to check if present of nprocs
    integer, optional,intent(in)          :: mrnk   !myrank tmp name to allow to check if present
    
    integer, dimension(:), allocatable    :: npts_per_pe
    integer, dimension(:), allocatable    :: basins_pe_start
    integer, dimension(:), allocatable    :: basins_pe_end    
    integer                               :: ideal_npts_per_pe
    
    integer :: nprocs,myrank,nworkers
    
    integer :: bg,ed,i,j,k
    integer :: current_basin
    integer :: start_basin
    integer :: end_basin
    integer :: npts_tmp
    logical :: keep_searching
    
    if (present(np)) then
      nprocs = np
    else
      nprocs = 1
    end if
    
    if (present(mrnk)) then
      myrank = mrnk
    else
      myrank = 0
    end if
    
    nworkers = nprocs - 1   !master does no work following CABLE mpi model
    
    allocate(npts_per_pe(0:nworkers))
    npts_per_pe(:) = 0   
    
    allocate(basins_pe_start(0:nworkers)) 
    basins_pe_start(:) = 0
    
    allocate(basins_pe_end(0:nworkers)) 
    basins_pe_end(:) = 0
    
    if (myrank .eq. 0) then
    
      if (nworkers .gt. 1) then   !2 procs means 1 worker 1 master.  worker does all computations
    
        ideal_npts_per_pe = ceiling(real(grid_var%npts)/real(nworkers))
      
        current_basin = 1
      
        do i=1,nworkers
      
          keep_searching = .true.
          npts_tmp = 0
          start_basin = current_basin
        
          do while (keep_searching)
            npts_tmp = npts_tmp + basins(current_basin)%n_basin_cells
            if (npts_tmp .gt. 0.95*ideal_npts_per_pe) then
              keep_searching = .false.
            elseif (current_basin .lt. grid_var%nbasins -1) then
              current_basin = current_basin + 1
            else
              keep_searching = .false.
            end if
          end do
          end_basin = current_basin
          npts_per_pe(i) = npts_tmp
        
          basins_pe_start(i) = start_basin
          basins_pe_end(i)   = end_basin
        
          current_basin = current_basin + 1
          
        end do
      
        basins_pe_start(nworkers) = current_basin    !set last process to the rest
        basins_pe_end(nworkers)   = grid_var%nbasins
      
        !need to send this info to mpi workers
        do i=1,nworkers
          call MPI_SEND(basins_pe_start, nprocs, MPI_INTEGER, i, 0, comm, ierr)
          call MPI_SEND(basins_pe_end, nprocs, MPI_INTEGER, i, 0, comm, ierr)
        end do
      
      
      else
    
        basins_pe_start(1) = 1
        basins_pe_end(1)   = grid_var%nbasins
      
      end if
    

      write(*,*) 'The number of river cells per mpi proc is:'
      
      do i=1,nworkers
        write(*,*) 'proc: ',i,' number of cells: ',npts_per_pe(i)
      end do
      
      
    else
    
      call MPI_RECV(basins_pe_start,nprocs,MPI_INTEGER,0,0,comm,ierr)
      call MPI_RECV(basins_pe_end,nprocs,MPI_INTEGER,0,0,comm,ierr)
      
    end if  !rank 0 if
    
    my_basin_start = basins_pe_start(myrank)  !note myrank=0 is the master.  so set to 0
    my_basin_end   = basins_pe_end(myrank)
    
    global_index_start = basins(my_basin_start)%begind   !start and endind indices in the global river array
    global_index_end   = basins(my_basin_end)%endind
    
    
      
    deallocate(npts_per_pe)
    deallocate(basins_pe_start)
    deallocate(basins_pe_end)      
    

  end subroutine calculate_basins_per_pe
  
!----------------------------------------------------------------------------! 

  subroutine send_lsm_runoff_to_procs()
  end subroutine send_lsm_runoff_to_procs
  
  subroutine send_river_vars_to_procs()
  end subroutine send_river_vars_to_procs
  
  subroutine collect_river_vars_from_procs()
  end subroutine collect_river_vars_from_procs
  
  
  subroutine river_routing_main()   !called from cable.
  !if no mpi or rank=0
  !compact basins
  !dtermine basins per pe
  !start substepping
  !
  end subroutine river_routing_main
  
!end preproc check for mpi  
#endif

#endif

end module cable_routing