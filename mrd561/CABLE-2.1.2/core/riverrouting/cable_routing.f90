module cable_routing

  use netcdf
  use cable_types
  use cable_common_module, only : cable_user
  use cable_def_types_mod, only : r_2, ms, mp,mlat,mlon

  use cable_IO_vars_module, only : landpt, patch, max_vegpatches, metGrid,  &
                                  land_x,land_y,latitude,longitude,         &
                                  lat_all, lon_all,                         &
                                  logn,output,xdimsize,ydimsize,check,mask

  implicit none
  
  type map_grid_type
    integer, pointer, dimension(:,:)   :: ind_river_to_land
    real(r_2), pointer, dimension(:,:) :: weight_river_to_land
    integer, pointer, dimension(:)     :: n_ovrlap   !number of rivers cells that over lap the land cell
  end type map_grid_type

  
  type river_grid_type    !will move to cable types eventually    
    integer                            :: nlat
    integer                            :: nlon
    integer                            :: npts
    integer                            :: nbasins
    integer                            :: nrr_cells
    
    real(r_2), pointer, dimension(:) :: lat   !1dimensional array stored as lat=0,lon={1,nlon}; lat=1,lon={1,nlon}, etc
    real(r_2), pointer, dimension(:) :: lon
    real(r_2), pointer, dimension(:) :: length
    real(r_2), pointer, dimension(:) :: slope
    real(r_2), pointer, dimension(:) :: elev
    
    integer,   pointer, dimension(:) :: land_mask
    integer,   pointer, dimension(:) :: dwnstrm_index
    integer,   pointer, dimension(:) :: ocean_outlet   !index of the ending grid cell
    integer,   pointer, dimension(:) :: upstrm_number  !number of grid cells upstream from the current

    type(map_grid_type) :: maps
    
  end type river_grid_type
  
  
  type river_flow_type     !one global.  make 2D and have more than one tracer??????
    real(r_2), pointer, dimension(:) :: vol     !river water volume
    real(r_2), pointer, dimension(:) :: mass    !river water mass
    real(r_2), pointer, dimension(:) :: vel     !river water speed
    real(r_2), pointer, dimension(:) :: hgt     !river water height 
    
  end type river_flow_type
  
  type basin_type
    integer                        :: begind         !basin index start in reordered global 1d array
    integer                        :: endind         !basin end index in global reordered array
    integer                        :: n_basin_cells  !number of river cells in the basin
    integer, pointer, dimension(:) :: river_points   !the indices for the basin in the unordered global 1d array
    
  end type basin_type
  
  
  integer,   parameter :: n_river_tracers = 1   !number of tracers...liq  -->carbon,nitrogen and ice??
  integer,   parameter :: max_n_ovrlap = 324   !assume at most 1x1 to 0.0625x0.0625 -->16x16  do 18x18=324 to be safe
  real(r_2), parameter :: re = 6371000.0              !radius of the earth (km)
  real(r_2), parameter :: deg2rad = 3.14159/180.0     !constant converts degrees to radians
  real(r_2), parameter :: eps=1e-5                    !tolerance parameter.  UNUSED NOW 

  real(r_2), parameter :: fNaN = -1e36
  integer, parameter   :: iNaN = -999999

  type(river_grid_type), TARGET, SAVE :: river_grid
  type(river_flow_type), TARGET, SAVE :: river
  type(basin_type), pointer, save, dimension(:) :: basins
  
  
  !outline
  !if timestep==0 --> 
  !  call river_route_init(filename,river_grid)
  !    call read_river_route_file(filename,river_grid)
  !                                          => grid%{nlat,nlon,npts,lat,lon,direction,length,slope,elevation,land_mask}
  !       !need to only read lat/lon for river over same area area as lsm lat/lon
  !
  !       call alloc_river_route_vars()  !-->call within
  !
  !    call find_downstream_index(river_grid)
  !                            => river_grid%{dstrm_index,river_dirc}
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
    !determine runoff from lsm to river cells:
    !Q_river_runoff(:) = 0._r_2
    !do k=1,mp
    !  do i=1,n_ovrlap(k)
    !    Q_river_runoff(map_hi_to_lo(k,i)) = Q_river_runoff(map_hi_to_lo(k,i))  + Q_runoff(k)*weight_hi_to_lo(k,i)
    !  end do
    !end do
    
    !map the river water to each grid cell
    !wb_river(:) = 0._r_2
    !do k=1,mp
    !  do i=1,n_ovrlap(k)
    !    wb_river_land(k) = wb_river_land(k) + river%wb(map_hi_to_lo(k,i))*weight_hi_to_lo(k,i)
    !  end do
    !end do    


  

contains
!-----------------------------------------------------------------------------
  subroutine read_river_route_data(filename,river_grid)
    use netcdf
    implicit none
    
    character(len=250), intent(in)       :: filename
    type(river_grid_type), intent(out) :: river_grid
    
    integer :: ncid_river
    integer :: mask_id, rdir_id, length_id, slope_id, elev_id
    integer :: rr_lat_dim_id,rr_lon_dim_id,rr_lat_var_id,rr_lon_var_id
    integer :: nlat_rr_file, nlon_rr_file, npts_rr_file
    
    
    character(len=*), parameter :: mask_name   = "land_mask"
    character(len=*), parameter :: length_name = "length"
    character(len=*), parameter :: slope_name  = "slope"
    character(len=*), parameter :: elev_name   = "elevation"
    character(len=*), parameter :: rdir_name   = "direction"
    character(len=*), parameter :: rr_lat_dim_name = "lat"
    character(len=*), parameter :: rr_lon_dim_name = "lon"
    character(len=*), parameter :: rdir_name   = "direction"
    character(len=*), parameter :: rr_lat_var_name = "latitude"
    character(len=*), parameter :: rr_lon_var_name = "longitude"
    
    integer  , dimension(:,:), allocatable :: int_data
    real(r_2), dimension(:,:), allocatable :: flt_data
    
    real(r_2), dimension(:), allocatable :: lat_data
    real(r_2), dimension(:), allocatable :: lon_data
    
    integer :: i,j,k     !integers to loop through lsm data
    integer :: ri,rj,rk  !integers to loop through river
    integer :: nc_check
    
    integer :: lon_start,lon_end,lat_start,lat_end
    
    real :: dlat_lsm,dlon_lsm
    
    
    nc_check = nf90_open(trim(filename), nf90_nowrite, ncid_river)
    
    !get dim ids
    nc_check = nf90_inq_dimid(ncid_river, rr_lat_dim_name, rr_lat_dim_id)
    nc_check = nf90_inq_dimid(ncid_river, rr_lon_dim_name, rr_lon_dim_id)
    
    !get dim lengths
    nc_check = nf90_inquire_dimension(ncid_river,rr_lat_dim_id,len=nlat_rr_file) 
    nc_check = nf90_inquire_dimension(ncid_river,rr_lon_dim_id,len=nlon_rr_file) 
    
    allocate(lat_data(nlat_rr_file))
    allocate(lon_data(nlon_rr_file))
    
    !read the latitude and longitudes
    nc_check = nf90_inq_varid(ncid_river,rr_lat_var_name,rr_lat_var_id)
    nc_check = nf90_inq_varid(ncid_river,rr_lon_var_name,rr_lon_var_id)
    
    !
    nc_check = nf90_get_var(ncid_river,rr_lat_var_id,lat_data(:))
    nc_check = nf90_get_var(ncid_river,rr_lon_var_id,lon_data(:))
    
    !determine regions covered by the lsm forcing/data/simulation
    !assume always all goes from W->E and S->N
    dlat_lsm = 0.5*abs(latitude(2)-latitude(1))
    dlon_lsm = 0.5*abs(longitude(2)-longitude(1))
    
    min_lsm_lat = minval(latitude(:)) - dlat_lsm
    max_lsm_lat = maxval(latitude(:)) + dlat_lsm
    min_lsm_lon = minval(longitude(:)) - dlon_lsm
    max_lsm_lon = maxval(longitue(:)) + dlon_lsm
    
    found_edge = .false.
    i=1
    do while (.not. found_edge)
      if (lon_data(i) .le. min_lsm_lon .and. lon_data(i+1) .ge. min_lsm_lon) then
        found_edge = .true.
      elseif (i .le. nlon_rr_file-1) then
        i=i+1
      else
         stop
         found_edge = .true.
    end do
    lon_start = i
    

    found_edge = .false.
    i=nlon_rr_file
    do while (.not. found_edge)
      if (lon_data(i-1) .le. max_lsm_lon .and. lon_data(i) .ge. max_lsm_lon) then
        found_edge = .true.
      elseif (i .ge. 2) then
        i=i-1
      else
       stop
       found_edge = .true.
    end do
    lon_end = i    
       


    
    
  end subroutine read_river_data
    


!-----------------------------------------------------------------------------
!             !current data set doesn't use these numbers as it is: 1,2,4,8,16,32,64,128,256
!                       32 64  128
!                       16     1
!                        8  4  2
!
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
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
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
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!  call find_downstream_index(river_grid%dstrm_index,river_dirc)
  subroutine find_downstream_index(river_grid)
    implicit none
    
    type(river_grid_type), intent(inout) :: river_grid
    
    integer :: i,j,k,ii,jj,kk
    
    river_grid%dstrm_index(:) = 0
    
    do j=1,river_grid%nlat
      do i=1,river_grid%nlon
      
        k = i + (j-1)*river_grid%nlon
        
        if (river_grid%direction(k) /= 0) then
        
          ii = i + dirc2latindex(river_grid%direction(k))
          jj = j + dirc2lonindex(rriver_grid%direction(k))
          
          if (ii .lt. 1     )          ii = ii + river_grid%nlon  
          if (ii .gt. river_grid%nlon) ii = ii - river_grid%nlon
          if (jj .lt. 1 .or. jj .gt. river_grid%nlat .or. ii .lt. 1 .or. ii .gt. river_grid%nlon) then
            stop
          endif
          
          kk = ii + (jj-1)*river_grid%nlat
          
          river_grid%dstrm_index(k) = kk

        endif
      enddo
    enddo

  end subroutine find_downstream_index
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------  
  subroutine associate_ocean_outlet_points(river_grid)
    implicit none
    
    type(river_grid_type), intent(inout) :: river_grid
  
    !local variables
    integer :: i,j,k
  
    river_grid%ocean_outlet(:) = -1
    
    do i=1,river_grid%npts
      j = i
      if (river_grid%land_mask(j) .eq. 1) then
        k = 0
         
        do while ((river_grid%land_mask(k) .eq. 1) .and. k < river_grid%npts)
           j = river_grid%dstrm_index(j)
           k = k + 1
        end do
         
        if (river_grid%land_mask(j) .eq. 2) then
          river_grid%ocean_outlet(i) = j
          river_grid%upstrm_number(j) = river_grid%upstrm_number(j) + 1
        elseif (land_mask(j) .eq. 1) then
          river_grid%ocean_outlet(i) = j
          river_grid%upstrm_number(j) = river_grid%upstrm_number(j) + 1
        else
          stop
        end if
         
      else
        river_grid%ocean_outlet(i) = j
        river_grid%upstrm_number(j) = river_grid%upstrm_number(j) + 1
      end if
    end do
         
    !count the total number of basins.
    river_grid%nrr_cells = 0
    river_grid%nbasins   = 0
    do i=1,river_grid%npts
      if (river_grid%upstrm_number(j) .ne. 0) then
         river_grid%nbasins = river_grid%nbasins + 1
         river_grid%nrr_cells = river_grid%nrr_cells + 1
      end if
    end do
    
  end subroutine associate_ocean_outlet_points    
!-----------------------------------------------------------------------------    
!-----------------------------------------------------------------------------  
  subroutine reorder_grid_by_basin(river_var,grid_var,basin)
    implicit none
    type(river_flow_type), dimension(:), intent(inout)   :: river_var
    type(river_grid_type), dimension(:), intent(inout)   :: grid_var
    type(basin_type), pointer, dimension(:), intent(out) :: basins
    !loop: basins
    !loop : find all for basin i
    !  create vector of indices for basin i --> output: basin_indices(nbasins,npts).
    !             basins won't have contiguous memory...hmmmm

    integer :: cnt, i, kk
    integer, allocatable, dimension(:) :: tmp_indices
    
    type(river_grid_type) :: ord_rid_var   !grid variable ordered so basins are continuous
    
    call alloc_river_grid(ord_grid_var)
    
    ord_grid_var%npts = grid_var%npts
    ord_grid_var%nlon = grid_var%nlon
    ord_grid_var%nlat = grid_var%nlat        
    ord_grid_var%nrr_cells = grid_var%nrr_cells
    ord_grid_var%nbasins   = grid_var%nbasins    
    
    allocate(tmp_indices(grid_var%npts))
    
    allocate(basins(grid_var%nbasins))

    do i=1,grid_var%nbasins
       tmp_indices(:) = 0
       cnt = 0
       do kk=1,grid_var%npts
          if (grid_var%ocean_outlet(kk) .eq. i) then
             cnt = cnt + 1
             tmp_indices(cnt) = kk
          end if
        end do
        
        ncells = cnt
        basins(i)%n_basin_cells = cnt
        allocate(basins(i)%river_points(cnt))
        basins(i)%river_points(:) = tmp_indices(1:cnt)  !i like this solution.  simply pass basin indices to loop over. 
                                                        !will need to put these in contiguous array to pass back to master.
    end do
        ! I need to reorder the basin cells so they are contiguous.
        !then each loop can go from i=pe_start,pe_end; wat_mass(i) = dt*(Fin-Fout) + wat_mass(i)
    cnt=1
    do i=1,grid_var%nbasins
      basins(i)%begind = cnt
      basins(i)%endind = cnt + basins(i)%n_basin_cells
      do kk=1,basins(i)%n_basin_cells
        k = basins(i)%river_points(kk)
        
        ord_grid_var%dwnstrm_index(cnt) = grid_var%dwnstrm_index(k)   !set ordered by basin to the non ordered values
        ord_grid_var%ocean_outlet(cnt)  = grid_var%ocean_outlet(k)
        ord_grid_var%upstrm_number(cnt) = grid_var%upstrm_number(k)
        ord_grid_var%maps%ind_river_to_land(cnt,:) = grid_var%maps%ind_river_to_land(k,:)
        ord_grid_var%maps%weight_river_to_land(cnt,:) = grid_var%maps%weight_river_to_land(k,:)
        ord_grid_var%maps%n_ovrlap(cnt,:) = grid_var%maps%n_ovrlap(k,:)
        ord_grid_var%lat(cnt)             = grid_var%lat(k)
        ord_grid_var%lon(cnt)             = grid_var%lon(k)
        ord_grid_var%slope(cnt)           = grid_var%slope(k)
        ord_grid_var%length(cnt)          = grid_var%length(k)
        ord_grid_var%land_mask(cnt)       = grid_var%land_mask(k)
        
        !flow variables yet to be defined.  no need to remap
        ! wat_mass, wat_vol, wat_hgt, wat_length
        cnt = cnt + 1
      end do
    end do                                                    
    deallocate(tmp_indices)
    
    !then set to grid_var which is now continuous in terms of basins.  map same as before reordered
    grid_var = ord_grid_var
    call dealloc_river_grid(ord_grid_var)

  end subroutine reorder_grid_by_basin
!----------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------  
  !Determine the wieghts and mapping from 1D global land grid to the 1D river grid
  !river_grid%{lat,lon} are assumed higher resolution than lat_out lon_out?
  subroutine determine_hilo_res_mapping(river_grid,lat_lo,lon_lo)
  
    implicit none
    
    type(river_grid_type), intent(inout) :: river_grid
    real(r_2), dimension(:),     intent(in)  :: lat_lo,lon_lo  !the lo resolution grid

    !local variables
    integer   :: i,j,ii,jj,k,kk,k_tmp,kk_tmp  !integer counters for the loops

    real(r_2) :: dlat_lo,dlon_lo   !grid cell size (degrees) of lo res grid
    real(r_2) :: dlat_hi,dlon_hi !grid cell size (degrees) of hi res grid

    real(r_2) :: dlone,dlonw,dx,dlats,dlatn,dy   !cell size (degrees) and m (dx) in hi-lo  cells that overlap
    
    real(r_2) :: dy_lo,area_lo  !size of the north-south side of the lo resolution grid (lo), area of the lo res grid
    
    real(r_2) :: Eedge_hi, Eedge_lo   !eastern  edge (lon) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Wedge_hi, Wedge_lo   !western  edge (lon) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Sedge_hi, Sedge_lo   !southern edge (lat) of high resolution (hi) and lo res (lo) grid cells
    real(r_2) :: Nedge_hi, Nedge_lo   !northern edge (lat) of high resolution (hi) and lo res (lo) grid cells

    !are lat and lon global 1D vectors or only vectors with lat lon values?
    if (size(river_grid%lat_hi) .gt. river_grid%nlat .and. size(river_grid%lon) .gt. river_grid%nlon) then
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
    river_grid%map%n_ovrlap(:)               = 0
    river_grid%map%weight_river_to_land(:,:) = 0._r_2
    river_grid%map%ind_river_to_land(:,:)    = -1

    !map each input lat/lon to output grid and calc the fraction of cell within (1 cell can map to upto four output)
    do k=1,mp   !loop over all the land grid points
    
      Sedge_lo = lat_lo(k)  - dlat_lo/2._r_2
      Nedge_lo = lat_lo(k)  + dlat_lo/2._r_2    
      
      Wedge_lo = lon_lo(k)  - dlon_lo/2._r_2
      Eedge_lo = lon_lo(k)  + dlon_lo/2._r_2
    
      dy_lo = sin(deg2rad*Nedge_lo) - sin(deg2rad*Sedge_lo)
      area_lo = dy_lo * dlon_lo * re * re  !adjust for fraction of grid cells??
      
      do jj=1,river_grid%nlat   !by using two loops Sedge_hi,Nedge_hi only recalced when needed

        kk_tmp = 1 + jj*river_grid%nlon

        Sedge_hi = river_grid%lat(kk_tmp) - dlat_hi/2._r_2  !avoid some calcs dooing it here
        Nedge_hi = river_grid%lat(kk_tmp) + dlat_hi/2._r_2
        
        dy_hi = sin(deg2rad*Nedge_hi) - sin(deg2rad*Sedge_hi)
        area_hi = dy_hi * dlon_hi * re * re  !adjust for fraction of grid cells??

        if ((Sedge_hi .le. Nedge_lo) .and. (Nedge_hi .ge. Sedge_lo)) then  !overlap in lat

          do ii=1,river_grid%nlon

            kk = ii + jj*river_grid%nlon

            Wedge_hi = river_grid%lon(kk) - dlon_hi/2._r_2
            Eedge_hi = river_grid%lon(kk) + dlon_hi/2._r_2

            if ((Wedge_hi .le. Eedge_lo) .and. (Eedge_hi .ge. Wedge_lo)) then
              
              river_grid%map%n_ovrlap(k) = river_grid%map%n_ovrlap(k) + 1  !number of overlapping cells

              dlone = min(Eedge_lo,Eedge_hi)*deg2rad !determine area of input cell within the output grid cell
              dlonw = max(Wedge_lo,Wedge_hi)*deg2rad 
              dx = max(0.0,(dlone-dlonw))

              dlatn = min(Nedge_lo,Nedge_hi)*deg2rad 
              dlats = max(Sedge_lo,Sedge_hi)*deg2rad 
              dy = max(0.0,(sin(dlatn)-sin(dlats)))

              river_grid%map%ind_river_to_land( k , river_grid%map%n_ovrlap(k) )    = kk              !all river points for a given land point
              river_grid%map%weight_river_to_land( k , river_grid%map%n_ovrlap(k) ) = dx*dy / area_hi !frac of hi res cell within this lo res cell

            end if  !test lon overlap
            
          end do  !ii lon loop over input
          
        end if  !test lat overlap
        
      end do !jj loop over input lat
      
    end do  !loop over mp land points
    
  end subroutine determine_hilo_res_mappings 
!----------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------   
  subroutine step_river_routing(river,river_grid,basins(:),basins_pe_start,basins_pe_end)
     implicit none
     
     type(river_flow_type),          intent(inout) :: river   !contains mass,flow variables
     type(river_grid_type),          intent(in)    :: river_grid
     type(basin_type), dimension(:), intent(in)    :: basins          !contains info on each basin     
     
     real(r_2), pointer, dimension(:) :: river_water
     real(r_2), pointer, dimension(:) :: river_vel
     real(r_2), pointer, dimension(:) :: river_hgt
     real(r_2), pointer, dimension(:) :: river_distance
     
     integer :: npoints_local_basin
     integer :: i
     
     river_water => river%mass(:)                            !pointers to full arrays.  each pe has full copy of array?
     river_vel   => river%vel(:)
     river_hgt   => river%hgt(:)
     river_distance => river_grid%distance(:)
     
     do i=basins_pe_start,basins_pe_end                       !loop over basins assinged to this pe
       npoints_local_basin = basins(i)%n_basin_cells
       do k=1,npoints_local_basin
         kk = basins(i)%river_points(k)
         Fin(kk) = 0._r_2
       end do
       
       do k=1,npoints_local_basin
         kk = basins(i)%river_points(k)
         Fout(kk) = river_water(kk)*river_vel(kk)/river_disctance(kk)/river_distance(kk)*dels  !made up
         Fin(grid_var%dwnstrm_index(kk)) = Fin(grid_var%dwnstrm_index(kk)) + Fout(kk)
       end do
       
       do k=1,npoints_local_basin
         kk = basins(i)%river_points(k)
         river_water(kk) = river_water(kk) + (Fin(kk)-Fout(kk))*dels
         if (grid_var%is_outlet .eq. 1) then
           total_outflow(i) = Fout(kk)*dels
         end if
       end do
     
     end do
     
   end subroutine step_river_routing
!-----------------------------------------------------------------------------  
!-----------------------------------------------------------------------------
  subroutine alloc_river_route_vars(river_var,grid_var,npts)
    implicit none
    type(river_flow_type), intent(inout) :: river_var
    type(river_grid_type), intent(inout) :: grid_var
    integer, intent(in)                  :: npts
    
    call alloc_river_flow(river_var,npts)
    call alloc_river_grid(gird_var,npts)
  end subroutine alloc_river_route_vars  
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine alloc_river_flow(var,npts)

    type(river_flow_type), intent(inout) :: var
    integer, intent(in) :: npts

    allocate(var%vol(npts))
    var%vol(:) = fNaN
    
    allocate(var%mass(npts))
    var%mass(:,:) = fNaN
    
    allocate(var%vel(npts))
    var%vel(:) = fNaN
    
    allocate(var%hgt(npts))
    var%hgt(:) = fNaN    

  end subroutine alloc_river_flow
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine dealloc_river_flow(var)

    type(river_flow_type), intent(inout) :: var

    deallocate(var%vol)
    deallocate(var%mass)
    deallocate(var%vel)
    deallocate(var%hgt)

  end subroutine dealloc_river_flow
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine alloc_river_map_grid(var,npts)
    type(map_grid_type), intent(inout) :: var
    integer, intent(in)                :: npts

    allocate(var%ind_river_to_land(npts,max_n_ovrlap))
    var%ind_river_to_land(:,:) = iNaN
    
    allocate(var%weight_river_to_land(npts,max_n_ovrlap))
    var%weight_river_to_land(:,:) = fNaN
    
    allocate(var%n_ovrlap(npts))
    var%n_ovrlap(:,:) = iNaN

  end subroutine alloc_river_map_grid
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine dealloc_river_map_grid(var)
    type(map_grid_type), intent(inout) :: var
    
    deallocate(var%ind_river_to_land)    
    deallocate(var%weight_river_to_land)    
    deallocate(var%n_ovrlap)
    
  end subroutine dealloc_river_map_grid
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine alloc_river_grid(var,npts)

    type(river_grid_type), intent(inout) :: var
    integer, intent(in)             :: npts

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

    call alloc_river_map_grid(var%maps,npts)


  end subroutine alloc_river_grid
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine dealloc_river_grid(var,npts)

    type(river_grid_type), intent(inout) :: var

    deallocate(var%lat)
    deallocate(var%lon)
    deallocate(var%slope)
    deallocate(var%length)
    deallocate(var%dwnstrm_index)
    deallocate(var%elev)
    deallocate(var%land_mask)
    
    call dealloc_river_map_grid(var%maps)

  end subroutine alloc_river_grid

!-----------------------------------------------------------------------------  


end module cable_routing

    
    
    
    
    
  








end module cable_routing
