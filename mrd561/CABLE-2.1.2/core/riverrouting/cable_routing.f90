module cable_routing

  use routing_constants, only :: nlat_river,nlon_river,npts_river
  use cable_types
  USE cable_common_module, ONLY : cable_user
  use cable_def_types_mod, only : r_2, ms, mp,mlat,mlon

  use cable_IO_vars_module, only : landpt, patch, max_vegpatches, metGrid,  &
                                  land_x,land_y,latitude,longitude,         &
                                  lat_all, lon_all,                         &
                                  logn,output,xdimsize,ydimsize,check,mask

  implicit none
  
  type map_grid_type
    integer, pointer, dimension(:,:)   :: ind_river_to_land
    real(r_2), pointer, dimension(:,:) :: weight_river_to_land
    integer, pointer, dimension(:)     :: n_overlap   !number of rivers cells that over lap the land cell
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
    
    integer,   pointer, dimension(:) :: dwnstrm_index
    integer,   pointer, dimension(:) :: ocean_outlet   !index of the ending grid cell
    integer,   pointer, dimension(:) :: upstrm_number  !number of grid cells upstream from the current

    type(map_grid_type) :: maps
    
  end type river_grid_type
  
  
  type river_flow_type     !will need one for each mpi basin?    
    real(r_2), pointer, dimension(:) :: wat_vol     !river water volume
    real(r_2), pointer, dimension(:) :: wat_mass    !river water mass
    real(r_2), pointer, dimension(:) :: wat_vel     !river water speed
    real(r_2), pointer, dimension(:) :: wat_hgt     !river water height 
  end type river_flow_type
  
  integer,   parameter :: max_n_overlap = 169   !assume at most 3x3 to 0.25x0.25 -->12x12  do 13x13=169 to be safe
  real(r_2), parameter :: re = 6371000.0              !radius of the earth (km)
  real(r_2), parameter :: deg2rad = 3.14159/180.0     !constant converts degrees to radians
  real(r_2), parameter :: eps=1e-5                    !tolerance parameter.  UNUSED NOW 

  real(r_2), parameter :: fNaN = -1e36
  integer, parameter   :: iNaN = -999999

  type(river_grid_type), SAVE :: river_grid
  type(river_flow_type), SAVE :: river
  
    


contains

!-----------------------------------------------------------------------------

  subroutine alloc_river_flow(var,npts)

    type(river_flow_type), intent(inout) :: var
    integer, intent(in) :: npts

    allocate(var%wat_vol(npts))
    var%wat_vol(:) = fNaN
    
    allocate(var%wat_mass(npts))
    var%wat_mass(:) = fNaN
    
    allocate(var%wat_vel(npts))
    var%wat_vel(:) = fNaN
    
    allocate(var%wat_hgt(npts))
    var%wat_hgt(:) = fNaN    

  end subroutine alloc_river_flow

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine dealloc_river_flow(var)

    type(river_flow_type), intent(inout) :: var

    deallocate(var%wat_vol)
    deallocate(var%wat_mass)
    deallocate(var%wat_vel)
    deallocate(var%wat_hgt)

  end subroutine dealloc_river_flow

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine alloc_river_map_grid(var,npts)
    type(map_grid_type), intent(inout) :: var
    integer, intent(in)                :: npts

    allocate(var%ind_river_to_land(npts,max_n_overlap))
    var%lat_overlap(:,:) = iNaN
    
    allocate(var%weight_river_to_land(npts,max_n_overlap))
    var%lon_overlap(:,:) = fNaN
    
    allocate(var%n_overlap(npts))
    var%frac_overlap(:,:) = iNaN

  end subroutine alloc_river_map_grid

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine dealloc_river_map_grid(var)
    type(map_grid_type), intent(inout) :: var
    
    deallocate(var%ind_river_to_land)    
    deallocate(var%weight_river_to_land)    
    deallocate(var%n_overlap)
    
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
    
    call dealloc_river_map_grid(var%maps)


  end subroutine alloc_river_grid

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!             !current data set doesn't use these numbers as it is: 1,2,4,8,16,32,64,128,256
!                8   1   2      128  1   2
!                7   0   3       64      4
!                6   5   4       32  16  8
!
!

  integer function dirc2latindex(in_dirc) 
    implicit none
    integer, intent(in) :: in_dirc

    dirc2latindex = 0
    if ((in_dirc .le. 2) .or.  (in_dirc .eq. 128)) dirc2latindex  = 1
    if ((in_dirc .ge. 8) .and. (in_dirc .le.  32)) dirc2latindex = -1

  end function dirc2latindex
!-----------------------------------------------------------------------------
  

!-----------------------------------------------------------------------------
!
  integer function dirc2lonindex(in_dirc) 
    implicit none
    integer, intent(in) :: in_dirc

    dirc2lonindex = 0
    if ((in_dirc .ge.  2) .and. (in_dirc .le.   8)) dirc2lonindex  = 1
    if ((in_dirc .ge. 32) .and. (in_dirc .le. 128)) dirc2lonindex = -1


  end function dirc2lonindex
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!
!  call find_downstream_index(river_grid%dstrm_index,river_dirc)

  subroutine find_downstream_index(dstrm_index,river_dirc)

    integer, dimension(:), intent(out) :: dstrm_index
    integer, dimension(:), intent(in ) :: river_dirc
    
    integer :: i,j,k,ii,jj,kk
    
    dstrm_index(:) = 0
    
    do j=1,river_grid%nlat
      do i=1,river_grid%nlon
      
        k = i + (j-1)*river_grid%nlon
        
        if (rdirc(k) /= 0) then
        
          ii = i + dirc2latindex(river_dirc(k))
          jj = j + dirc2lonindex(river_dirc(k))
          
          if (ii .lt. 1     )          ii = ii + river_grid%nlon  
          if (ii .gt. river_grid%nlon) ii = ii - river_grid%nlon
          if (jj .lt. 1 .or. jj .gt. river_grid%nlat .or. ii .lt. 1 .or. ii .gt. river_grid%nlon) then
            stop
          endif
          
          kk = ii + (jj-1)*river_grid%nlat
          
          dstrm_index(k) = kk

        endif
      enddo
    enddo

  end subroutine find_downstream_index
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------

  subroutine read_stream_length(dstrm_index)
  end subroutine find_stream_length

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!  
  subroutine map_cable2river(grid_i,grid_o,data_i,data_o,weights_i2o,mapindex_i2o)
    type(grid_type), intent(in)   :: grid_i
    type(grid_type), intent(in)   :: grid_o
    real(kind=r_2), intent(in)    :: data_i(:,:)
    real(kind=r_2), intent(out)   :: data_o(:,:)
    real(kind=r_2), intent(inout) :: weights_i2o(:,:)
    integer, intent(inout)        :: mapindex_i2o(:,:)

  end subroutine map_cable2river
!-----------------------------------------------------------------------------
  
!-----------------------------------------------------------------------------

  subroutine step_riverflow()

  end subroutine step_riverflow
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------

  subroutine map_river2cable()

  end subroutine map_river2cable
!-----------------------------------------------------------------------------

 
!-----------------------------------------------------------------------------

  subroutine mass_balance()
  end subroutine mass_balance

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------

  subroutine calc_lsm_to_river_mapping()

    type(grid_type), intent(inout)   :: river_grid



  end subroutine calc_lsm_to_river_mapping
  
  
!-----------------------------------------------------------------------------  
  
!Need to determine basins.  They are independnet and can be used to divide up 
!computations when using MPI  
! call associate_ocean_outlet_points(river_grid%ocean_outlet,river_grid%upstrm_number,river_grid%nbasins,river_grid%nrr_cells,land_mask)
  
  subroutine associate_ocean_outlet_points(ocean_outlet,upstrm_number,nbasins,nrr_cells,land_mask)
  
    integer,   dimension(:), intent(out) :: ocean_outlet     !number identifying the river end point
    integer,   dimension(:), intent(out) :: upstrm_number    !the number of cells upstream of a grid cell
    integer,                 intent(out) :: nbasins          !the total number of river basins.  these are independnet
    integer,                 intent(out) :: nrr_cells        !the total number of river routing cells
    integer,   dimension(:), intent(in)  :: land_mask        !mask that defines if a cell is land or ocean
  
    integer :: i,j,k
  
    ocean_outlet(:) = -1
    
    do i=1,npts
      j = i
      if (land_mask(j) .eq. 1) then
        k = 0
         
        do while ((land_mask(k) .eq. 1) .and. k < npts)
           j = dstrm_index(j)
           k = k + 1
        end do
         
        if (land_mask(j) .eq. 2) then
          ocean_outlet(i) = j
          upstrm_number(j) = upstrm_number(j) + 1
        elseif (land_mask(j) .eq. 1) then
          ocean_outlet(i) = j
          upstrm_number(j) = upstrm_number(j) + 1
        else
          stop
        end if
         
      else
        ocean_outlet(i) = j
        upstrm_number(j) = upstrm_number(j) + 1
      end if
    end do
         
    !count the total number of basins.
    nrr_cells = 0
    nbasins   = 0
    do i=1,npts
      if (upstrm_number(j) .ne. 0) then
         nbasins = nbasins + 1
         nrr_cells = nrr_cells + 1
      end if
    end do
    
  end subroutine associate_ocean_outlet_points    
    
  
    !or I can simply pass around parts of the large array
  subroutine partition_basins_for_mpi(river_var,grid_var,nbasins,basin)
    

    type(river_flow_type), dimension(:), intent(inout)   :: river_var
    type(river_grid_type), dimension(:), intent(inout)   :: grid_var
    integer,                             intent(in)      :: nbasins
    type(basin_type), pointer, dimension(:), intent(out) :: basins
    !loop: basins
    !loop : find all for basin i
    !  create vector of indices for basin i --> output: basin_indices(nbasins,npts).
    !             basins won't have contiguous memory...hmmmm

    integer :: cnt, i, kk
    integer, allocatable, dimension(:) :: tmp_indices
    
    allocate(tmp_indices(grid_var%npts))
    
    allocate(basins(nbasins))

    do i=1,nbasins
       tmp_indices(:) = 0
       cnt = 0
       do kk=1,grid_var%npts
          if (grid_var%ocean_outlet(kk) .eq. i) then
             cnt = cnt + 1
             tmp_indices(cnt) = kk
          end if
        end do
        
        ncells = cnt
        allocate(basins(i)%river_points(cnt))
        basins(i)%river_points(:) = tmp_indices(1:cnt)  !i like this solution.  simply pass basin indices to loop over. 
                                                        !will need to put these in contiguous array to pass back to master.
        
        !when doint the routing:
        !below works in the global space.  should work in local.  pass to master.  put back into global vector???
        !  unless I change how cable does mpi,  must gather all lsm results, send out for river, than gather back for next lsm step
        !
        !  river_water_by_basin(npts_largest_basin,nbasins)
        !  do i=1,nbasins
        !     do k=1,size(basins(i)%river_points(:))
        !        river_by_basin(k,i) = basins(i)%river_points(k)
        !     end do
        !   end do
        !   GIVES me 2d array.  if I sort basins.  then can send round robin style
        
        
        !do i=nbasins_pe_srt,nbasins_pe_end
        !   do k=1,size(basins(i)%river_points(:))
        !      kk = basins(i)%river_points(k)
        !      Fin(kk) = 0._r_2
        !   end do
        !
        !   do k=1,size(basins(i)%river_points(:))
        !      kk = basins(i)%river_points(k)  !current index
        !      !calc the output flux.  send to downstream cell
        !      Fout(kk) = river_var%vel * ...blah blah
        !      Fin(grid_var%dwnstrm_index(kk)) = Fin() + Fout(kk)
        !    end do
        !
        !   do k=1,size(basins(i)%river_points(:))
        !     kk = basins(i)%river_points(k) 
        !     river_var%water_mass(kk) = river_var%water_mass(kk)  + (Fin(kk) - Fout(kk))*dels
        !     if (grid_var%is_outlet .eq. 1) then
        !        total_outflow(i) = Fout(kk)*dels
        !     end if
        !   end do
        !end do

    end do

  end subroutine partition_basins_for_mpi


!-----------------------------------------------------------------------------  

  !Determine the wieghts and mapping from 1D global land grid to the 1D river grid
  !lat_in and lon_in are assumed higher resolution than lat_out lon_out
  
  !call determine_hilo_res_mappings(river_grid%lat(:),river_grid%lon(:),latitude(:),longitude(:),                      &
  !                                 river_grid%maps%weight_river_to_land(:,:),river_grid%maps%ind_river_to_land(:,:))

  subroutine determine_hilo_res_mappings(lat_hi,lon_hi,lat_lo,lon_lo,  &
                                         weight_hi_to_lo,map_hi_to_lo)

    real(r_2), dimension(:),     intent(in)  :: lat_lo,lon_lo  !the lo resolution grid
    real(r_2), dimension(:),     intent(in)  :: lat_hi,lon_hi  !the hi resolution grid

    integer,   dimension(:,:), intent(inout) :: map_hi_to_lo     !indices on the hi res grid for a given lo res point
    real(r_2), dimension(:,:), intent(inout) :: weight_hi_to_lo  !fraction of the hi res in grid the lo res

    integer   :: n_ovrlap(mp)   !the number of hi res grid cells at each lo res point that overlap to the lo res cell
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
    if (size(lat_hi) .gt. nlat_river .and. size(lon_hi) .gt. nlon_river) then
       dlat_hi  = lat_hi(nlon_river+1)-lat_hi(1)
       dlon_hi  = lon_hi(2) - lon_hi(1)
    else
       dlat_hi  = lat_hi(2) - lat_hi(1)
       dlon_hi  = lon_hi(2) - lon_hi(1)
    end if

    if (size(lat_lo) .gt. mlat .and. size(lon_lo) .gt. mlon) then
       dlat_lo  = lat_all(1,2)-lat_all(1,1)  !array is compacted,  must use lat_all lon_all brought in from io_vars module
       dlon_lo  = lon_all(2,1) - lon_all(1,1)
    else
       dlat_lo  = lat_lo(2) - lat_lo(1)
       dlon_lo  = lon_lo(2) - lon_lo(1)
    end if

    !init overlap mappings and areas
    n_ovrlap(:)          = 0
    weight_hi_to_lo(:,:) = 0._r_2
    map_hi_to_lo(:,:) = -1

    !map each input lat/lon to output grid and calc the fraction of cell within (1 cell can map to upto four output)
    do k=1,mp   !loop over all the land grid points
    
      Sedge_lo = lat_lo(k)  - dlat_lo/2._r_2
      Nedge_lo = lat_lo(k)  + dlat_lo/2._r_2    
      
      Wedge_lo = lon_lo(k)  - dlon_lo/2._r_2
      Eedge_lo = lon_lo(k)  + dlon_lo/2._r_2
    
      dy_lo = sin(deg2rad*Nedge_lo) - sin(deg2rad*Sedge_lo)
      area_lo = dy_lo * dlon_lo * re * re  !adjust for fraction of grid cells??
      
      do jj=1,ilat

        kk_tmp = 1 + jj*nlon_river

        Sedge_hi = lat_hi(kk_tmp)  - dlat_hi/2._r_2  !avoid some calcs dooing it here
        Nedge_hi = lat_hi(kk_tmp) + dlat_hi/2._r_2
        
        dy_hi = sin(deg2rad*Nedge_hi) - sin(deg2rad*Sedge_hi)
        area_hi = dy_hi * dlon_hi * re * re  !adjust for fraction of grid cells??

        if ((Sedge_hi .le. Nedge_lo) .and. (Nedge_hi .ge. Sedge_lo)) then  !overlap in lat

          do ii=1,ilon

            kk = ii + jj*nlon_river

            Wedge_hi = lon_hi(kk) - dlon_hi/2._r_2
            Eedge_hi = lon_hi(kk) + dlon_hi/2._r_2

            if ((Wedge_hi .le. Eedge_lo) .and. (Eedge_hi .ge. Wedge_lo)) then
              
              n_ovrlap(k) = n_ovrlap(k) + 1  !number of overlapping cells

              dlone = min(Eedge_lo,Eedge_hi)*deg2rad !determine area of input cell within the output grid cell
              dlonw = max(Wedge_lo,Wedge_hi)*deg2rad 
              dx = max(0.0,(dlone-dlonw))

              dlatn = min(Nedge_lo,Nedge_hi)*deg2rad 
              dlats = max(Sedge_lo,Sedge_hi)*deg2rad 
              dy = max(0.0,(sin(dlatn)-sin(dlats)))

              map_hi_to_lo(k,n_ovrlap(k))    = kk               !all river points for a given land point
              weight_hi_to_lo(k,n_ovrlap(k)) = dx*dy / area_hi  !frac of hi res cell within this lo res cell

            end if  !test lon overlap
            
          end do  !ii lon loop over input
          
        end if  !test lat overlap
        
      end do !jj loop over input lat
      
    end do  !loop over mp land points
    
    
    !determine runoff from lsm to river cells:
    !Q_river_runoff(:) = 0._r_2
    !do k=1,mp
    !  do i=1,n_overlap(k)
    !    Q_river_runoff(map_hi_to_lo(k,i)) = Q_river_runoff(map_hi_to_lo(k,i))  + Q_runoff(k)*weight_hi_to_lo(k,i)
    !  end do
    !end do
    
    !map the river water to each grid cell
    !wb_river(:) = 0._r_2
    !do k=1,mp
    !  do i=1,n_overlap(k)
    !    wb_river_land(k) = wb_river_land(k) + river%wb(map_hi_to_lo(k,i))*weight_hi_to_lo(k,i)
    !  end do
    !end do    


  end subroutine determine_hilo_res_mappings 

    
    
    
    
    
  








end module cable_routing
