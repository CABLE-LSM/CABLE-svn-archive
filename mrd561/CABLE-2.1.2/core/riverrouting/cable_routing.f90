module cable_routing

  use routing_constants, only :: rrlat,rrlon,rrlatlon
  use cable_types

  implicit none
  
  type map_grid_type
    private

    integer, pointer, dimension(:,:)   :: lat_overlap
    integer, pointer, dimension(:,:)   :: lon_overlap
    real(r_2), pointer, dimension(:,:) :: frac_overlap

  end type map_grid_type

  
  type river_grid_type    !will move to cable types eventually
    private
    
    integer                            :: nlat
    integer                            :: nlon
    integer                            :: npts
    
    real(r_2), pointer, dimension(:) :: lat   !1dimensional array stored as lat=0,lon={1,nlon}; lat=1,lon={1,nlon}, etc
    real(r_2), pointer, dimension(:) :: lon
    real(r_2), pointer, dimension(:) :: latS
    real(r_2), pointer, dimension(:) :: latN
    real(r_2), pointer, dimension(:) :: lonE
    real(r_2), pointer, dimension(:) :: lonW

    type(map_grid_type) :: maps
    
    
  end type river_grid_type
  
  
  type river_flow_type     !will need one for each mpi basin?
    private
    
    real(r_2), pointer, dimension(:) :: wat_vol
    real(r_2), pointer, dimension(:) :: wat_mass
    real(r_2), pointer, dimension(:) :: wat_vel
    
  end type river_flow_type
  
  integer, parameter    :: max_n_overlap = 2305   !assume at most 3x3 to 0.0625x0.0625

  real(r_2), parameter :: fNaN = -1e36
  integer, parameter :: iNaN = -999999

  type(river_grid_type) :: river_grid
  type(river_flow_type) :: river
    
    


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

  end subroutine alloc_river_flow

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------

  subroutine alloc_river_map_grid(var,npts)
    type(map_grid_type), intent(inout) :: var
    integer, intent(in)                :: npts

    allocate(var%lat_overlap(npts,max_n_overlap))
    var%lat_overlap(:,:) = iNaN
    allocate(var%lon_overlap(npts,max_n_overlap))
    var%lon_overlap(:,:) = iNaN
    allocate(var%frac_overlap(npts,max_n_overlap))
    var%frac_overlap(:,:) = fNaN

  end subroutine alloc_river_map_grid

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------

  subroutine alloc_river_grid(var,npts)

    type(river_grid_type), intent(inout) :: var
    integer, intent(in)             :: npts

    allocate(var%lat(npts))
    var%lat(:) = fNaN
    allocate(var%lon(npts))
    var%lon(:) = fNaN
    allocate(var%latS(npts))
    var%latS(:) = fNaN
    allocate(var%latN(npts))
    var%latN(:) = fNaN    
    allocate(var%lonE(npts))
    var%lonE(:) = fNaN    
    allocate(var%lonW(npts))
    var%lonW(:) = fNaN    

    call alloc_river_map_grid(var%maps,npts)


  end subroutine alloc_river_grid

!-----------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------
!             !current data set doesn't use these numbers as it is: 1,2,4,8,16,32,64,128,256
!                8   1   2
!                7   0   3
!                6   5   4
!
!

  integer function dirc2latindex(in_dirc) 
    implicit none
    integer, intent(in) :: in_dirc

    dirc2latindex = 0
    if ((in_dirc .le. 1) .or.  (in_dirc .eq. 9)) dirc2latindex  = 1
    if ((in_dirc .ge. 3) .and. (in_dirc .le. 5)) dirc2latindex = -1

  end function dirc2latindex
!-----------------------------------------------------------------------------
  

!-----------------------------------------------------------------------------
!
  integer function dirc2lonindex(in_dirc) 
    implicit none
    integer, intent(in) :: in_dirc

    dirc2lonindex = 0
    if ((in_dirc .ge. 1) .and. (in_dirc .le. 3)) dirc2lonindex  = 1
    if ((in_dirc .ge. 5) .and. (in_dirc .le. 7)) dirc2lonindex = -1


  end function dirc2lonindex
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!
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
  
  subroutine associate_ocean_outlet_points(ocean_outlet,upstrm_number,nbasins,nrr_cells,land_mask)
  
    integer,   dimension(:), intent(out) :: ocean_outlet
    integer,   dimension(:), intent(out) :: upstrm_number
    integer,                 intent(out) :: nbasins
    integer,                 intent(out) :: nrr_cells
    integer,   dimension(:), intent(in)  :: land_mask
  
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
    
    
    !need to allocate basins to nprocs for mpi
    !can do it this way....
    
  subroutine allocate_local_basins(river_var,grid_var,nbasins)
    type(river_flow_type), dimension(:), intent(out) :: river_var
    type(river_grid_type), dimension(:), intent(out) :: grid_var
    integer,                             intent(in)  :: nbasins
    
    
    allocate(river_var(nbasins))
    allocate(grid_var(nbasins))
    
    !call other allocs here??
    
  end subroutine allocate_local_basins
  
  
    !or I can simply pass around parts of the large array
    
    subroutine partition_basins_for_mpi(river_var,grid_var,nbasins,upstrm_number,ocean_outlet)
    
    !loop: basins
    !loop : find all for basin i
    !  create vector of indices for basin i --> output: basin_indices(nbasins,npts).
    !             basins won't have contiguous memory...hmmmm
    
    
    
    
    
  








end module cable_routing
