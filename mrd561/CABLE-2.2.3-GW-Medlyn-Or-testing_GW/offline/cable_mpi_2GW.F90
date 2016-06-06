MODULE cable_mpi_2dgw

  USE cable_mpicommon
  USE cable_def_types_mod
  USE cable_common_module
  USE cable_IO_vars_module, ONLY : landpt, patch,mask,land_x,land_y
  USE cable_2dgw_types


  IMPLICIT NONE

  SAVE


contains


subroutine alloc_halo_param_type(var,npts)
  implicit none
  type(gw_halo_param_type), intent(out) :: var
  integer,            intent(in)  :: npts

  var%ncells = npts

  allocate(var%hksat(npts))
  var%hksat(:) = 0._r_2

  
  allocate(var%slope(npts))
  var%slope(:) = 0._r_2

  allocate(var%elv(npts))
  var%elv(:) = 0._r_2

  allocate(var%latitude(npts))
  var%latitude(:) = 0._r_2

END SUBROUTINE alloc_halo_param_type

subroutine alloc_halo_var_type(var,npts)
  implicit none
  type(gw_halo_var_type), intent(out) :: var
  integer,            intent(in)  :: npts

  var%ncells = npts

  allocate(var%wtd(npts))
  var%wtd(:) = 0._r_2

  allocate(var%tgg_ms(npts))
  var%tgg_ms(:) = 0._r_2

END SUBROUTINE alloc_halo_var_type


!below is a subroutine to break up the 1d array to the workers
!ensureing that only whole rows of grid points are passed
!this makes the mpi passing for a 2d gw scheme much simpler
!likely hurts the load balancing for the last point though
SUBROUTINE master_decomp_2dgw(comm,mland,mp,wland)

  USE mpi
  IMPLICIT NONE


  INTEGER, INTENT(IN)   :: comm ! MPI communicator to talk to the workers
  INTEGER, INTENT(IN)   :: mland ! total number of landpoints in the global grid
  INTEGER, INTENT(IN)   :: mp ! total number of land patches in the global grid
  TYPE(lpdecomp_t), ALLOCATABLE, DIMENSION(:), intent(out) :: wland

  INTEGER :: lpw  ! average number of landpoints per worker
  INTEGER :: rank, rest, nxt, pcnt, ierr, i, tmp, wnp, iworker
  INTEGER :: patchcnt  ! sum of patches for a range landpoints

  real(r_2) :: total_points
  real(r_2) :: total_land_points
  real(r_2) :: ideal_fraction_land_per_worker
  real(r_2) :: current_worker_fraction

  integer :: j,k,i

  integer :: current_lat_index
  logical :: keep_looking
  integer :: current_mland_index
  integer :: j_index,i_index
  integer :: j_start,j_end
  integer :: send_south_halo_npoints,send_north_halo_npoints
  integer :: recv_south_halo_npoints,recv_north_halo_npoints
  integer, allocatable, dimension(:,:) :: map_index_worker
  integer :: worker_mlat

  ! how many workers do we have?
  CALL MPI_Comm_size (comm, wnp, ierr)
  wnp = wnp - 1

  ALLOCATE (wland(wnp), STAT=ierr)
  ALLOCATE (master_halo(wnp), STAT=ierr)

  IF (ierr /= 0) THEN
          ! TODO: print an error message
          CALL MPI_Abort(comm, 0, ierr)
  END IF

  total_points = real(mlon,r_2)*real(mlat,r_2)
  total_land_points = real(mland,r_2)
  ideal_fraction_land_per_worker = total_points/wnp

  current_mland_index = 0
  nxt = 1

  do iworker = 1,wnp

     current_mland_index = current_mland_index + 1

     current_worker_fraction = 0._r_2

     wland(rank)%landp0 = nxt   !first point for this worker in the global array

     j_start = land_y(current_mland_index)

     keep_looking = .true.
     do while (keep_looking) 


        j_index = land_y(current_mland_index)
        i_index = land_x(current_mland_index)

        current_worker_fraction = current_worker_fraction + 1._r_2

        if (((i_index .eq. mlon) .or. (sum(mask(i_index:mlon,j_index)) .eq. 0))  .and. &  !end of row or no more land points in row
           (current_worker_fraction .ge. 0.95*ideal_fraction_land_per_worker)) then            !have accumulated enough land points
              keep_looking = .false.
        else
           current_mland_index = current_mland_index + 1
        end if

     end do

     j_index = land_y(current_mland_index)
     i_index = land_x(current_mland_index)

     j_end = j_index

     !for the halo find the number of land points at j_start -1 and j_end + 1
     if (j_start .gt. 1) then
        recv_north_halo_npoints = sum(mask(:,j_start-1),dim=1)
        send_north_halo_npoints = sum(mask(:,j_start),dim=1)
     else
        recv_north_halo_npoints = 0
        send_north_halo_npoints = 0
     end if
     if (j_end .lt. mlat) then
        recv_south_halo_npoints = sum(mask(:,j_end+1),dim=1)
        send_south_halo_npoints = sum(mask(:,j_end),dim=1)
     else
        recv_south_halo_npoints = 0
        send_south_halo_npoints = 0
     end if

     worker_mlat = j_end - j_start + 1
     allocate(map_index_worker(mlon,0:worker_mlat+1),STAT=ierr)

     map_index_worker(:,:) = map_index(:,j_start:j_end) - nxt + 1  !adjust for all previous workers

     k=0
     do i=1,mlon
        if (mask(i,j_start-1) .eq. 1) then
           k = k + 1
           map_index_worker(i,0) = k
        else
           map_index_worker(i,0) = -1
        end if
     end do

     k=0
     do i=1,mlon
        if (mask(i,j_end+1) .eq. 1) then
           k = k + 1
           map_index_worker(i,worker_mlat+1) = k
        else
           map_index_worker(i,worker_mlat+1) = -1
        end if
     end do
         
     wland(iworker)%nland = current_mland_index
     nxt = nxt + current_mland_index

     ! MPI: let each worker know their assignement
     ! in this version of cable workers care only about the number of points
     ! CALL MPI_Send (nxt, 1, MPI_INTEGER, rank, 0, comm, ierr)
     CALL MPI_Send (current_mland_index, 1, MPI_INTEGER, iworker, 0, comm, ierr)

     ! MPI: should be the same as landpt(nxt)%cstart
     wland(iworker)%patch0 = landpt(nxt)%cstart
     ! MPI: calculate no of patches for pcnt landpoint starting from nxt
     ! MPI: TODO: workers only need to know the number of their patches
     ! or maybe not (if they need patch displacement in the input)

     ! MPI: find number of patches in all landpoints assigned to this
     ! worker (difference between last patch of last landpoint and first of
     ! first)
     patchcnt = landpt(nxt+current_mland_index-1)%cend - landpt(current_mland_index)%cstart + 1
     wland(iworker)%npatch = patchcnt

     ! MPI: development check
     tmp = 0
     DO i = 1, current_mland_index
        tmp = tmp + landpt(nxt+i-1)%nap
     END DO
     IF (patchcnt /= tmp) THEN
        WRITE (*,*) 'invalid patch number for worker ', &
        &           patchcnt,tmp,iworker
        CALL MPI_Abort (comm, 0, ierr)
     END IF

     master_halo(iworker)%npts_nothern = recv_north_halo_npoints
     master_halo(iworker)%npts_southern = recv_south_halo_npoints

     ! MPI: at this point workers can't determine patches on their own
     ! so we need to send it explicitely
     ! in this version of cable workers care only about the number of patches
     ! CALL MPI_Send (wland(rank)%patch0, 1, MPI_INTEGER, rank, 0, comm, ierr)
     CALL MPI_Send (wland(iworker)%npatch, 1, MPI_INTEGER, iworker, 0, comm, ierr)

     !send the number of send land points in the halo to the north of this worker
     CALL MPI_Send (send_north_halo_npoints, 1, MPI_INTEGER, iworker, 0, comm, ierr)
     !send the number of send land points in the halo to the south of this worker
     CALL MPI_Send (send_south_halo_npoints, 1, MPI_INTEGER, iworker, 0, comm, ierr)

     !send the number of recv land points in the halo to the north of this worker
     CALL MPI_Send (recv_north_halo_npoints, 1, MPI_INTEGER, iworker, 0, comm, ierr)
     !send the number of recv land points in the halo to the south of this worker
     CALL MPI_Send (recv_south_halo_npoints, 1, MPI_INTEGER, iworker, 0, comm, ierr)

     !send number of longitude points
     CALL MPI_Send (mlon, 1, MPI_INTEGER, iworker, 0, comm, ierr)
     CALL MPI_Send (worker_mlat, 1, MPI_INTEGER, iworker, 0, comm, ierr)

     !send the entire 2d mask to each worker
     CALL MPI_Send (map_index_worker(1,0), mlon*(worker_mlat+2), MPI_INTEGER, iworker, 0, comm, ierr)
     !also need to send the map_index array
     !this needs to me calculated for each worker
     !check this
     !nxt = nxt + pcnt
     !not dealing with tiling yet

     deallocate(map_index_local)

end do   !loop over the workers     

END SUBROUTINE master_decomp_2dgw


! MPI: receives grid decomposition info from the master
SUBROUTINE worker_TwoDGW_halo_sizes(comm)

  USE mpi

  USE cable_def_types_mod, ONLY: mland, mp,mlon


  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: comm ! MPI communicator to talk to the workers

  INTEGER :: stat(MPI_STATUS_SIZE), ierr
  INTEGER :: worker_mlat

  !get number of points to send to the north (1) and south (2)
  CALL MPI_Recv (worker_dims%npts_send(1), 1, MPI_INTEGER, 0, 0, comm, stat, ierr)
  CALL MPI_Recv (worker_dims%npts_send(2), 1, MPI_INTEGER, 0, 0, comm, stat, ierr)
  !get the number of points to recv from north (1) and south (2)
  CALL MPI_Recv (worker_dims%npts_recv(1), 1, MPI_INTEGER, 0, 0, comm, stat, ierr)
  CALL MPI_Recv (worker_dims%npts_recv(2), 1, MPI_INTEGER, 0, 0, comm, stat, ierr)
  !nubmer of lon points (same for each worker)
  CALL MPI_Recv (mlon, 1, MPI_INTEGER, 0, 0, comm, stat, ierr)
  CALL MPI_Recv (mlat, 1, MPI_INTEGER, 0, 0, comm, stat, ierr)
  !since decomp does all longitudes per proc just set the i_start and i_end
  !this var needs to be allocated first
  allocate(worker_dims%worker_map_index(mlon,0:mlat+1))
  worker_dims%worker_map_index(:,:) = -1
  CALL MPI_Recv (worker_dims%worker_map_index(1,0),(mlat+2)*mlon,MPI_INTEGER,0,0,comm,stat,ierr)
  !note row 0 is for data from northern worker
  ! row w_mlat+1 is from southern worker

  worker_dims%i_start = 1
  worker_dims%i_end   = mlon


  RETURN

END SUBROUTINE worker_TwoDGW_halo_sizes

SUBROUTINE master_send_2dgw_parameters(comm,soil,wland)


  implicit none

  INTEGER, INTENT(IN) :: comm ! MPI communicator
  TYPE (soil_parameter_type), INTENT(OUT)    :: soil
  TYPE(lpdecomp_t), DIMENSION(:), intent(in) :: wland

  !local
  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen 
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: stat(MPI_STATUS_SIZE), ierr 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: param_ts

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride
  integer :: bidx
  INTEGER :: tsize, localtotal, remotetotal
  INTEGER :: i, wnp

  wnp = size(wland(:),dim=1)

  ALLOCATE (param_ts(wnp))

  ALLOCATE (blen(ngw_worker_recv_param_types))
  ALLOCATE (displs(ngw_worker_recv_param_types))
  ALLOCATE (types(ngw_worker_recv_param_types))


  ! total size of input data sent to all workers
  localtotal = 0

  do rank=1,wnp

     !send the parameters
     !first to the northern halo
     !then to the southern halo
     ! starting patch and number for each worker rank - number to get from north
     off = wland(rank)%patch0 - master_halo(rank)%npts_nothern   !check this

     bidx = 0

     !create the nothern halo send type
     bidx = bidx + 1
     CALL MPI_Get_address (soil%hksat(off,ms), displs(bidx), ierr)
     blen(bidx) = extr2*master_halo%npts_nothern

     bidx = bidx + 1
     CALL MPI_Get_address (soil%slope(off), displs(bidx), ierr)
     blen(bidx) = extr2*master_halo%npts_nothern

     bidx = bidx + 1
     CALL MPI_Get_address (soil%elev(off), displs(bidx), ierr)
     blen(bidx) = extr2*master_halo%npts_nothern

     bidx = bidx + 1
     CALL MPI_Get_address (latitude(off), displs(bidx), ierr)
     blen(bidx) = extr1*master_halo%npts_nothern

     !now do the southern halo type
     off = wland(rank)%patch0 + wland(rank)%npatch + master_halo(rank)%npts_southern + 1  !check this

     bidx = bidx + 1
     CALL MPI_Get_address (soil%hksat(off,ms), displs(bidx), ierr)
     blen(bidx) = extr2*master_halo%npts_nothern

     bidx = bidx + 1
     CALL MPI_Get_address (soil%slope(off), displs(bidx), ierr)
     blen(bidx) = extr2*master_halo%npts_nothern

     bidx = bidx + 1
     CALL MPI_Get_address (soil%elev(off), displs(bidx), ierr)
     blen(bidx) = extr2*master_halo%npts_nothern

     bidx = bidx + 1
     CALL MPI_Get_address (latitude(off), displs(bidx), ierr)
     blen(bidx) = extr1*master_halo%npts_southern


     ! MPI: sanity check
     IF (bidx /= ngw_worker_recv_param_types) THEN
        WRITE (*,*) 'master: invalid number of param_t fields ',bidx,', fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blen, displs, types, param_ts(rank), ierr)
     CALL MPI_Type_commit (param_ts(rank), ierr)

     CALL MPI_Type_size (param_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (param_ts(rank), tmplb, text, ierr)

     write (*,*) 'master to rank param_t blocks, size, extent and lb: ',bidx,tsize,text,tmplb

     localtotal = localtotal + tsize

  end do ! rank

  WRITE (*,*) 'total cable params size sent to all workers: ', localtotal

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blen)

  ! MPI: check whether total size of send input data equals total
  ! data received by all the workers
  remotetotal = 0
  CALL MPI_Reduce (MPI_IN_PLACE, remotetotal, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  WRITE (*,*) 'total cable params size received by all workers: ', remotetotal

  IF (localtotal /= remotetotal) THEN
          WRITE (*,*) 'error: total length of cable params sent and received differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  ! so, now send all the parameters
  CALL twodim_master_send_input (comm, param_ts, 0)

  CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  ! finally free the MPI type
  DO rank = 1, wnp
     CALL MPI_Type_Free (param_ts(rank), ierr)
  END DO

END SUBROUTINE master_send_2dgw_parameters


SUBROUTINE worker_get_2dgw_parameters(northern_halo_params,southern_halo_params)
  implicit none

  type(gw_halo_param_type), intent(out) :: northern_halo_params
  type(gw_halo_param_type), intent(out) :: southern_halo_params

  ! temp arrays for marshalling all fields into a single struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen 
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  ! temp vars for verifying block number and total length of inp_t
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb
  INTEGER :: tsize, localtotal, remotetotal

  INTEGER :: stat(MPI_STATUS_SIZE), ierr 
  INTEGER :: param_t

  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride
  integer :: bidx
  INTEGER :: tsize, localtotal, remotetotal
  INTEGER :: i, wnp

  ntyp = ngw_worker_recv_param_types

  !allocate (northern_halo_params,worker_dims%npts_recv(1))
  !allocate (southern_halo_params,worker_dims%npts_recv(2))

  call alloc_halo_param_type(northern_halo_params,worker_dims%npts_recv(1))
  call alloc_halo_param_type(southern_halo_params,worker_dims%npts_recv(2))

  ALLOCATE (blen(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  bidx = 0

  ! the order of variables follows argument list
  ! the order of fields within follows alloc_*_type subroutines

  ! ----------- met --------------

  bidx = bidx + 1
  CALL MPI_Get_address (northern_halo_params%hksat(:), displs(bidx), ierr)
  blen(bidx) = extr2*northern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (northern_halo_params%slope(:), displs(bidx), ierr)
  blen(bidx) = extr2*northern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (northern_halo_params%elev(:), displs(bidx), ierr)
  blen(bidx) = extr2*northern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (northern_halo_params%latitude(:), displs(bidx), ierr)
  blen(bidx) = extr1*northern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (southern_halo_params%hksat(:), displs(bidx), ierr)
  blen(bidx) = extr2*southern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (southern_halo_params%slope(:), displs(bidx), ierr)
  blen(bidx) = extr2*southern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (southern_halo_params%elev(:), displs(bidx), ierr)
  blen(bidx) = extr2*southern_halo_params%ncells

  bidx = bidx + 1
  CALL MPI_Get_address (southern_halo_params%latitude(:), displs(bidx), ierr)
  blen(bidx) = extr1*northern_halo_params%ncells

  ! MPI: sanity check
  IF (bidx /= ntyp) THEN 
     WRITE (*,*) 'worker ',rank,' invalid number of param_t fields',bidx,', fix it!'
     CALL MPI_Abort (comm, 1, ierr)
  END IF

  CALL MPI_Type_create_struct (bidx, blen, displs, types, param_t, ierr)
  CALL MPI_Type_commit (param_t, ierr)

  CALL MPI_Type_size (param_t, tsize, ierr)
  CALL MPI_Type_get_extent (param_t, tmplb, text, ierr)

  WRITE (*,*) 'worker param_t blocks, size, extent and lb: ',rank,bidx,tsize,text,tmplb

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  CALL MPI_Reduce (tsize, MPI_DATATYPE_NULL, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blen)

  ! if anything went wrong the master will mpi_abort
  ! which mpi_recv below is going to catch...

  ! so, now receive all the parameters
  CALL MPI_Recv (MPI_BOTTOM, 1, param_t, 0, 0, comm, stat, ierr)

  ! finally free the MPI type
  CALL MPI_Type_Free (param_t, ierr)

  ! all CABLE parameters have been received from the master by now
  RETURN


END SUBROUTINE worker_get_2dgw_parameters


!could use routine in cable_master but that yields circular dependency hell
SUBROUTINE twodim_master_send_input (comm, dtypes, ktau)

  USE mpi

  IMPLICIT NONE 

  INTEGER, INTENT(IN) :: comm 
  INTEGER, DIMENSION(:), INTENT(IN) :: dtypes
  INTEGER, INTENT(IN) :: ktau    ! timestep

  INTEGER :: rank, ierr 

!  IF (.NOT. ALLOCATED(inp_req)) THEN
!     ALLOCATE (inp_req(wnp))
!  END IF

  DO rank = 1, wnp
     CALL MPI_Isend (MPI_BOTTOM, 1, dtypes(rank), rank, ktau, comm, &
     &               inp_req(rank), ierr)
  END DO

  !IF (.NOT. ALLOCATED(inp_stats)) THEN
  !   ALLOCATE (inp_stats(MPI_STATUS_SIZE, wnp))
  !END IF

  !CALL MPI_Waitall (wnp, inp_req, inp_stats, ierr)

  RETURN

END SUBROUTINE twodim_master_send_input


SUBROUTINE worker_pass_halos(comm,wnp,ssnow,northern_halo_var,southern_halo_var)
   use mpi

   implicit none

   integer, intent(in)  :: comm
   integer, intent(in)  :: wnp  !total number of worker processes
   type(soil_snow_type), intent(in) :: ssnow
   type(gw_halo_var_type), intent(inout) :: northern_halo_var, &
                                            southern_halo_var

   integer :: ierr
   integer :: statuses(MPI_STATUS_SIZE,8), requests(8) 
   integer :: ib,ie
   integer :: i_rq,j_rq,n_rq

   if (rank .gt. 2) then   !in current mpi version rank == 1 does only IO
      call MPI_Irecv ( northern_halo_var%wtd(:), northern_halo_var%ncells, MPI_DOUBLE_PRECISION ,rank-1,ring_comm, requests(2), ierr)
      call MPI_Irecv ( northern_halo_var%tgg_ms(:), northern_halo_var%ncells, MPI_SINGLE_PRECISION, rank-1,ring_comm, requests(4), ierr)
   end if

   if (rank .lt. wnp) then
      call MPI_Irecv ( southern_halo_var%wtd(:), southern_halo_var%ncells, MPI_DOUBLE_PRECISION ,rank+1,ring_comm, requests(6), ierr)
      call MPI_Irecv ( southern_halo_var%tgg_ms(:), southern_halo_var%ncells, MPI_SINGLE_PRECISION, rank+1,ring_comm, requests(8), ierr)
   end if

   ib = 1
   ie = ib - worker_dims(rank)%npts_send(1)

   if (rank .gt. 2) then
   
      call MPI_Isend( ssnow%wtd(ib:ie), m, MPI_DOUBLE_PRECISION , rank-1, ring_comm, requests(1), ierr)
      call MPI_Isend( ssnow%tgg(ib:ie,ms), m, MPI_SINLGE_PRECISION , rank-1, ring_comm, requests(3), ierr)
   end if

   ie = size(ssnow%wtd,dim=1)
   ib = ie - worker_dims(rank)%npts_send(2)

   if (rank .lt. wnp) then
      call MPI_Isend( ssnow%wtd(ib:ie), m, MPI_DOUBLE_PRECISION , rank+1, ring_comm, requests(5), ierr)
      call MPI_Isend( ssnow%tgg(ib:ie,ms), m, MPI_SINLGE_PRECISION , rank+1, ring_comm, requests(7), ierr)
   end if

   if (rank .gt. 2 .and. rank .lt. wnp) then
      i_rq = 1
      j_rq = 8
      n_rq = 8
   elseif (rank .gt. 2) then
      i_rq = 1
      j_rq = 4
      n_rq = 4
   else
      i_rq = 5
      j_rq = 8
      n_rq = 4
   end if

   call MPI_Waitall( n_rq, requests(i_rq:j_rq), statuses(:,i_rq:j_rq),ierr) 

END MODULE cable_mpi_2dgw
