MODULE mpi_2dgw_types

  USE mpi
  USE cable_mpicommon
  USE cable_def_types_mod

  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER :: ngw_types      = 6  !number of worker types
  INTEGER, PARAMETER :: ngw_halo_types = 1  !number of worker types
  INTEGER, PARAMETER :: ngw_master_recv_types = 2   !data master recv from workers
  INTEGER, PARAMETER :: ngw_master_send_types = 4   !data master sends to workers
  INTEGER, PARAMETER :: ngw_worker_send_types = 2
  INTEGER, PARAMETER :: ngw_worker_recv_types = 4 

  INTEGER, ALLOCATABLE, DIMENSION(:) :: worker_nx   !nx varies,  all do ny
  INTEGER, ALLOCATABLE, DIMENSION(:) :: worker_nx_nohalo !number of worker
                                                         !columns not including the halos

  INTEGER, ALLOCATABLE, DIMENSION(:) :: master_send_2dgw
  INTEGER, ALLOCATABLE, DIMENSION(:) :: master_recv_2dgw

  INTEGER                            :: worker_send_2dgw
  INTEGER                            :: worker_recv_2dgw
  INTEGER                            :: send_to_right,send_to_left
  INTEGER                            :: recv_from_right,recv_from_left

  INTEGER, ALLOCATABLE, DIMENSION(:) :: mrecv_req, msend_req

  INTEGER                             :: stat(MPI_STATUS_SIZE)

!code outline
! define twodim_gw_type_module
!  declare global and local arrays that are saved.  other wise can't
!  have constant derived mpi types.  addresses will change
!
!  prior to timestepping:
!  if (myrank .eq. 0) then
!     master_send_gw_types (comm,h,elv,hycond,poros)  !global variables
!     master_recv_gw_types (comm,h,gwconv)
!  else
!     worker_send_gw_types (comm,myrank,h,gwconv)  !local variables
!     worker_recv_gw_types (comm,myrank,h,elv,hycond,poros_local)  !all are
!                                                            local variables
!     worker_2_worker_halo_gw_type(comm,myrank,h_local)
!  end if
!
!  start time stepping so withing actual calc 2d subroutine call 

!  if (myrank .eq. 0) then
!     master_send_data(comm,master_send_2dgw,ktau)
!  else
!     worker_recv_data(comm,myrank,worker_recv_2dgw,ktau)
!  end if
!
!  if (myrank .ne. 0)
!  begin the adi loop to solve for convgw and hn
!     if (iter .gt. 1) then
!        comm,myrank,iter,send_to_left,send_to_right,&
!                                                 recv_from_left,recv_from_right)
!     end if
!  end adi loop
!  end if
!
!  if (myrank .eq. 0)
!    master_recv_data()
!    mpi_wait_all
!  else
!    worker_send_data()
!  end if

! done with timestep



contains

SUBROUTINE decompose_global_to_workers(comm,worker_nx,worker_nx_nohalo)
  implicit none
  integer, intent(in) :: comm
  integer, intent(in) :: nworkers
  integer, dimension(:), intent(out) :: worker_nx
  integer, dimension(:), intent(out) :: worker_nx_nohalo

  integer :: myrank
  integer :: npts
  integer :: npts_leftover
  integer :: i,irank

  CALL MPI_Comm_rank (comm, myrank, ierr)

  npts = mlon / nworkers
  npts_leftover = mlon - npts*nworkers

  do i=1,nworkers
     npts_leftover = npts_leftover - 1
     if (npts_leftover .ge. 0) then
       worker_nx_nohalo(i) = npts + 1
     else
       worker_nx_nohalo(i) = npts
     end if
     
     if ((i .eq. 1) .or. (i .eq. nworkers)) then
       worker_nx(i) = worker_nx_nohalo(i) + 1
     else
       worker_nx(i) = worker_nx_nohalo(i) + 2
     end if

  end do
  

END SUBROUTINE decompose_global_to_workers

SUBROUTINE worker_2_worker_halo_gw_type(comm,myrank,h)

  INTEGER, INTENT(IN) :: comm
  INTEGER, INTENT(IN) :: myrank
  REAL(r_2), DIMENSION(:,:), INTENT(IN) :: h

  INTEGER, ALLOCATABLE, DIMENSION(:)   :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  INTEGER :: off, cnt
  INTEGER :: bidx, midx, vidx, ierr
  INTEGER :: tsize
  INTEGER :: i



  ALLOCATE (blocks(ngw_halo_types))
  ALLOCATE (displs(ngw_halo_types))
  ALLOCATE (types(ngw_halo_types))

  r2len = cnt * extr2

  !create send_to_right to send 2nd to last column to the right
  off = worker_nx(myrank) - 1     !off set should be worker dependant for what master sends, 1 for what workers send back
  cnt = 1     !processing 1 column at a time
  bidx = 0
  bidx = bidx + 1
  CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat

  types = MPI_BYTE
  CALL MPI_Type_create_struct (bidx, blocks, displs, types, send_to_right(myrank),ierr)
  CALL MPI_Type_commit (send_to_right(myrank), ierr)


  !now create the send to the left
  !note 1st worker won't send anything but create type anyways?
  off = 2     !2nd from left is sent to the left
  cnt = 1     !processing 1 column at a time
  bidx = 0
  bidx = bidx + 1
  CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat
  CALL MPI_Type_create_struct (bidx, blocks, displs, types, send_to_left(myrank),ierr)
  CALL MPI_Type_commit (send_to_left(myrank), ierr)


  !create the recv from the right
  off = worker_nx(myrank)     !last colunm in the local array
  cnt = 1     !processing 1 column at a time
  bidx = 0
  bidx = bidx + 1
  CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat
  types = MPI_BYTE
  CALL MPI_Type_create_struct (bidx, blocks, displs, types,recv_from_right(myrank),ierr)
  CALL MPI_Type_commit (recv_from_right(myrank), ierr)

  !create recieve from the left
  off = 1   !first column of the local array
  cnt = 1     !processing 1 column at a time
  bidx = 0
  bidx = bidx + 1
  CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat
  types = MPI_BYTE
  CALL MPI_Type_create_struct (bidx, blocks, displs,types,recv_from_left(myrank),ierr)
  CALL MPI_Type_commit (recv_from_left(myrank), ierr)


  deallocate(blocks)
  deallocate(bidx)
  deallocate(displs)

END SUBROUTINE worker_2_worker_halo_gw_type

SUBROUTINE worker_recv_gw_types (comm,myrank,h,elv,hycond,poros)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: comm ! MPI communicator to talk to the workers
  INTEGER, INTENT(IN) :: myrank
  !these are the local worker arrays
  REAL(r_2), DIMENSION(:,:), INTENT(IN) :: h,       &
                                           elv,     &
                                           hycond,  &
                                           poros


  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  INTEGER :: off, cnt
  INTEGER :: bidx, midx, vidx, ierr
  INTEGER :: tsize
  INTEGER :: i


  ! MPI: calculate the sizes/extents of Fortran types used by
  ! CABLE
  ! do I need custom find extents??  No
  !ALL find_extents

  ! MPI: allocate temp vectors used for marshalling
  ALLOCATE (blocks(ngw_worker_recv_types))
  ALLOCATE (displs(ngw_worker_recv_types))
  ALLOCATE (types(ngw_worker_recv_types))

  off = 1     !off set should be worker dependant for what master sends, 1 for what workers send back
  ! so master will make a type for each worker starting at
  ! worker_nx(worker_rank)
  cnt = worker_nx(myrank)  !npts not set of declared  !includes left and right halo columns

  r2len = cnt * extr2

  bidx = 0

  bidx = bidx + 1
  CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat

  bidx = bidx + 1
  CALL MPI_Get_address (evl(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat

  bidx = bidx + 1
  CALL MPI_Get_address (poros(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat

  bidx = bidx + 1
  CALL MPI_Get_address (hycond(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat


  ! MPI: sanity check
  IF (bidx .ne. ngw_types) THEN
     WRITE (*,*) 'worker ',rank,': invalid outtype nmat, nvec or n3d constant, fix it!'
     CALL MPI_Abort (comm, 1, ierr)
  END IF

  types = MPI_BYTE
  CALL MPI_Type_create_struct (bidx, blocks, displs, types, worker_recv_2dgw, ierr)
  CALL MPI_Type_commit (worker_recv_2dgw, ierr)


  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)


END SUBROUTINE worker_recv_gw_types


SUBROUTINE worker_send_gw_types (comm,myrank,h,gwconv)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: comm ! MPI communicator to talk to the workers
  INTEGER, INTENT(IN) :: myrank
  !these are the local worker arrays
  REAL(r_2), DIMENSION(:,:), INTENT(IN) :: h,       &
                                           gwconv


  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  INTEGER :: off, cnt
  INTEGER :: bidx, midx, vidx, ierr
  INTEGER :: i


  ! MPI: calculate the sizes/extents of Fortran types used by
  ! CABLE
  ! do I need custom find extents??  No
  !ALL find_extents

  ! MPI: allocate temp vectors used for marshalling
  ALLOCATE (blocks(ngw_worker_send_types))
  ALLOCATE (displs(ngw_worker_send_types))
  ALLOCATE (types(ngw_worker_send_types))

  if (myrank .ne. 1) then  !first local array send 1-end-1, last its 2->end, else sending 2->end-1
     off = 2     !off set should be worker dependant for what master sends, 1 for what workers send back
  else
     off = 1
  end if
  ! so master will make a type for each worker starting at
  ! worker_nx(worker_rank)
  cnt = worker_nx_nohalo(myrank)  !npts not set of declared  !includes left and right halo columns

  r2len = cnt * extr2

  bidx = 0

  bidx = bidx + 1
  CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat

  bidx = bidx + 1
  CALL MPI_Get_address (gwconv(off,1), displs(bidx), ierr)
  blocks(bidx) = r2len * mlat


  ! MPI: sanity check
  IF (bidx .ne. ngw_worker_send_types) THEN
     WRITE (*,*) 'worker ',rank,': invalid outtype nmat, nvec or n3d constant, fix it!'
     CALL MPI_Abort (comm, 1, ierr)
  END IF

  types = MPI_BYTE
  CALL MPI_Type_create_struct (bidx, blocks, displs, types, worker_send_2dgw, ierr)
  CALL MPI_Type_commit (worker_send_2dgw, ierr)

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)


END SUBROUTINE worker_send_gw_types


SUBROUTINE master_send_gw_types (comm,h,elv,hycond,poros)
! must be called with if (myrank .eq. 0) master_send_gw_types(....)

  USE mpi
  USE cable_mpicommon

  USE cable_def_types_mod

  IMPLICIT NONE
  INTEGER :: comm ! MPI communicator to talk to the workers
  !these are the global arrays
  REAL(r_2), DIMENSION(:,:), INTENT(IN) :: h,       &
                                           elv,     &
                                           hycond,  &
                                           poros


  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types

  INTEGER :: off, cnt
  INTEGER :: bidx, ierr
  INTEGER :: i


  ! MPI: calculate the sizes/extents of Fortran types used by
  ! CABLE
  ! do I need custom find extents??  No
  !ALL find_extents

  ! MPI: allocate temp vectors used for marshalling
  ALLOCATE (blocks(ngw_types))
  ALLOCATE (displs(ngw_types))
  ALLOCATE (types(ngw_types))

  ! MPI: should work because worker versions of CABLE variables
  ! are allocated with starting index set to patch0
  !off = wpatch%patch0
  !cnt = wpatch%npatch
  ! MPI: new version, all arrays indices run 1:mp


  do irank=1,nworkers   !need to define nworkers (nprocs -1 )

     if (irank .eq. 1) then  !offset in the global array
        off = 1
     else
        !need to correct for the halos
        off = sum(worker_nx_nohalo(1:irank-1))
     end if

     ! so master will make a type for each worker starting at
     ! worker_nx(worker_rank)
     cnt = worker_nx(irank)  !includes left and right halos  data without halo
                             !is send back from worker though
     r2len = cnt * extr2
 
     bidx = 0

  ! MPI: an hvector type for each vector, maddr contains displacements
     bidx = bidx + 1
     CALL MPI_Get_address (hs(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     bidx = bidx + 1
     CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     bidx = bidx + 1
     CALL MPI_Get_address (sf2(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     bidx = bidx + 1
     CALL MPI_Get_address (evl(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     bidx = bidx + 1
     CALL MPI_Get_address (poros(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     bidx = bidx + 1
     CALL MPI_Get_address (hycond(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat


     ! MPI: sanity check
     IF (bidx .ne. ngw_types) THEN
        WRITE (*,*) 'worker ',irank,': invalid outtype nmat, nvec or n3d constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     types = MPI_BYTE

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, master_send_2dgw(irank), ierr)
     CALL MPI_Type_commit (master_send_2dgw(irank), ierr)

  end do  !loop over all worker ranks


  !to send these to the workers
  !if myrank == 0
  !do irank=1,nworkers
  !  mpi_send(mpi_bottom,1,master_send_2dgw(irank),irank,ktau,comm,send_recq(irank),ierr)
  !end do

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

END SUBROUTINE master_send_gw_types

SUBROUTINE master_recv_gw_types (comm,h,gwconv)
!called with if (myrank .eq. 0) master_recv_gw_types(...)

  USE mpi
  USE cable_mpicommon

  USE cable_def_types_mod

  IMPLICIT NONE
  INTEGER :: comm ! MPI communicator to talk to the workers
  !these are the global arrays
  REAL(r_2), DIMENSION(:,:), INTENT(IN) :: h,       &
                                           gwconv


  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: off, cnt
  INTEGER :: bidx, ierr
  INTEGER :: i

  ! MPI: allocate temp vectors used for marshalling
  ALLOCATE (blocks(ngw_master_recv_types))
  ALLOCATE (displs(ngw_master_recv_types))
  ALLOCATE (types(ngw_master_recv_types))

  do irank=1,nworkers   !need to define nworkers (nprocs -1 )

     if (irank .eq. 1) then  !offset in the global array
        off = 1
     else
        !need to correct for the halos
        off = sum(worker_nx_nohalo(1:irank-1))
     end if

     ! so master will make a type for each worker starting at
     ! worker_nx(worker_rank)
     cnt = worker_nx_nohalo(irank)  !doesn't include left and right halos no need for master to gather it
     r2len = cnt * extr2
 
     bidx = 0

  ! MPI: an hvector type for each vector, maddr contains displacements
     bidx = bidx + 1
     CALL MPI_Get_address (h(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     bidx = bidx + 1
     CALL MPI_Get_address (gwconv(off,1), displs(bidx), ierr)
     blocks(bidx) = r2len * mlat

     ! MPI: sanity check
     IF (bidx .ne. ngw_master_recv_types) THEN
        WRITE (*,*) 'worker ',irank,': invalid outtype nmat, nvec or n3d constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     types = MPI_BYTE

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, master_recv_2dgw(irank), ierr)
     CALL MPI_Type_commit (master_recv_2dgw(irank), ierr)

  end do  !loop over all worker ranks


  !to send these to the workers
  !if myrank == 0
  !do irank=1,nworkers
  !  mpi_send(mpi_bottom,1,master_send_2dgw(irank),irank,ktau,comm,send_recq(irank),ierr)
  !end do

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)


END SUBROUTINE master_recv_gw_types


SUBROUTINE master_send_data(comm,data_types,ktau) 
!called with if (myrank .eq. 0) master_send_data(...)
  implicit none
  integer, intent(in)               :: comm
  integer, dimension(:), intent(in) :: data_types  !master_send_2dgw
  integer, intent(in)               :: ktau  !time step labels message

  integer :: i,irank


  do irank=1,nworkers

     call mpi_isend(MPI_BOTTOM,1,data_types(irank),irank,ktau,&
                    comm,msend_req(irank),ierr)

  end do

END SUBROUTINE master_send_data(comm,data_types,ktau)


SUBROUTINE master_recv_data(comm,data_types,ktau)
!called with if (myrank .eq. 0) master_recv_data(...)

  implicit none
  integer, intent(in)               :: comm
  integer, dimension(:), intent(in) :: data_types  !master_recv_2dgw  
  integer, intent(in)               :: ktau  !time step labels message

  integer :: i,irank


  do irank=1,nworkers

     call mpi_irecv(MPI_BOTTOM,1,data_types(irank),irank,ktau,&
                    comm,mrecv_req(irank),ierr)

  end do

END SUBROUTINE master_recv_data

SUBROUTINE worker_send_data(comm,myrank,data_types,ktau)
  implicit none
  integer, intent(in)               :: comm
  integer, intent(in)               :: myrank
  integer, intent(in)               :: data_type  !worker_send_2dgw
  integer, intent(in)               :: ktau  !time step labels message

  integer :: i,irank


  call mpi_send(MPI_BOTTOM,1,data_type,0,ktau,comm,stat,ierr)  !stat and ierr not defined?

END SUBROUTINE worker_send_data


SUBROUTINE worker_recv_data(comm,myrank,data_types,ktau)
  implicit none
  integer, intent(in)               :: comm
  integer, intent(in)               :: myrank
  integer, intent(in)               :: data_type  !worker_recv_2dgw
  integer, intent(in)               :: ktau  !time step labels message

  integer :: i,irank


  call mpi_recv(MPI_BOTTOM,1,data_type,0,ktau,comm,stat,ierr)  !stat and ierr not defined?

END SUBROUTINE worker_recv_data

SUBROUTINE worker_to_worker_transfer(comm,myrank,lp_cnt,send_to_left,send_to_right,&
                                                 recv_from_left, recv_from_right)
  implicit none
  integer, intent(in) :: comm
  integer, intent(in) :: myrank
  integer, intent(in) :: lp_cnt
  integer, dimension(:), intent(in) :: send_to_left,send_to_right,&
                                       recv_from_left, recv_from_right

  if (myrank .lt. nworkers) then
    call mpi_send(MPI_BOTTOM,1,send_to_right(myrank),myrank+1,lp_cnt,comm,stat,ierr)
    call mpi_recv(MPI_BOTTOM,1,recv_from_right(myrank),myrank+1,lp_cnt,comm,stat,ierr)
  end if

  if (myrank .gt. 1) then
    call mpi_send(MPI_BOTTOM,1,send_to_left(myrank),myrank-1,lp_cnt,comm,stat,ierr) 
    call mpi_recv(MPI_BOTTOM,1,recv_from_left(myrank),myrank-1,lp_cnt,comm,stat,ierr)
  end if

END SUBROUTINE worker_to_worker_transfer



END MODULE mpi_2dgw_types
