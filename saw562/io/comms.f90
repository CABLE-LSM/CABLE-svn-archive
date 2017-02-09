!> 
!! Copyright 2017 ARC Centre of Excellence for Climate Systems Science
!!
!! \author  Scott Wales <scott.wales@unimelb.edu.au>
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.

module comms_mod
    use, intrinsic :: iso_c_binding
    use :: mpi_f08
    use :: cable_mpicommon, only: lpdecomp_t
    use log_mod

    type(err_code_t), parameter :: err_field_size = err_code_t(-100, "Field sizes don't match")

    type :: mode_t
        integer :: val
    end type

    type(mode_t), parameter :: mode_scatter = mode_t(1)
    type(mode_t), parameter :: mode_bcast   = mode_t(2)

    type :: field_t
        integer :: dim, shape(5)
        type(MPI_DATATYPE) :: type
        type(MPI_DATATYPE) :: basetype
        character(len=:), allocatable :: name
        type(mode_t) :: mode = mode_scatter

        ! Pointers to the field data (just the first value is enough for MPI)
        ! Only one will be valid, depending on field%type
        integer(kind=4), pointer :: i4_ptr => null()
        real(kind=4),    pointer :: r4_ptr => null()
        real(kind=8),    pointer :: r8_ptr => null()
        logical,         pointer :: l_ptr => null()
    contains
        ! Choose between scatter and bcast
        procedure :: distribute
        ! Scatter field from rank 0 to others
        procedure :: scatter => scatter_field
        ! Send identical data to all ranks
        procedure :: bcast => bcast_field
        procedure :: gather => gather_field
    end type

    type :: comms_t
        type(MPI_Comm) :: comm

        type(field_t), allocatable :: field(:)
        integer :: field_count

        type(lpdecomp_t), allocatable :: decomp(:)

        integer :: nland, npatch
    contains
        procedure :: init
        procedure :: free
        procedure :: resize

        ! MPI_Scatter all fields
        procedure :: scatter => scatter_comms
        procedure :: gather  => gather_comms

        ! Register a field_t
        !     comms%register_field(name, field)
        procedure :: register_field_t

        ! Register a field array
        !     comms%register_field(name, array)
        procedure :: register_field_i4_1
        procedure :: register_field_i4_2
        procedure :: register_field_i4_3
        procedure :: register_field_r4_1
        procedure :: register_field_r4_2
        procedure :: register_field_r4_3
        procedure :: register_field_r8_1
        procedure :: register_field_r8_2
        procedure :: register_field_r8_3
        procedure :: register_field_land
        procedure :: register_field_patch
        procedure :: register_field_l_1

        generic :: register_field &
            => register_field_t &
            , register_field_i4_1 &
            , register_field_i4_2 &
            , register_field_i4_3 &
            , register_field_r4_1 &
            , register_field_r4_2 &
            , register_field_r4_3 &
            , register_field_r8_1 &
            , register_field_r8_2 &
            , register_field_r8_3 &
            , register_field_land &
            , register_field_patch &
            , register_field_l_1
    end type

contains

    subroutine init(self, comm, decomp)
        class(comms_t), intent(inout) :: self
        integer :: comm
        type(lpdecomp_t), optional, intent(in) :: decomp(:)
        integer :: i

        self%comm = MPI_Comm(comm)

        ! Only needed on master
        if (present(decomp)) self%decomp = decomp

        self%nland = 0
        self%npatch = 0

        if (present(decomp)) then
            do i=1,size(decomp)
                self%nland  = self%nland + decomp(i)%nland
                self%npatch = self%npatch + decomp(i)%npatch
            end do
        end if
    end subroutine

    subroutine free(self)
        class(comms_t), intent(inout) :: self

    end subroutine

    subroutine register_field_t(self, field)
        ! Register a generic field object
        class(comms_t), intent(inout) :: self
        type(field_t), intent(in) :: field
        integer :: total, rank, comm_size
        integer(kind=8) :: hash, rhash

        call MPI_Comm_rank(self%comm, rank)
        call MPI_Comm_size(self%comm, comm_size)
        rhash = -1

        ! Check field names match across ranks
        hash = djb2(field%name)
        !write(*,*) rank, hash, rhash
        CALL MPI_Reduce(hash, rhash, 1, MPI_INTEGER8, MPI_BAND, 0, self%comm)
        if (rank == 0 .and. rhash /= hash) then
            write(*,*) "Error: Fields don't match, expected ", field%name
            call MPI_Abort(self%comm, -2)
        end if
       
        if (.not. allocated(self%field)) then
            allocate(self%field(20))
        end if

        if (size(self%field) == self%field_count) then
            call self%resize
        end if

        self%field_count = self%field_count + 1
        self%field(self%field_count) = field

        ! Do some sanity checks

        if (rank == 0) then
            total = 0
            call MPI_Reduce(MPI_IN_PLACE, total, 1, MPI_INTEGER, &
                MPI_SUM, 0, self%comm)

            if (field%mode%val == mode_scatter%val .and. total /= field%shape(1)) then
                write(*,*) rank, "Incompatible scatter field ", field%name, field%shape(1), total
                call log_error(err_field_size, field%name)

            else if (field%mode%val == mode_bcast%val .and. total /= (comm_size-1) * field%shape(1)) then
                write(*,*) rank, "Incompatible bcast field ", field%name, field%shape(1), total
                call log_error(err_field_size, field%name)

            end if
        else
            call MPI_Reduce(field%shape, total, 1, MPI_INTEGER, &
                MPI_SUM, 0, self%comm)
        end if
    end subroutine

    subroutine resize(self)
        ! Resize the field array
        class(comms_t), intent(inout) :: self
        type(field_t), allocatable :: temp(:)

        allocate(temp(size(self%field)*2))
        temp(1:size(self%field)) = self%field
        call move_alloc(temp, self%field)
    end subroutine

    subroutine scatter_comms(self)
        class(comms_t), intent(inout) :: self
        integer :: i

        DO i=1, self%field_count
            call self%field(i)%distribute(self, self%decomp)
            call MPI_Barrier(self%comm)
        END DO
    end subroutine

    subroutine gather_comms(self)
        class(comms_t), intent(inout) :: self
        integer :: i

        DO i=1, self%field_count
            call self%field(i)%gather(self, self%decomp)
        END DO
    end subroutine

    subroutine distribute(self, comms, decomp)
        class(field_t), intent(inout) :: self
        type(comms_t), intent(in) :: comms
        type(lpdecomp_t), allocatable, intent(in) :: decomp(:)

        if (self%mode%val == mode_bcast%val) then
            call self%bcast(comms%comm)
        else
            call self%scatter(comms, decomp)
        end if
    end subroutine

    subroutine gather_field(self, comms, decomp)
        ! Gather a single field
        ! Decomp is only required on comm_rank 0, to know what to send where
        class(field_t), intent(inout) :: self
        type(comms_t), intent(in) :: comms
        type(lpdecomp_t), allocatable, intent(in) :: decomp(:)
        integer, allocatable :: recvcounts(:), displs(:)
        integer :: comm_rank, comm_size, i
        integer :: sendcount

        sendcount  = self%shape(1) !product(self%shape(1:self%dim))

        call MPI_Comm_rank(comms%comm, comm_rank)
        call MPI_Comm_size(comms%comm, comm_size)
        allocate(recvcounts(comm_size), displs(comm_size)) 

        if (comm_rank == 0) then
            ! Send no data to rank 0
            recvcounts(1) = 0
            displs(1) = sendcount + 1

            ! Are we distributing land or patches?
            if (self%shape(1) == comms%npatch) then
                do i=2, comm_size
                    recvcounts(i) = decomp(i-1)%npatch
                    displs(i)     = (decomp(i-1)%patch0 -1)
                end do
            else if (self%shape(1) == comms%nland) then
                do i=2, comm_size
                    recvcounts(i) = decomp(i-1)%nland
                    displs(i)     = (decomp(i-1)%landp0 -1)
                end do
            else
                call MPI_Abort(comms%comm, -110)
            end if

            ! Rank 0 doesn't recieve any data
            sendcount = 0
        end if

        if (self%basetype == MPI_INTEGER4) then
            call MPI_Gatherv(self%i4_ptr, sendcount, self%type, &
                self%i4_ptr, recvcounts, displs, self%type, 0, comms%comm)
        else if (self%basetype == MPI_REAL4) then
            call MPI_Gatherv(self%r4_ptr, sendcount, self%type, &
                self%r4_ptr, recvcounts, displs, self%type, 0, comms%comm)
        else if (self%basetype == MPI_REAL8) then
            call MPI_Gatherv(self%r8_ptr, sendcount, self%type, &
                self%r8_ptr, recvcounts, displs, self%type, 0, comms%comm)
        else if (self%basetype == MPI_LOGICAL) then
            call MPI_Gatherv(self%l_ptr, sendcount, self%type, &
                self%l_ptr, recvcounts, displs, self%type, 0, comms%comm)
        else
            call MPI_Abort(comms%comm, -100)
        end if
    end subroutine

    subroutine scatter_field(self, comms, decomp)
        ! Scatter a single field
        ! Decomp is only required on comm_rank 0, to know what to send where
        class(field_t), intent(inout) :: self
        type(comms_t), intent(in) :: comms
        type(lpdecomp_t), allocatable, intent(in) :: decomp(:)
        integer, allocatable :: sendcounts(:), displs(:)
        integer :: comm_rank, comm_size, i
        integer :: patch_size, recvcount

        recvcount  = self%shape(1) ! product(self%shape(1:self%dim))

        call MPI_Comm_rank(comms%comm, comm_rank)
        call MPI_Comm_size(comms%comm, comm_size)
        allocate(sendcounts(comm_size), displs(comm_size)) 

        if (comm_rank == 0) then
            ! Send no data to rank 0
            sendcounts(1) = 0
            displs(1) = recvcount + 1

            ! Are we distributing land or patches?
            if (self%shape(1) == comms%npatch) then
                do i=2, comm_size
                    sendcounts(i) = decomp(i-1)%npatch
                    displs(i)     = (decomp(i-1)%patch0 -1)
                end do
            else if (self%shape(1) == comms%nland) then
                do i=2, comm_size
                    sendcounts(i) = decomp(i-1)%nland
                    displs(i)     = (decomp(i-1)%landp0 -1)
                end do
            else
                call MPI_Abort(comms%comm, -110)
            end if

            ! Rank 0 doesn't recieve any data
            recvcount = 0
        end if

        if (self%basetype == MPI_INTEGER4) then
            call MPI_Scatterv(self%i4_ptr, sendcounts, displs, self%type, &
                self%i4_ptr, recvcount, self%type, 0, comms%comm)
        else if (self%basetype == MPI_REAL4) then
            call MPI_Scatterv(self%r4_ptr, sendcounts, displs, self%type, &
                self%r4_ptr, recvcount, self%type, 0, comms%comm)
        else if (self%basetype == MPI_REAL8) then
            call MPI_Scatterv(self%r8_ptr, sendcounts, displs, self%type, &
                self%r8_ptr, recvcount, self%type, 0, comms%comm)
        else if (self%basetype == MPI_LOGICAL) then
            call MPI_Scatterv(self%l_ptr, sendcounts, displs, self%type, &
                self%l_ptr, recvcount, self%type, 0, comms%comm)
        else
            call MPI_Abort(comms%comm, -100)
        end if
    end subroutine

    subroutine bcast_field(self, comm)
        class(field_t), intent(inout) :: self
        type(MPI_Comm), intent(in) :: comm

        if (self%basetype == MPI_INTEGER4) then
            call MPI_Bcast(self%i4_ptr, self%shape(1), self%type, 0, comm)
        else if (self%basetype == MPI_REAL4) then
            call MPI_Bcast(self%r4_ptr, self%shape(1), self%type, 0, comm)
        else if (self%basetype == MPI_REAL8) then
            call MPI_Bcast(self%r8_ptr, self%shape(1), self%type, 0, comm)
        else if (self%basetype == MPI_LOGICAL) then
            call MPI_Bcast(self%l_ptr, self%shape(1), self%type, 0, comm)
        else
            call MPI_Abort(comm, -100)
        end if

    end subroutine

    subroutine test_scatter(send, sendcounts, displs, sendtype, &
            recv, recvcount, recvtype, root, comm)
        use mpi_f08
        real :: send, recv
        type(MPI_DATATYPE) :: sendtype, recvtype
        integer :: sendcounts(:), displs(:), recvcount
        integer :: root
        type(MPI_COMM) :: comm
            call MPI_Scatterv(send, sendcounts, displs, sendtype, &
                recv, recvcount, recvtype, root, comm)
    end subroutine

    function mpi_array_2d(source, shape) result(res)
        type(MPI_DATATYPE), intent(in) :: source
        integer :: shape(2)
        type(MPI_DATATYPE) :: res, tmp_type
        integer(kind=MPI_COUNT_KIND) :: lb, extent

        call MPI_Type_get_extent(source, lb, extent)
        call MPI_Type_vector(shape(2), 1, shape(1), source, tmp_type)
        call MPI_Type_create_resized(tmp_type, lb, extent, res)
        call MPI_Type_commit(res)
    end function

    function mpi_array_3d(source, shape) result(res)
        type(MPI_DATATYPE), intent(in) :: source
        integer :: shape(3)
        type(MPI_DATATYPE) :: res, tmp_type1, tmp_type2
        integer(kind=MPI_COUNT_KIND) :: lb, extent

        call MPI_Type_get_extent(source, lb, extent)
        call MPI_Type_vector(shape(2), 1, shape(1), source, tmp_type1)
        call MPI_Type_create_resized(tmp_type1, lb, extent, res)
        call MPI_Type_vector(shape(3), 1, shape(1)*shape(2), res, tmp_type2)
        call MPI_Type_create_resized(tmp_type2, lb, extent, res)
        call MPI_Type_commit(res)
    end function
    
    ! Register specific types of field by creating a generic object then registering that

    subroutine register_field_i4_1(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(kind=4), intent(in), pointer, contiguous :: ptr(:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_INTEGER4
        field%type = field%basetype
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%i4_ptr => ptr(1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_i4_2(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(kind=4), intent(in), pointer, contiguous  :: ptr(:,:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_INTEGER4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%i4_ptr => ptr(1,1)

        field%type = mpi_array_2d(field%basetype, field%shape)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_i4_3(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(kind=4), intent(in), pointer, contiguous  :: ptr(:,:,:)
        type(field_t) :: field
        type(MPI_DATATYPE) :: tmp_type1, tmp_type2
        integer(kind=MPI_COUNT_KIND) :: lb, extent

        field%name = name
        field%basetype = MPI_INTEGER4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%i4_ptr => ptr(1,1,1)

        field%type = mpi_array_3d(field%basetype, field%shape)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r4_1(self, name, ptr, mode)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), target, contiguous  :: ptr(:)
        type(field_t) :: field
        type(mode_t), intent(in), optional :: mode

        field%name = name
        field%basetype = MPI_REAL4
        field%type = field%basetype
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1)

        if (present(mode)) field%mode = mode

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r4_2(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), pointer, contiguous  :: ptr(:,:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_REAL4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1,1)

        field%type = mpi_array_2d(field%basetype, field%shape)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r4_3(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), pointer, contiguous  :: ptr(:,:,:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_REAL4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1,1,1)

        field%type = mpi_array_3d(field%basetype, field%shape)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r8_1(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=8), intent(in), pointer, contiguous  :: ptr(:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_REAL8
        field%type = field%basetype
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r8_ptr => ptr(1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r8_2(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=8), intent(in), pointer, contiguous  :: ptr(:,:)
        type(field_t) :: field
        type(MPI_DATATYPE) :: tmp_type
        integer(kind=MPI_COUNT_KIND) :: lb, extent

        field%name = name
        field%basetype = MPI_REAL8
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r8_ptr => ptr(1,1)

        field%type = mpi_array_2d(field%basetype, field%shape)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r8_3(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=8), intent(in), pointer, contiguous  :: ptr(:,:,:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_REAL8
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r8_ptr => ptr(1,1,1)

        field%type = mpi_array_3d(field%basetype, field%shape)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_l_1(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        logical, intent(in), pointer, contiguous  :: ptr(:)
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_LOGICAL
        field%type = field%basetype
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%l_ptr => ptr(1)

        call self%register_field_t(field)
    end subroutine

    subroutine register_field_r4t_1(self, name, ptr, mode)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), target, contiguous  :: ptr(:)
        type(field_t) :: field
        type(mode_t), intent(in), optional :: mode

        field%name = name
        field%basetype = MPI_REAL4
        field%type = field%basetype
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1)

        if (present(mode)) field%mode = mode

        call self%register_field_t(field)
    end subroutine

    subroutine register_field_patch(self, name, ptr)
        use cable_io_vars_module, only: patch_type
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(patch_type), intent(in), target, contiguous :: ptr(:)
        integer, parameter :: elements = 3
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_REAL4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1)%frac

        ! Sanity check
        if (SIZEOF(ptr(1)) /= SIZEOF(ptr(1)%frac) * elements) call MPI_Abort(MPI_COMM_WORLD, -200)
        call MPI_Type_vector(1, elements, 0, field%basetype, field%type)
        call MPI_Type_commit(field%type)

        call self%register_field_t(field)
    end subroutine

    subroutine register_field_land(self, name, ptr)
        use cable_io_vars_module, only: land_type
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        type(land_type), intent(in), target, contiguous :: ptr(:)
        integer, parameter :: elements = 5
        type(field_t) :: field

        field%name = name
        field%basetype = MPI_INTEGER4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%i4_ptr => ptr(1)%nap

        ! Sanity check
        if (SIZEOF(ptr(1)) /= SIZEOF(ptr(1)%nap) * elements) call MPI_Abort(MPI_COMM_WORLD, -200)
        call MPI_Type_vector(1, elements, 0, field%basetype, field%type)
        call MPI_Type_commit(field%type)

        call self%register_field_t(field)
    end subroutine

    ! djb2 hash function http://www.cse.yorku.ca/~oz/hash.html
    pure function djb2(string)
        character(len=*), intent(in) :: string
        integer(kind=8) :: djb2
        integer :: i, c

        djb2 = 5381

        do i=1,len(string)
            c = ichar(string(i:i))
            djb2 = 33 * djb2 + c
        end do
    end function

end module
