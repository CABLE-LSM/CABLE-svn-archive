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

    type :: field_t
        integer :: dim, shape(5)
        type(MPI_DATATYPE) :: type
        character(len=:), allocatable :: name

        ! Pointers to the field data (just the first value is enough for MPI)
        ! Only one will be valid, depending on field%type
        integer(kind=4), pointer :: i4_ptr => null()
        real(kind=4),    pointer :: r4_ptr => null()
        real(kind=8),    pointer :: r8_ptr => null()
    contains
        ! Scatter field from rank 0 to others
        procedure :: scatter => scatter_field
    end type

    type :: comms_t
        type(MPI_Comm) :: comm

        type(field_t), allocatable :: field(:)
        integer :: field_count

        type(lpdecomp_t), allocatable :: decomp(:)
    contains
        procedure :: init
        procedure :: free
        procedure :: resize

        ! MPI_Scatter all fields
        procedure :: scatter

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
            , register_field_r8_3
    end type

contains

    subroutine init(self, comm, decomp)
        class(comms_t), intent(inout) :: self
        integer :: comm
        type(lpdecomp_t), optional, intent(in) :: decomp(:)

        self%comm = MPI_Comm(comm)

        ! Only needed on master
        if (present(decomp)) self%decomp = decomp
    end subroutine

    subroutine free(self)
        class(comms_t), intent(inout) :: self

    end subroutine

    subroutine register_field_t(self, field)
        ! Register a generic field object
        class(comms_t), intent(inout) :: self
        type(field_t), intent(in) :: field
        integer :: total, rank
       
        if (.not. allocated(self%field)) then
            allocate(self%field(20))
        end if

        if (size(self%field) == self%field_count) then
            call self%resize
        end if

        self%field_count = self%field_count + 1
        self%field(self%field_count) = field

        ! Do some sanity checks
        call MPI_Comm_rank(self%comm, rank)

        if (rank == 0) then
            total = 0
            call MPI_Reduce(MPI_IN_PLACE, total, 1, MPI_INTEGER, &
                MPI_SUM, 0, self%comm)
            if (total /= field%shape(1)) then
                write(*,*) rank, "Incompatible field ", field%name, field%shape(1), total
                call MPI_Abort(self%comm, -100)
            end if
        else
            call MPI_Reduce(field%shape, total, 1, MPI_INTEGER, &
                MPI_SUM, 0, self%comm)
            write(*,*) rank, "field", field%name, field%shape
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

    subroutine scatter(self)
        class(comms_t), intent(inout) :: self
        integer :: i

        DO i=1, self%field_count
            call self%field(i)%scatter(self%comm, self%decomp)
        END DO
    end subroutine

    subroutine scatter_field(self, comm, decomp)
        ! Scatter a single field
        ! Decomp is only required on rank 0, to know what to send where
        class(field_t), intent(inout) :: self
        type(MPI_Comm), intent(in) :: comm
        type(lpdecomp_t), allocatable, intent(in) :: decomp(:)
        integer, allocatable :: sendcounts(:), displs(:)
        integer :: rank, size, i
        integer :: patch_size, recvcount

        patch_size = product(self%shape(2:self%dim))
        recvcount  = product(self%shape(1:self%dim))

        call MPI_Comm_rank(comm, rank)
        call MPI_Comm_size(comm, size)
        allocate(sendcounts(size), displs(size)) 

        if (rank == 0) then
            recvcount = 0
            sendcounts(1) = 0
            displs(1) = 0

            do i=2, size
                sendcounts(i) = patch_size * decomp(i-1)%npatch
                displs(i)     = patch_size * decomp(i-1)%patch0
            end do
        end if
        write(*,*) 'scatter ', self%name, rank, recvcount

        if (self%type == MPI_INTEGER4) then
            call MPI_Scatterv(self%i4_ptr, sendcounts, displs, self%type, &
                self%i4_ptr, recvcount, self%type, 0, comm)
        else if (self%type == MPI_REAL4) then
            call MPI_Scatterv(self%r4_ptr, sendcounts, displs, self%type, &
                self%r4_ptr, recvcount, self%type, 0, comm)
        else if (self%type == MPI_REAL8) then
            call MPI_Scatterv(self%r8_ptr, sendcounts, displs, self%type, &
                self%r8_ptr, recvcount, self%type, 0, comm)
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

        
    
    ! Register specific types of field by creating a generic object then registering that

    subroutine register_field_i4_1(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(kind=4), intent(in), pointer, contiguous :: ptr(:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_INTEGER4
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
        field%type = MPI_INTEGER4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%i4_ptr => ptr(1,1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_i4_3(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        integer(kind=4), intent(in), pointer, contiguous  :: ptr(:,:,:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_INTEGER4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%i4_ptr => ptr(1,1,1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r4_1(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), pointer, contiguous  :: ptr(:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_REAL4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r4_2(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), pointer, contiguous  :: ptr(:,:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_REAL4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1,1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r4_3(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=4), intent(in), pointer, contiguous  :: ptr(:,:,:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_REAL4
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r4_ptr => ptr(1,1,1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r8_1(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=8), intent(in), pointer, contiguous  :: ptr(:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_REAL8
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

        field%name = name
        field%type = MPI_REAL8
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r8_ptr => ptr(1,1)

        call self%register_field_t(field)
    end subroutine
    subroutine register_field_r8_3(self, name, ptr)
        class(comms_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=8), intent(in), pointer, contiguous  :: ptr(:,:,:)
        type(field_t) :: field

        field%name = name
        field%type = MPI_REAL8
        field%dim = size(shape(ptr))
        field%shape(1:field%dim) = shape(ptr) 
        field%r8_ptr => ptr(1,1,1)

        call self%register_field_t(field)
    end subroutine

end module
