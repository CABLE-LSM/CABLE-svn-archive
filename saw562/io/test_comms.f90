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

program test_comms
    use comms_mod
    use mpi_f08
    use :: cable_mpicommon, only: lpdecomp_t
    use cable_def_types_mod, only: met_type, alloc_cbm_var
    
    type(comms_t) :: comms
    type(lpdecomp_t), allocatable :: decomp(:)
    integer :: rank, size
    integer :: npatch, i
    integer, allocatable :: field(:)
    type(met_type) :: met

    call MPI_Init()
    call MPI_Comm_rank(MPI_COMM_WORLD, rank)
    call MPI_Comm_size(MPI_COMM_WORLD, size)
    npatch = 100
        call alloc_cbm_var(met, npatch)

    if (rank == 0) then
        allocate(decomp(size-1))

        do i=1, size-1
            decomp(i)%npatch = npatch
            decomp(i)%patch0 = npatch*(i-1)
        end do
        allocate(field(npatch*(size-1)))
        field = 999
        met%tk = 123
        call comms%init(MPI_COMM_WORLD%mpi_val, decomp)
    else
        allocate(field(npatch))
        field = 0
        call comms%init(MPI_COMM_WORLD%mpi_val)
    end if


    call comms%register_field('met%tk', met%tk)
    if (.not. associated(comms%field(1)%r4_ptr, met%tk(1))) then
        write(*,*) 'not assoc'
    end if
    
    !call comms%register_field('foo', field)
    call comms%scatter()

    !if (.not. all(field == 999)) then
    !    write(*,*) rank, "Error"
    !    write(*,*) field
    !end if

    if (.not. all(met%tk == 123)) then
        write(*,*) rank, "Error met%tk", count(met%tk /= 123)
    end if

    call MPI_Finalize()
end program


