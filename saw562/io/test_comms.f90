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
    npatch = 15238
    call alloc_cbm_var(met, npatch)

    if (rank == 0) then
        allocate(decomp(size-1))

        do i=1, size-1
            decomp(i)%npatch = npatch
            decomp(i)%patch0 = npatch*(i-1) + 1
        end do
        allocate(field(npatch*(size-1)))
        met%fsd       = 1
        met%tk        = 2
        met%pmb       = 3
        met%qv        = 4
        met%ua        = 5
        met%precip    = 6
        met%precip_sn = 7
        met%fld       = 8
        met%ca        = 9
        met%coszen    = 10
        met%Ndep      = 11
        met%year      = 12
        met%moy       = 13
        met%doy       = 14
        met%hod       = 15
        call comms%init(MPI_COMM_WORLD%mpi_val, decomp)
    else
        allocate(field(npatch))
        field = 0
        call comms%init(MPI_COMM_WORLD%mpi_val)
    end if


        call comms%register_field('met%fsd', met%fsd)
        call comms%register_field('met%tk', met%tk)
        call comms%register_field('met%pmb', met%pmb)
        call comms%register_field('met%qv', met%qv)
        call comms%register_field('met%ua', met%ua)
        call comms%register_field('met%precip', met%precip)
        call comms%register_field('met%precip_sn', met%precip_sn)
        call comms%register_field('met%fld', met%fld)
        call comms%register_field('met%ca', met%ca)
        call comms%register_field('met%coszen', met%coszen)
        call comms%register_field('met%Ndep', met%Ndep)
        call comms%register_field('met%year', met%year)
        call comms%register_field('met%moy', met%moy)
        call comms%register_field('met%doy', met%doy)
        call comms%register_field('met%hod', met%hod)
    
    !call comms%register_field('foo', field)
    call comms%scatter()


        if (any(met%fsd /= 1)) write(*,*) "Error met%fsd ", count(met%fsd /= 1), met%fsd(1:5,1)
        if (any(met%tk /= 2)) write(*,*) "Error met%tk ", count(met%tk /= 2), met%tk(1:5)
        if (any(met%pmb /= 3)) write(*,*) "Error met%pmb ", count(met%pmb /= 3), met%pmb(1:5)
        if (any(met%qv /= 4)) write(*,*) "Error met%qv ", count(met%qv /= 4), met%qv(1:5)
        if (any(met%ua /= 5)) write(*,*) "Error met%ua ", count(met%ua /= 5), met%ua(1:5)
        if (any(met%precip /= 6)) write(*,*) "Error met%precip ", count(met%precip /= 6), met%precip(1:5)
        if (any(met%precip_sn /= 7)) write(*,*) "Error met%precip_sn ", count(met%precip_sn /= 7), met%precip_sn(1:5)
        if (any(met%fld /= 8)) write(*,*) "Error met%fld ", count(met%fld /= 8), met%fld(1:5)
        if (any(met%ca /= 9)) write(*,*) "Error met%ca ", count(met%ca /= 9), met%ca(1:5)
        if (any(met%coszen /= 10)) write(*,*) "Error met%coszen ", count(met%coszen /= 10), met%coszen(1:5)
        if (any(met%Ndep /= 11)) write(*,*) "Error met%Ndep ", count(met%Ndep /= 11), met%Ndep(1:5)
        if (any(met%year /= 12)) write(*,*) "Error met%year ", count(met%year /= 12), met%year(1:5)
        if (any(met%moy /= 13)) write(*,*) "Error met%moy ", count(met%moy /= 13), met%moy(1:5)
        if (any(met%doy /= 14)) write(*,*) "Error met%doy ", count(met%doy /= 14), met%doy(1:5)
        if (any(met%hod /= 15)) write(*,*) "Error met%hod ", count(met%hod /= 15), met%hod(1:5)

    call MPI_Finalize()
end program


