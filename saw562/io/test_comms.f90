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

module test_comms_mod

contains

    subroutine setup_decomp(decomp, land_points, patches_per_land, comm)
        use MPI_f08
        use cable_mpicommon, only: lpdecomp_t

        type(lpdecomp_t), allocatable, intent(out) :: decomp(:)
        integer, intent(in) :: land_points, patches_per_land
        type(MPI_Comm), intent(in) :: comm

        integer :: rank, comm_size

        call MPI_Comm_rank(comm, rank)
        call MPI_Comm_size(comm, comm_size)
        allocate(decomp(comm_size-1))

        decomp(:)%nland = land_points / size(decomp)
        decomp(1)%nland = decomp(1)%nland + mod(land_points, size(decomp))

        decomp(:)%npatch = patches_per_land * decomp(:)%nland
    end subroutine

    subroutine test_return()
        use MPI_f08
        use cable_mpicommon, only: lpdecomp_t
        use comms_mod

        type(lpdecomp_t), allocatable :: decomp(:)
        integer, target, allocatable :: field(:), orig(:)
        integer, pointer :: ptr(:)
        integer :: land_points
        integer :: rank, i
        type(comms_t) :: comms

        call MPI_Comm_rank(MPI_COMM_WORLD, rank)

        land_points = 11
        call setup_decomp(decomp, land_points, 1, MPI_COMM_WORLD)
        if (rank/=0) write(*,*) rank, decomp(rank)%nland

        if (rank/=0) allocate(field(decomp(rank)%nland))
        if (rank==0) allocate(field(land_points))

        do i=1,size(field)
            field(i) = i
        end do
        orig = field
        ptr => field

        call comms%init(MPI_COMM_WORLD%mpi_val, decomp)
        call comms%register_field('field', ptr)
        call comms%scatter()
        if (rank == 0) field = -999
        call comms%gather()

        if (rank == 0) then
            do i=1,land_points
                if (field(i) /= orig(i)) then
                    write(*,*) 'Return error at',i, 'expect', orig(i), 'got',field(i)
                end if
            end do
        end if
    end subroutine

end module

program test_comms
    use test_comms_mod
    use comms_mod
    use mpi_f08
    use :: cable_mpicommon, only: lpdecomp_t
    use cable_def_types_mod, only: met_type, alloc_cbm_var

    type(comms_t) :: comms
    type(lpdecomp_t), allocatable :: decomp(:)
    integer :: rank, comm_size
    integer :: npatch, i, j, k
    type(met_type) :: met
    !integer , parameter :: npatch_total = 15238
    integer, parameter :: land_points = 9
    integer, parameter :: patches_per_land = 1
    integer, parameter :: npatch_total = land_points * patches_per_land
    integer, pointer :: field(:,:,:)

    call MPI_Init()
    call MPI_Comm_rank(MPI_COMM_WORLD, rank)
    call MPI_Comm_size(MPI_COMM_WORLD, comm_size)

    call setup_decomp(decomp, land_points, patches_per_land, MPI_COMM_WORLD)

    if (rank == 0) then
        write(*,*) rank, "patches", npatch_total
        call alloc_cbm_var(met, npatch_total)
        allocate(field(npatch_total, 2, 2))

        do k=1,size(field,3)
            do j=1,size(field,2)
                field(:,j,k) = j * size(field,3) + i
            end do
        end do

        met%fsd(:,1)  = 1
        met%fsd(:,2)  = 2
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
        write(*,*) rank, "patches", decomp(rank)%npatch
        call alloc_cbm_var(met, decomp(rank)%npatch)
        allocate(field(decomp(rank)%npatch, 2, 2))

        field = -999
        met%fsd       = -999
        met%tk        = -999
        met%pmb       = -999
        met%qv        = -999
        met%ua        = -999
        met%precip    = -999
        met%precip_sn = -999
        met%fld       = -999
        met%ca        = -999
        met%coszen    = -999
        met%Ndep      = -999
        met%year      = -999
        met%moy       = -999
        met%doy       = -999
        met%hod       = -999
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
    call comms%register_field('field', field)

    call comms%scatter()

    if (any(met%fsd(:,1) /= 1)) write(*,*) rank, "Error met%fsd ", count(met%fsd(:,1) /= 1), met%fsd(1:5,1)
    if (any(met%fsd(:,2) /= 2)) write(*,*) rank, "Error met%fsd ", count(met%fsd(:,2) /= 2), met%fsd(1:5,2)
    if (any(met%tk /= 2)) write(*,*) rank, "Error met%tk ", count(met%tk /= 2), met%tk(1:5)
    if (any(met%pmb /= 3)) write(*,*) rank, "Error met%pmb ", count(met%pmb /= 3), met%pmb(1:5)
    if (any(met%qv /= 4)) write(*,*) rank, "Error met%qv ", count(met%qv /= 4), met%qv(1:5)
    if (any(met%ua /= 5)) write(*,*) rank, "Error met%ua ", count(met%ua /= 5), met%ua(1:5)
    if (any(met%precip /= 6)) write(*,*) rank, "Error met%precip ", count(met%precip /= 6), met%precip(1:5)
    if (any(met%precip_sn /= 7)) write(*,*) rank, "Error met%precip_sn ", count(met%precip_sn /= 7), met%precip_sn(1:5)
    if (any(met%fld /= 8)) write(*,*) rank, "Error met%fld ", count(met%fld /= 8), met%fld(1:5)
    if (any(met%ca /= 9)) write(*,*) rank, "Error met%ca ", count(met%ca /= 9), met%ca(1:5)
    if (any(met%coszen /= 10)) write(*,*) rank, "Error met%coszen ", count(met%coszen /= 10), met%coszen(1:5)
    if (any(met%Ndep /= 11)) write(*,*) rank, "Error met%Ndep ", count(met%Ndep /= 11), met%Ndep(1:5)
    if (any(met%year /= 12)) write(*,*) rank, "Error met%year ", count(met%year /= 12), met%year(1:5)
    if (any(met%moy /= 13)) write(*,*) rank, "Error met%moy ", count(met%moy /= 13), met%moy(1:5)
    if (any(met%doy /= 14)) write(*,*) rank, "Error met%doy ", count(met%doy /= 14), met%doy(1:5)
    if (any(met%hod /= 15)) write(*,*) rank, "Error met%hod ", count(met%hod /= 15), met%hod(1:5)

        do k=1,size(field,3)
            do j=1,size(field,2)
                if (any(field(:,j,k) /= j * size(field,3) + i)) write(*,*) rank, "Error field ", field(:,j,k)
            end do
        end do

        call test_return()

    call MPI_Finalize()
end program


