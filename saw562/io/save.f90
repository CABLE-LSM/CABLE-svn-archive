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
module save_mod
        integer, parameter :: land_points = 15238

    contains

    subroutine save_fields(met, veg)
        use mpi_f08
        use cable_def_types_mod, only: met_type, veg_parameter_type
        implicit none
        type(met_type), intent(in) :: met
        type(veg_parameter_type), intent(in) :: veg
        integer :: unit

        open(newunit=unit, file='tmp', form='UNFORMATTED', status='REPLACE')
        write(unit) met%fsd
        write(unit) met%tk
        write(unit) met%pmb
        write(unit) met%qv
        write(unit) met%ua
        write(unit) met%precip
        write(unit) met%precip_sn
        write(unit) met%fld
        write(unit) met%ca
        write(unit) met%coszen
        write(unit) met%Ndep
        write(unit) veg%vlai
        write(unit) met%year
        write(unit) met%moy
        write(unit) met%doy
        write(unit) met%hod

        call MPI_Abort(MPI_COMM_WORLD, -1)
    end subroutine

    subroutine load_fields(met, veg, filename)
        use cable_def_types_mod, only: met_type, veg_parameter_type, alloc_cbm_var
        implicit none
        type(met_type), intent(out) :: met
        type(veg_parameter_type), intent(out) :: veg
        character(len=*), intent(in) :: filename
        integer :: unit

        call alloc_cbm_var(met, 15238)
        call alloc_cbm_var(veg, 15238)

        open(newunit=unit, file=filename, form='UNFORMATTED', status='OLD')
        read(unit) met%fsd
        read(unit) met%tk
        read(unit) met%pmb
        read(unit) met%qv
        read(unit) met%ua
        read(unit) met%precip
        read(unit) met%precip_sn
        read(unit) met%fld
        read(unit) met%ca
        read(unit) met%coszen
        read(unit) met%Ndep
        read(unit) veg%vlai
        read(unit) met%year
        read(unit) met%moy
        read(unit) met%doy
        read(unit) met%hod

    end subroutine

    subroutine diff_fields(file_a, file_b)
        use cable_def_types_mod, only: met_type, veg_parameter_type

        character(len=*), intent(in) :: file_a, file_b
        type(met_type):: met_a, met_b
        type(veg_parameter_type):: veg_a,veg_b
        integer :: i

        call load_fields(met_a, veg_a, file_a)
        call load_fields(met_b, veg_b, file_b)

        ! call diff_r(met_a%fsd,met_b%fsd,'met%fsd')
        call diff_r(met_a%tk,met_b%tk,'met%tk')
        call diff_r(met_a%pmb,met_b%pmb,'met%pmb')
        call diff_r(met_a%qv,met_b%qv,'met%qv')
        call diff_r(met_a%ua,met_b%ua,'met%ua')
        call diff_r(met_a%precip,met_b%precip,'met%precip')
        call diff_r(met_a%precip_sn,met_b%precip_sn,'met%precip_sn')
        call diff_r(met_a%fld,met_b%fld,'met%fld')
        call diff_r(met_a%ca,met_b%ca,'met%ca')
        call diff_r(met_a%coszen,met_b%coszen,'met%coszen')
        call diff_r(met_a%Ndep,met_b%Ndep,'met%Ndep')
        call diff_r(veg_a%vlai,veg_b%vlai,'veg%vlai')
        ! call diff_r(met_a%year,met_b%year,'met%year')
        ! call diff_r(met_a%moy,met_b%moy,'met%moy')
        call diff_r(met_a%doy,met_b%doy,'met%doy')
        call diff_r(met_a%hod,met_b%hod,'met%hod')
    end subroutine

    subroutine diff_r(a, b, name)
        real, intent(in) :: a(:), b(:)
        character(len=*), intent(in) :: name

        integer c
        integer i(1)

        c = count(a /= b)
        if (c /= 0) then
            i = maxloc(abs(a - b), a /= b)
            write(*,*) name, c, i, a(i(1)), b(i(1))
        end if
    end subroutine
end module
