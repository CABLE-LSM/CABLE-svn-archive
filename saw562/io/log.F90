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

module log_mod

    type :: err_code_t
        integer :: code
        character(len=1024) :: message
    end type

    type(err_code_t), parameter :: err_success = err_code_t(0, "Success")

    type(err_code_t), parameter :: err_missing_optional = err_code_t(500, "OPTIONAL not present")

contains

    subroutine traceback()
#ifdef __INTEL_COMPILER
        use ifcore
        call tracebackqq(user_exit_code=-1)
#else
#message "Error log traceback not supported"
#endif
    end subroutine

    subroutine log_warning(code, message)
        use MPI_f08
        type(err_code_t), intent(in) :: code
        character(len=*), intent(in), optional :: message
        integer :: rank

        call MPI_Comm_rank(MPI_COMM_WORLD, rank)
        if (present(message)) then
            write(*,*) rank, 'WARNING', code%code, '(', trim(code%message), ') ', message
        else
            write(*,*) rank, 'WARNING', code%code, '(', trim(code%message), ') '
        end if

        call traceback()
    end subroutine

    subroutine log_error(code, message)
        use MPI_f08
        type(err_code_t), intent(in) :: code
        character(len=*), intent(in), optional :: message
        integer :: rank

        call MPI_Comm_rank(MPI_COMM_WORLD, rank)
        if (present(message)) then
            write(*,*) rank, 'ERROR', code%code, '(', trim(code%message), ') ', message
        else
            write(*,*) rank, 'ERROR', code%code, '(', trim(code%message), ') '
        end if

        call traceback()
        call MPI_Abort(MPI_COMM_WORLD, code%code)
    end subroutine
end module
