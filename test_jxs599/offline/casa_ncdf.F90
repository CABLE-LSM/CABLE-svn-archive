!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: handles additional, dynamically decided diagnostic output from model.
!          permanently used for bitwise identical testing. more applications
!          will follow.
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Currently stripped down version of cable_diag here. will be
!          re-implemented in time.
!
! ==============================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ USE this module in any subr. you wish to write vars from.             +++!
!+++ x is typically the number of landpoints(tiles). binary file is        +++!
!+++ then appended every timestep with the new foo(x_i)                    +++!
!+++                                                                       +++!
!+++ CALL syntax:                                                          +++!
!+++                                                                       +++!
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, var1 )     +++!
!+++                                                                       +++!
!+++ output binaries can be interpreted from the command line              +++!
!+++ using a suite of tools. Currently, only zero_diff.ksh is supported.   +++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!CABLE_LSM:This has to be commented for offline
!#define UM_BUILD YES
MODULE casa_ncdf_module
   
   IMPLICIT NONE

#ifndef UM_BUILD
  interface put_var_nc
     module procedure put_var_ncr1, put_var_ncr2, put_var_ncr3
  end interface put_var_nc

  interface get_var_nc
     module procedure get_var_ncr2, get_var_ncr3
  end interface get_var_nc

#endif

CONTAINS

#ifndef UM_BUILD
  subroutine def_dims(nd, ncid, dimID, dim_len, dim_name )
    use netcdf
    implicit none
    integer, intent(in) :: nd, ncid
    character(len=*), dimension(:), intent(in) :: dim_name
    integer, dimension(:), intent(out) :: dimID
    integer, dimension(:), intent(in) :: dim_len
    integer :: j, ncok

    do j=1, nd
       ncok = NF90_DEF_DIM(ncid, trim(dim_name(j)), dim_len(j), dimID(j) )
       if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def dim ', dim_name(j))
    enddo

    return
  end subroutine def_dims




  subroutine def_vars(nv, ncid,  xtype, dimID, var_name,varID )
    use netcdf
    implicit none
    integer, intent(in) :: nv, ncid, xtype
    integer, dimension(:), intent(in) :: dimID
    integer, dimension(:), intent(inout) :: varID
    character(len=*), dimension(:), intent(in) :: var_name
    integer :: j, ncok

    ! lat
    ncok = NF90_DEF_VAR( ncid, trim(var_name(1)), xtype, &
         (/ dimID(1) /), varID(1))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(1))

    ! lon
    ncok = NF90_DEF_VAR(ncid, trim(var_name(2)), xtype, &
         (/ dimID(1) /), varID(2))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(2))

    ! tairk
    ncok = NF90_DEF_VAR(ncid, trim(var_name(3)), xtype, &
         (/ dimID(1), dimID(3) /), varID(3))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(3))

    !tsoil
    ncok = NF90_DEF_VAR(ncid, trim(var_name(4)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(4))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(4))

    ! moist
    ncok = NF90_DEF_VAR(ncid, trim(var_name(5)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(5))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(5))

    !cgpp
    ncok = NF90_DEF_VAR(ncid, trim(var_name(6)), xtype, &
         (/ dimID(1), dimID(3)/), varID(6))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(6))

    !crmplant
    ncok = NF90_DEF_VAR(ncid, trim(var_name(7)), xtype, &
         (/ dimID(1), dimID(2),dimID(3)/), varID(7))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(7))

    !phenphase
    ncok = NF90_DEF_VAR(ncid, trim(var_name(8)), xtype, &
         (/ dimID(1), dimID(3)/), varID(8))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(8))

    !doyphase1
    ncok = NF90_DEF_VAR(ncid, trim(var_name(9)), xtype, &
         (/ dimID(1), dimID(3)/), varID(9))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(9))

    !doyphase2
    ncok = NF90_DEF_VAR(ncid, trim(var_name(10)), xtype, &
         (/ dimID(1), dimID(3)/), varID(10))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(10))

    !doyphase3
    ncok = NF90_DEF_VAR(ncid, trim(var_name(11)), xtype, &
         (/ dimID(1), dimID(3)/), varID(11))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(11))

    !doyphase4
    ncok = NF90_DEF_VAR(ncid, trim(var_name(12)), xtype, &
         (/ dimID(1), dimID(3)/), varID(12))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(12))


    !mtemp
    ncok = NF90_DEF_VAR(ncid, trim(var_name(13)), xtype, &
         (/ dimID(1),dimID(3)/), varID(13))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(13))

    !Ndep
    ncok = NF90_DEF_VAR(ncid, trim(var_name(14)), xtype, &
         (/ dimID(1),dimID(3)/), varID(14))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def var ', var_name(14))

    return
  end subroutine def_vars

  subroutine def_var_atts( ncfile_in, ncid, varID )
    use netcdf
    implicit none
    character(len=*), intent(in) :: ncfile_in
    integer, intent(in):: ncid       ! netcdf file ID
    integer, dimension(:), intent(in) :: varID ! (1) ~ tvair, (2) ~ pmb
    integer :: j, ncok
    character(len=10) dummy

    write(dummy,11) varID(1)
11  format(i2)
    ncok = NF90_PUT_ATT(ncid, nf90_global, "Title", "Forcing for define_air subroutine")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', ncfile_in)
    ncok = NF90_PUT_ATT(ncid, varID(3), "longname", "air temperature within canopy")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', dummy)
    ncok = NF90_PUT_ATT(ncid, varID(3), "units", "K")
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'def att ', dummy)

    write(dummy,11) varID(2)


    return
  end subroutine def_var_atts


  subroutine put_var_ncr1(ncid, var_name, var )
    use netcdf
    use cable_def_types_mod, only : mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real, dimension(:),intent(in) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID,j

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1/), &
         count=(/mp/) )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

  end subroutine put_var_ncr1


  subroutine put_var_ncr2(ncid, var_name, var, n_call )
    use netcdf
    use cable_def_types_mod, only : r_2, mp
    implicit none
    character(len=*), intent(in) ::  var_name
    real(r_2), dimension(:),intent(in) :: var
    integer, intent(in) :: ncid, n_call
    integer :: ncok, varID

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,n_call /), &
         count=(/mp,1/) )

    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

  end subroutine put_var_ncr2

  !soil vars
  subroutine put_var_ncr3(ncid, var_name, var, n_call, nl)
    use netcdf
    use cable_def_types_mod, only : r_2, mp, ms
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:,:),intent(in) :: var
    integer, intent(in) :: ncid, n_call, nl
    integer :: ncok, varID

    ncok = NF90_INQ_VARID( ncid, var_name, varId )
    IF( ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name )

    ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,1,n_call /), &
         count=(/mp,nl,1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)

    return
  end subroutine put_var_ncr3



  subroutine get_var_ncr2(ncid, var_name, var, n_call )
    use netcdf
    use cable_def_types_mod, only : r_2,mp
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:),intent(out) :: var
    integer, intent(in) :: ncid
    integer :: ncok, varID, n_call
    real, dimension(mp) :: temp

    temp = 0.

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)
    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,n_call/), &
         count=(/mp,1/) )

    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'getting var ', var_name)

    var = real( temp, r_2 )
  end subroutine get_var_ncr2

  subroutine get_var_ncr3(ncid, var_name, var, n_call, nl )
    use netcdf
    use cable_def_types_mod, only : r_2, mp, ms
    implicit none
    character(len=*), intent(in) :: var_name
    real(r_2), dimension(:,:),intent(out) :: var
    integer, intent(in) :: ncid, n_call, nl
    integer :: ncok, varID
    real, dimension(mp,1:nl) :: temp

    ncok = NF90_INQ_VARID(ncid, var_name, varId )
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'inquire var ', var_name)

    ncok = NF90_GET_VAR(ncid, varId, temp, start=(/1,1,n_call /), &
         count=(/mp, nl, 1/))
    if (ncok /= nf90_noerr ) call stderr_nc(ncok,'putting var ', var_name)
    var = real( temp, r_2 )
  end subroutine get_var_ncr3



  subroutine stderr_nc(status,message, var)
    use netcdf
    character(len=*), intent(in) :: message, var
    INTEGER, INTENT(IN) :: status
    character(len=7) :: err_mess
    err_mess = 'ERROR:'
    print *, (err_mess//message), var
    PRINT*,NF90_STRERROR(status)
    stop
  end subroutine stderr_nc
#endif

END MODULE casa_ncdf_module



