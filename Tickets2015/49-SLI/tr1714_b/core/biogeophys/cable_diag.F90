!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
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


MODULE cable_diag_module
  use cable_def_types_mod, only : r_2
   IMPLICIT NONe
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1
  
   !--- subrs overloaded to respond to call cable_diag 
   INTERFACE cable_diag
      MODULE PROCEDURE cable_diag1
   END INTERFACE cable_diag
   
  interface put_var_nc
     module procedure put_var_ncr1, put_var_ncr2, put_var_ncr3
  end interface put_var_nc

  interface get_var_nc
     module procedure get_var_ncr2, get_var_ncr3
  end interface get_var_nc

CONTAINS

!==========================================================================!
! cable_diag1/2/3 call subrs to write filename.dat which contains description
! of data and format etc., and filename.bin containing the data   
!==========================================================================!

SUBROUTINE cable_diag1( Nvars, basename, dimx, dimy, timestep, node, &
                        vname1, var1 )
   integer, intent(in) :: Nvars,dimx, dimy, timestep,node
   real, intent(in), dimension(:) :: var1
   integer :: i=0
   character(len=*), intent(in) :: basename, vname1
   character(len=30) :: filename, chnode
  
      write(chnode,10) node
   10 format(i2.2)   
      filename=trim(trim(basename)//trim(chnode))
      
      if (timestep == 1) & 
         call cable_diag_desc1( Nvars, trim(filename), dimx, dimy, vname1 )
      
      call cable_diag_data1( Nvars, trim(filename), dimx, timestep, dimy, &
                             var1 )
END SUBROUTINE cable_diag1

!=============================================================================!
!=============================================================================!

SUBROUTINE cable_diag_desc1( Nvars, filename, dimx, dimy, vname1 )

   integer, intent(in) :: Nvars,dimx,dimy 
   character(len=*), intent(in) :: filename, vname1
   integer, save :: gopenstatus = 1

     open(unit=713941,file=filename//'.dat', status="replace", &
          action="write", iostat=gopenstatus )
     
      if(gopenstatus==gok) then
            write (713941,*) 'Number of var(s): '
            write (713941,*) Nvars
            write (713941,*) 'Name of var(s): '
            write (713941,7139) vname1 
 7139       format(a)            
            write (713941,*) 'dimension of var(s) in x: '
            write (713941,*) dimx 
            write (713941,*) 'dimension of var(s) in y: '
            write (713941,*) dimy 
      else
         write (*,*), filename//'.dat',' Error: unable to write'
      endif
      
   close(713941)
  
END SUBROUTINE cable_diag_desc1


SUBROUTINE cable_diag_data1( Nvars, filename, dimx, timestep, kend, var1  )

   integer, intent(in) :: Nvars, dimx, timestep, kend
   real, intent(in), dimension(:) :: var1
   character(len=*), intent(in) :: filename
   integer, save :: gopenstatus = 1

   if (timestep == 1)  then 
      open(unit=713942,file=filename//'.bin',status="unknown", &
           action="write", iostat=gopenstatus, form="unformatted", &
           position='append' )
   endif   
 
   if(gopenstatus==gok) then
         write (713942) var1
   else
      write (*,*) filename//'.bin',' NOT open for write. Error'
   endif

   if (timestep == kend) & 
      close(713942)

END SUBROUTINE cable_diag_data1


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

!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   use cable_common_module, only : ktau_gl, knode_gl

   character(len=*), intent(in) :: routname
      if(knode_gl==1) & 
         write(6,*) 'CABLE@  ', routname, ktau_gl

END SUBROUTINE cable_stat


END MODULE cable_diag_module



