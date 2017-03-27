
   subroutine read_casa_dump(ncfile, casamet, casaflux, ktau, kend)
      use netcdf
      use cable_def_types_mod,   only : r_2,ms
      use casadimension,         only : mplant,mdyear
      USE casa_cnp_module
      use cable_ncdf_module,     only : get_var_nc, stderr_nc

      TYPE (casa_flux), intent(inout) :: casaflux
      TYPE (casa_met), intent(inout)  :: casamet
      integer, intent(in) :: kend, ktau

      !netcdf IDs/ names
      character(len=*)  ncfile
      integer, parameter :: num_vars=5
      integer, parameter :: num_dims=4
      integer:: ncid       ! netcdf file ID

      !vars
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/  "casamet_tairk", &
                            "tsoil        ", &
                            "moist        ", &
                            "cgpp         ", &
                            "crmplant     " /)

      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb

      real(r_2), dimension(mp) :: &
         tairk,  &
         cgpp

      real(r_2), dimension(mp,ms) :: &
         tsoil, &
         moist

      real(r_2), dimension(mp,mplant) :: &
         crmplant

!      write(89,*)'opening file'
      ncok = NF90_OPEN(ncfile, nf90_nowrite, ncid)
         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile)

      do idoy=1,mdyear

!         write(89,*)'get tairk'
         call get_var_nc(ncid, var_name(1), tairk   , idoy, kend )
!         write(89,*)'get tsoil'
         call get_var_nc(ncid, var_name(2), tsoil   , idoy, kend ,ms)
!         write(89,*)'get moist'
         call get_var_nc(ncid, var_name(3), moist   , idoy, kend ,ms)
!         write(89,*)'get cgpp'
         call get_var_nc(ncid, var_name(4), cgpp    , idoy, kend )
!         write(89,*)'get crmplant'
         call get_var_nc(ncid, var_name(5), crmplant, idoy, kend ,mplant)


         casamet%Tairkspin(:,idoy) = tairk
         casamet%cgppspin (:,idoy) = cgpp
         casamet%crmplantspin_1(:,idoy) = crmplant(:,1)
         casamet%crmplantspin_2(:,idoy) = crmplant(:,2)
         casamet%crmplantspin_3(:,idoy) = crmplant(:,3)
         casamet%Tsoilspin_1(:,idoy)    = tsoil(:,1)
         casamet%Tsoilspin_2(:,idoy)    = tsoil(:,2)
         casamet%Tsoilspin_3(:,idoy)    = tsoil(:,3)
         casamet%Tsoilspin_4(:,idoy)    = tsoil(:,4)
         casamet%Tsoilspin_5(:,idoy)    = tsoil(:,5)
         casamet%Tsoilspin_6(:,idoy)    = tsoil(:,6)
         casamet%moistspin_1(:,idoy)    = moist(:,1)
         casamet%moistspin_2(:,idoy)    = moist(:,2)
         casamet%moistspin_3(:,idoy)    = moist(:,3)
         casamet%moistspin_4(:,idoy)    = moist(:,4)
         casamet%moistspin_5(:,idoy)    = moist(:,5)
         casamet%moistspin_6(:,idoy)    = moist(:,6)

      end do

      ncok = NF90_CLOSE(ncid)
         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile)


   end subroutine read_casa_dump


