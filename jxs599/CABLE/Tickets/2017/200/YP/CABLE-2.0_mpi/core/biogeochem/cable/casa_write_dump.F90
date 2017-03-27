SUBROUTINE write_casa_dump( ncfile, casamet, casaflux, phen, climate, n_call, kend )
      USE netcdf
      USE cable_def_types_mod
      use casadimension, only : mdyear, mplant
      USE casavariable
#     ifndef UM_BUILD
  USE cable_diag_module,     ONLY : def_dims, def_vars, def_var_atts, &
       put_var_ncr1, put_var_ncr2,       &
       put_var_ncr3, stderr_nc
#     endif

      use cable_ncdf_module, only : def_dims, def_vars, def_var_atts, &
                                   put_var_nc, stderr_nc
      USE cable_io_vars_module, only : patch

  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
    n_call, &         ! this timestep #
    kend              ! final timestep of run

   !var (type) to write
      TYPE (casa_met),            INTENT(INOUT) :: casamet


      !number of instances. dummied here and so=1
      !integer :: inst =1

      !netcdf IDs/ names
      character(len=*),intent(in):: ncfile
      integer, parameter :: num_vars=7
      integer, parameter :: num_dims=4

      integer, dimension(ms)     :: &
         soil

      integer, dimension(mvtype) :: &
         tile

      integer, dimension(mplant) :: &
         plant

      integer, save :: ncid       ! netcdf file ID

      !vars
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/  "latitude     ", &
                           "longitude    ", &
                            "casamet_tairk", &
                            "tsoil        ", &
                            "moist        ", &
                            "cgpp         ", &
                            "crmplant     " /)

      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb

      !dims
      character(len=*), dimension(num_dims), parameter :: &
            dim_name =  (/ "tile", &
                           "soil", &
                           "plant",&
                           "time" /)

      integer, dimension(num_dims)  :: &
            dimID   ! (1) mp, (2) ms, (3) mplant (4) time

      integer, dimension(num_dims)  :: &
            !x,y generally lat/lon BUT for single site = 1,1
            dim_len
      !local only
     integer :: ncok      !ncdf return status
     integer :: i,j


      real(r_2), dimension(1:mp) :: &
         cgpp,   &
         tairk

      real(r_2), dimension(1:mp,1:ms) :: &
         tsoil, &
         moist

      real(r_2), dimension(1:mp,1:mplant) :: &
         crmplant

      ! END header
!      tairk = 0
!      cgpp  = 0
!      tsoil = 0
!      moist = 0
!      crmplant = 0
!      write(89,*) 'tsoil,gpp',casamet%Tsoilspin_1(1,:),casamet%cgppspin(1,:),mplant,ms
!      write(89,*) 'creating dump file'
!      n_call = 1
!      kend   = mdyear


      print *, 'yp wang: calling ncdf_dump'
!      print *, 'latitude= ', patch(:)%latitude
!      print *, 'longitude= ', patch(:)%longitude
      print *, 'filename= ', ncfile
      print *, 'constants= ', mp,ms,mplant,mdyear,ncid

      dim_len(1) = mp
      dim_len(2) = ms
      dim_len(3) = mplant
      dim_len(4) = mdyear
         ! create netCDF dataset: enter define mode
      ncok = nf90_create(path = ncfile, cmode = nf90_clobber, ncid = ncid)

!      print *, 'ncok =', ncok

!      ncok = nf90_create(path = 'dump_casamet.nc', cmode = nf90_noclobber, ncid = ncid)
      if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile)
!      if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', 'dump_casamet.nc')

!      print *, 'here 1' ,ncid

      ! define dimensions: from name and length
!      write(89,*) 'defining dims'
      call def_dims(num_dims, ncid, dimID, dim_len, dim_name )

!      print *, 'here 2',varID,num_dims
      ! define variables: from name, type, dims
!      write(89,*) 'defining vars'
      call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )

!      print *, 'here 3',varID,num_vars
      ! define variable attributes
!      write(89,*) 'defining attribution'
      call def_var_atts(ncfile, ncid, varID )
!      call def_var_atts('dump_casamet.nc', ncid, varID )

!      print *, 'here 4', varID
      ncok = nf90_enddef(ncid)

!      print *, 'here 5', var_name(1), size(patch(:)%latitude)
!      write(89,*) 'writing latitude'
      call put_var_nc(ncid, var_name(1), patch(:)%latitude )

!      print *, 'here 6',var_name(2) , size(patch(:)%longitude)
!      write(89,*) 'writing longitude'
      call put_var_nc(ncid, var_name(2), patch(:)%longitude )

      write(*,901)  mdyear 
901   format(' yp wang at ncdf_dump', I6)
!      write(*,*) casamet%cgppspin(10,:)

      do i=1,mdyear
         tairk(:)      = casamet%Tairkspin(:,i)
         tsoil(:,1)    = casamet%Tsoilspin_1(:,i)
!      tairk           = casamet%Tsoilspin_1
         tsoil(:,2)    = casamet%Tsoilspin_2(:,i)
         tsoil(:,3)    = casamet%Tsoilspin_3(:,i)
         tsoil(:,4)    = casamet%Tsoilspin_4(:,i)
         tsoil(:,5)    = casamet%Tsoilspin_5(:,i)
         tsoil(:,6)    = casamet%Tsoilspin_6(:,i)
         moist(:,1)    = casamet%moistspin_1(:,i)
         moist(:,2)    = casamet%moistspin_2(:,i)
         moist(:,3)    = casamet%moistspin_3(:,i)
         moist(:,4)    = casamet%moistspin_4(:,i)
         moist(:,5)    = casamet%moistspin_5(:,i)
         moist(:,6)    = casamet%moistspin_6(:,i)
         cgpp (:)      = casamet%cgppspin   (:,i)
         crmplant(:,1) = casamet%crmplantspin_1(:,i)
         crmplant(:,2) = casamet%crmplantspin_2(:,i)
         crmplant(:,3) = casamet%crmplantspin_3(:,i)

         call put_var_nc(ncid, var_name(3), tairk, i,kend )
         call put_var_nc(ncid, var_name(4), tsoil, i,kend ,ms)
         call put_var_nc(ncid, var_name(5), moist, i,kend ,ms)
         call put_var_nc(ncid, var_name(6), cgpp , i,kend )
         call put_var_nc(ncid, var_name(7), crmplant, i,kend,mplant )
     end do

  IF (n_call == kend ) &
       ncok = nf90_close(ncid)            ! close: save new netCDF dataset

#endif
END SUBROUTINE write_casa_dump


