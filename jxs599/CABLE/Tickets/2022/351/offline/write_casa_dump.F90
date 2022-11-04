MODULE write_dump_casa_mod

IMPLICIT NONE

#define  UM_BUILD YES

#ifndef UM_BUILD

CONTAINS

SUBROUTINE write_casa_dump( ncfile, casamet, casaflux, phen, climate, n_call, kend )
  USE netcdf
  USE cable_def_types_mod,   ONLY : r_2,ms,mp, climate_type
  USE cable_common_module,   ONLY : kend_gl
  USE casa_ncdf_module,     ONLY : def_dims, def_vars, def_var_atts, &
       put_var_ncr1, put_var_ncr2,       &
       put_var_ncr3, stderr_nc
  USE casavariable,          ONLY : CASA_MET, CASA_FLUX
  USE casadimension,         ONLY : mplant
  USE phenvariable
  USE cable_common_module,  ONLY:  CABLE_USER

  IMPLICIT NONE

  INTEGER, INTENT(in) :: &
       n_call, &         ! this timestep #
       kend              ! final timestep of run

  TYPE (casa_flux),             INTENT(IN) :: casaflux
  TYPE (casa_met),              INTENT(IN) :: casamet
  TYPE (phen_variable),         INTENT(IN) :: phen
  TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables

  !number of instances. dummied here and so=1
  !integer :: inst =1

  !netcdf IDs/ names
  CHARACTER(len=*)   :: ncfile
  INTEGER            :: num_vars
  INTEGER, PARAMETER :: num_dims=3
  INTEGER, SAVE :: ncid       ! netcdf file ID

  !vars
  CHARACTER, DIMENSION(:), POINTER :: var_name*15


  INTEGER, DIMENSION(:), POINTER :: varID ! (1) tvair, (2) pmb

  !dims
  CHARACTER(len=*), DIMENSION(num_dims), PARAMETER :: &
       dim_name =  (/ "pnt ", &
       "soil", &
       "time" /)

  INTEGER, PARAMETER :: soil_dim = 6

  INTEGER, DIMENSION(soil_dim), PARAMETER  :: soil = (/ 1,2,3,4,5,6 /)

  INTEGER, DIMENSION(num_dims)  :: &
       dimID   ! (1) x, (2) y, (3) time

  INTEGER, DIMENSION(num_dims)  :: &
                                !x,y generally lat/lon BUT for single site = 1,1
       dim_len = (/-1,soil_dim,-1/)  ! (1) mp, (2) soil, (3) time [re-set]



  !local only
  INTEGER :: ncok      !ncdf return status

  ! END header
  dim_len(1)        = mp
  dim_len(num_dims) = NF90_unlimited

  num_vars = 14

  !Add extra mtemp variable when running with climate
  IF (cable_user%CALL_climate) THEN
    num_vars=num_vars+1
  ENDIF

  allocate(var_name(num_vars))
  allocate(varID(num_vars))

  var_name =  (/"lat          ", &
                "lon          ", &
                "casamet_tairk", &
                "tsoil        ", &
                "moist        ", &
                "cgpp         ", &
                "crmplant     ", &
                "phenphase    ", &
                "phendoyphase1", &
                "phendoyphase2", &
                "phendoyphase3", &
                "phendoyphase4", &
                "Ndep         ", &
                "Pdep         "/)

  !Add extra mtemp variable when running with climate
  IF (cable_user%CALL_climate) THEN
    var_name(num_vars)="mtemp"
  ENDIF

  IF (n_call == 1) THEN

     ! create netCDF dataset: enter define mode
     ncok = nf90_create(path = TRIM(ncfile), cmode = nf90_clobber, ncid = ncid)
     IF (ncok /= nf90_noerr) CALL stderr_nc(ncok,'ncdf creating ', ncfile)

     !ncok = nf90_redef(ncid)
     !if (ncok /= nf90_noerr) call stderr_nc(ncok,'enter def mode', ncfile)

     ! define dimensions: from name and length
     CALL def_dims(num_dims, ncid, dimID, dim_len, dim_name )

     ! define variables: from name, type, dims
     CALL def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )

     ! define variable attributes
     !CLN LATER!             CALL def_var_atts( ncfile, ncid, varID )

     ncok = nf90_enddef(ncid)
     if (ncok /= nf90_noerr) call stderr_nc(ncok,'end def mode', ncfile)


     CALL put_var_ncr1(ncid, var_name(1), REAL(casamet%lat)  )
     CALL put_var_ncr1(ncid, var_name(2), REAL(casamet%lon)  )


  ENDIF


  CALL put_var_ncr2(ncid, var_name(3), casamet%tairk    ,n_call )
  CALL put_var_ncr3(ncid, var_name(4), casamet%tsoil    ,n_call, ms )
  CALL put_var_ncr3(ncid, var_name(5), casamet%moist    ,n_call, ms )
  CALL put_var_ncr2(ncid, var_name(6), casaflux%cgpp    ,n_call )
  CALL put_var_ncr3(ncid, var_name(7), casaflux%crmplant,n_call, mplant )
  CALL put_var_ncr2(ncid, var_name(8), real(phen%phase , r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(9), real(phen%doyphase(:,1), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(10), real(phen%doyphase(:,2), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(11), real(phen%doyphase(:,3), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(12), real(phen%doyphase(:,4), r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(13), real(casaflux%Nmindep,r_2)    ,n_call )
  CALL put_var_ncr2(ncid, var_name(14), real(casaflux%Pdep,r_2)    ,n_call )
  if (cable_user%CALL_climate) then
     CALL put_var_ncr2(ncid, var_name(15), real(climate%mtemp_max,r_2)    ,n_call )
  endif

  deallocate(var_name)

  IF (n_call == kend ) &
       ncok = nf90_close(ncid)            ! close: save new netCDF dataset

#endif
END SUBROUTINE write_casa_dump

END MODULE write_dump_casa_mod
