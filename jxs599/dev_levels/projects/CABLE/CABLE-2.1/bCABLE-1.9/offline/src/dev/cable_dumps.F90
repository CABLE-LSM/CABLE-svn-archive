
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++ Specific read/write routines to subr +++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module cable_dumps_module
   IMPLICIT NONE
   PUBLIC 

   !interface 
   !   module procedure  
   !end interface 
   
CONTAINS

   !this subr is reasonably generic following some decs EXCEPT for call to _put_var**
   subroutine air_in_dump(i, n_call, instances, ncfile,inst)
      use netcdf
      use cable_data_module, only : air_in
      use cable_diag_module, only : def_dims, def_vars, def_var_atts, & 
                                    init_dims, put_var, stderr_nc

      implicit none  
      !var (type) to write 
      type (air_in), intent(in) :: i

      integer, intent(in) :: & 
            n_call, & ! this timestep # 
            instances, inst

      !netcdf IDs/ names 
      character(len=*), intent(in) :: ncfile 
      integer, parameter :: num_vars=2 
      integer, parameter :: num_dims=4 
      integer :: ncid       ! netcdf file ID
      
      !vars 
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/ "met_tvair", &
                           "met_pmb  " /)
      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb 
      
      integer, dimension(num_dims)  :: &
            dimID   !(1) instance, (1) x, (2) y, (3) time
      
      !local only
      integer :: ncok      !ncdf return status
      

      !on first call & first instance set up netcdf file
      if (n_call == 1 .AND. inst == 1 ) then
         ! create netCDF dataset: enter define mode
         ncok = NF90_CREATE(path = ncfile, cmode = nf90_noclobber, ncid = ncid)
            if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile) 
      
            ! define dimensions: from name and length
            call init_dims(ncid, dimID, num_dims, instances)
     
            ! define variables: from name, type, dims
            call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )
       
            ! define variable attributes
            call def_var_atts(ncfile, ncid, varID )
               
            ncok = nf90_enddef(ncid) 
         
         ncok = nf90_close(ncid)            ! close: save new netCDF dataset
      endif 
      
      call put_var(ncfile, var_name(1), i%met_tvair, n_call, inst)
      call put_var(ncfile, var_name(2), i%met_pmb, n_call, inst)

      return 
   end subroutine air_in_dump



   subroutine air_in_read(i, n_call)
      use netcdf
      use cable_data_module, only : air_in
      use cable_diag_module, only : get_var
      implicit none   

      !var (type) to write 
      type (air_in), intent(out) :: i

      integer, intent(in) :: n_call

      !netcdf IDs/ names 
      character(len=*), parameter :: ncfile = "air_in.nc"
      integer, parameter :: num_vars=2 
      integer, parameter :: num_dims=3 
      integer:: ncid       ! netcdf file ID
      
      !vars 
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/ "met_tvair", &
                           "met_pmb  " /)
      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb 

      call get_var(ncfile, var_name(1), i%met_tvair, n_call)
      call get_var(ncfile, var_name(2), i%met_pmb, n_call)
      
      return 
   end subroutine air_in_read



   subroutine air_out_dump(o, n_call)
      use netcdf
      use cable_common_module, only : kend_gl
      use cable_data_module, only : air_out
      use cable_diag_module, only : def_dims, def_vars, def_var_atts, & 
                                    put_var, stderr_nc

      implicit none  
      !var (type) to write 
      type (air_out), intent(in) :: o
integer :: inst =1
      ! this timestep # 
      integer, intent(in) :: n_call   
      
      !netcdf IDs/ names 
      character(len=*), parameter :: ncfile = "air_out.nc"
      integer, parameter :: num_vars=7 
      integer, parameter :: num_dims=3 
      integer:: ncid       ! netcdf file ID
      
      !vars 
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/    "air_cmolar", & 
                              "air_rho   ", & 
                              "air_rlam  ", &
                              "air_epsi  ", & 
                              "air_visc  ", & 
                              "air_psyc  ", & 
                              "air_dsatdk" /)
      integer, dimension(num_vars) :: varID ! (1) tvair, (2) pmb 
      
      !dims 
      character(len=*), dimension(num_dims), parameter :: & 
            dim_name =  (/ "lat ", &
                           "lon ", &
                           "time" /)
      
      integer, dimension(num_dims)  :: &
            dimID   ! (1) x, (2) y, (3) time
      
      integer, dimension(num_dims)  :: &
            !x,y generally lat/lon BUT for single site = 1,1       
            dim_len = (/1,1,-1/)  ! (1) x, (2) y, (3) time [re-set] 
      
      
      !local only
      integer :: ncok      !ncdf return status
      
      dim_len(3) = kend_gl 

      if (n_call == 1) then
         ! create netCDF dataset: enter define mode
         ncok = nf90_create(path = ncfile, cmode = nf90_noclobber, ncid = ncid)
            if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile) 
      
            ! define dimensions: from name and length
            call def_dims(num_dims, ncid, dimID, dim_len, dim_name )
     
            ! define variables: from name, type, dims
            call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )
       
            ! define variable attributes
            call def_var_atts(ncfile, ncid, varID )
               
            ncok = nf90_enddef(ncid) 
         
         ncok = nf90_close(ncid)            ! close: save new netCDF dataset
      endif 
      
      call put_var(ncfile, var_name(1), o%air_cmolar, n_call, inst)
      call put_var(ncfile, var_name(2), o%air_rho, n_call, inst)
      call put_var(ncfile, var_name(3), o%air_rlam, n_call, inst)
      call put_var(ncfile, var_name(4), o%air_epsi, n_call, inst)
      call put_var(ncfile, var_name(5), o%air_visc, n_call, inst)
      call put_var(ncfile, var_name(6), o%air_psyc, n_call, inst)
      call put_var(ncfile, var_name(7), o%air_dsatdk, n_call, inst)

      return 
   end subroutine air_out_dump




end module cable_dumps_module
