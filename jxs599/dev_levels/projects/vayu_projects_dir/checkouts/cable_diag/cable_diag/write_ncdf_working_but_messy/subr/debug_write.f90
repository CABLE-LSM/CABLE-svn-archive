

module debug_write_mod
   implicit none
   public

   interface put_var
      module procedure fput_var, iput_var, fput_var4
   end interface

   contains

   subroutine write_ncdf_template( ncfile, nlat, nlon, ntile, ntime, nvar, lat, lon, tile, timestep, newvar )
      use netcdf

      implicit none  

      character(len=*), intent(in) :: ncfile 
      integer, intent(in) :: nlat, nlon, ntile, ntime, nvar
      real, dimension(:,:,:,:), intent(in) :: newvar 
      real, dimension(:), intent(in) :: lat, lon
      integer, dimension(:), intent(in) ::  tile, timestep 

      integer :: ncid       ! netcdf file ID

      !dims
      integer, parameter :: num_dims=4 
      integer, dimension(num_dims)  :: dimID   
      
      character(len=*), dimension(num_dims), parameter :: &
         dim_name =  (/  "lon ", &
                         "lat ", &
                         "tile", &
                         "time" /)

      integer, dimension(num_dims) :: dim_len 
      
      !vars 
      integer, parameter :: num_vars=5 
      integer, dimension(num_vars) :: varID 
      
      character(len=*), dimension(num_vars), parameter :: &
            var_name =  (/  "lon     ", &
                            "lat     ", &
                            "tile     ", &
                            "time     ", &
                            "var" /)
      
      !local only
      integer :: ncok      !ncdf return status
      

         dim_len =  (/  nlon, &
                        nlat, &
                        ntile, &
                        ntime /)

         ! create netCDF dataset: enter define mode
         ncok = NF90_CREATE(path = ncfile, cmode = nf90_noclobber, ncid = ncid)
            if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile) 
      
            call def_dims(num_dims, ncid, dimID, dim_len, dim_name )
      
            ! define variables: from name, type, dims
            !call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )
       
            ncok = NF90_DEF_VAR(ncid, trim(var_name(1)), nf90_float, &
                            dimID(1), varID(1))
               if (ncok /= nf90_noerr) call stderr_nc('ncdf def var', trim(var_name(1))) 
            
            ncok = NF90_DEF_VAR(ncid, trim(var_name(2)), nf90_float, &
                            dimID(2), varID(2))
               if (ncok /= nf90_noerr) call stderr_nc('ncdf def var', trim(var_name(2))) 
            
            ncok = NF90_DEF_VAR(ncid, trim(var_name(3)), nf90_int, &
                            dimID(3), varID(3))
               if (ncok /= nf90_noerr) call stderr_nc('ncdf def var', trim(var_name(3))) 
            
            ncok = NF90_DEF_VAR(ncid, trim(var_name(4)), nf90_int, &
                            dimID(4), varID(4))
               if (ncok /= nf90_noerr) call stderr_nc('ncdf def var', trim(var_name(4))) 
            
            ncok = NF90_DEF_VAR(ncid, trim(var_name(5)), nf90_float, &
                            !(/dimID(1),dimID(2),dimID(4),dimID(3) /), varID(5))
                            dimID, varID(5))
               if (ncok /= nf90_noerr) call stderr_nc('ncdf def var', trim(var_name(5))) 
            
         ncok = nf90_enddef(ncid)
         ncok = nf90_close(ncid)            ! close: save new netCDF dataset

      call put_var(ncfile, var_name(1), lon)
      call put_var(ncfile, var_name(2), lat)
      call put_var(ncfile, var_name(3), tile)
      call put_var(ncfile, var_name(4), timestep )
      call put_var(ncfile, var_name(5), newvar)

      return 
   end subroutine write_ncdf_template


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
            if (ncok /= nf90_noerr ) call stderr_nc('def air_in dim ', dim_name(j))      
      enddo

      return
   end subroutine def_dims

   subroutine iput_var(ncfile_in, var_name, var)
      use netcdf
      implicit none
      character(len=*), intent(in) :: ncfile_in, var_name
      integer, dimension(:),intent(in) :: var
      integer :: ncid, ncok, varID
       
      ncok = NF90_OPEN(ncfile_in, nf90_write, ncid)           
         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile_in)      
         ncok = NF90_INQ_VARID(ncid, var_name, varId )
            if (ncok /= nf90_noerr ) call stderr_nc('inquire var ', var_name)      
         ncok = NF90_PUT_VAR(ncid, varId, var )
         !ncok = NF90_PUT_VAR(ncid, varId, var, start=(/1,1,1,1/) )
            if (ncok /= nf90_noerr ) call stderr_nc('putting var ', var_name)      
      ncok = NF90_CLOSE(ncid)            
         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile_in)      
     
      return
   end subroutine iput_var
   
   
   
   subroutine fput_var(ncfile_in, var_name, var)
      use netcdf
      implicit none
      character(len=*), intent(in) :: ncfile_in, var_name
      real, dimension(:),intent(in) :: var
      integer :: ncid, ncok, varID
       
      ncok = NF90_OPEN(ncfile_in, nf90_write, ncid)           
         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile_in)      
         ncok = NF90_INQ_VARID(ncid, var_name, varId )
            if (ncok /= nf90_noerr ) call stderr_nc('inquire var ', var_name)      
         ncok = NF90_PUT_VAR(ncid, varId, var )
            if (ncok /= nf90_noerr ) call stderr_nc('putting var ', var_name)      
      ncok = NF90_CLOSE(ncid)            
         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile_in)      
     
      return
   end subroutine fput_var



   
   
   subroutine fput_var4(ncfile_in, var_name, var)
      use netcdf
      implicit none
      character(len=*), intent(in) :: ncfile_in, var_name
      real, dimension(:,:,:,:),intent(in) :: var
      integer :: ncid, ncok, varID
       
      ncok = NF90_OPEN(ncfile_in, nf90_write, ncid)           
         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile_in)      
         ncok = NF90_INQ_VARID(ncid, var_name, varId )
            if (ncok /= nf90_noerr ) call stderr_nc('inquire var ', var_name)      
         ncok = NF90_PUT_VAR(ncid, varId, var )
            if (ncok /= nf90_noerr ) call stderr_nc('putting var ', var_name)      
      ncok = NF90_CLOSE(ncid)            
         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile_in)      
     
      return
   end subroutine fput_var4





   subroutine stderr_nc(message, var)      
      character(len=*), intent(in) :: message, var
      character(len=7) :: err_mess
         err_mess = 'ERROR:'
         print *, (err_mess//message), var
      stop
   end subroutine stderr_nc      



 
!   subroutine def_vars(nv, ncid,  xtype, dimID, var_name,varID )
!      use netcdf
!      implicit none 
!      integer, intent(in) :: nv, ncid, xtype 
!      integer, dimension(:), intent(in) :: dimID
!      integer, dimension(:), intent(inout) :: varID
!      character(len=*), dimension(:), intent(in) :: var_name
!      integer :: j, ncok
!      
!      do j=1, nv
!         ncok = NF90_DEF_VAR(ncid, trim(var_name(j)), xtype, &
!                            dimID, varID(j))
!            if (ncok /= nf90_noerr ) call stderr_nc('def var ', var_name(j))      
!      enddo
!
!      return
!   end subroutine def_vars
!
!
!
!!   subroutine def_var_atts(ncfile_in, ncid, varID )
!!     use netcdf
!!     implicit none   
!!     character(len=*), intent(in) :: ncfile_in
!!     integer, intent(in):: ncid       ! netcdf file ID
!!     integer, dimension(:), intent(in) :: varID ! (1) ~ tvair, (2) ~ pmb 
!!     integer :: j, ncok
!!     character(len=10) dummy
!!
!!     write(dummy,11) varID(1)
!!  11 format(i2)   
!!     ncok = NF90_PUT_ATT(ncid, nf90_global, "Title", "Forcing for define_air subroutine")
!!         if (ncok /= nf90_noerr ) call stderr_nc('def att ', ncfile_in)      
!!     ncok = NF90_PUT_ATT(ncid, varID(1), "longname", "air temperature within canopy")
!!         if (ncok /= nf90_noerr ) call stderr_nc('def att ', dummy)      
!!     ncok = NF90_PUT_ATT(ncid, varID(1), "units", "Kelvin")
!!         if (ncok /= nf90_noerr ) call stderr_nc('def att ', dummy)      
!!     
!!     write(dummy,11) varID(2)
!!     
!!
!!     return
!!   end subroutine def_var_atts
!! 



!!   subroutine get_var(ncfile_in, var_name, var, n_call)
!!      use netcdf
!!      implicit none
!!      character(len=*), intent(in) :: ncfile_in, var_name
!!      real, dimension(:),intent(out) :: var
!!      integer, intent(in) :: n_call
!!      integer :: ncid, ncok, varID
!!       
!!      ncok = NF90_OPEN(ncfile_in, nf90_nowrite, ncid)           
!!         if (ncok /= nf90_noerr ) call stderr_nc('re-opening ', ncfile_in)      
!!         ncok = NF90_INQ_VARID(ncid, var_name, varId )
!!            if (ncok /= nf90_noerr ) call stderr_nc('inquire var ', var_name)      
!!         ncok = NF90_GET_VAR(ncid, varId, var, start=(/1,1,n_call/) )
!!            if (ncok /= nf90_noerr ) call stderr_nc('getting var ', var_name)      
!!      ncok = NF90_CLOSE(ncid)            
!!         if (ncok /= nf90_noerr ) call stderr_nc('closing ', ncfile_in)      
!!     
!!      return
!!   end subroutine get_var
!! 

end module debug_write_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!!   !this subr is reasonably generic following some decs EXCEPT for call to _put_var**
!!   !subroutine air_in_dump(i, n_call, instances, ncfile,inst)
!!   subroutine write_ncdf_file( ncfile, newdata)
!!      use netcdf
!!
!!      implicit none  
!!      real, intent(inout), dimension(:) :: newdata
!!
!!      !netcdf IDs/ names 
!!      character(len=*), intent(in) :: ncfile 
!!      character(len=*) :: wfile 
!!
!!      integer, parameter :: num_vars=3 
!!      integer, parameter :: num_dims=3 
!!      integer :: ncid       ! netcdf file ID
!!      
!!      !vars 
!!      character(len=*), dimension(num_vars), parameter :: &
!!            var_name =  (/  "lat     ", &
!!                            "lon     ", &
!!                            "can_flux" /)
!!      integer, dimension(num_vars) :: varID 
!!      
!!      integer, dimension(num_dims)  :: &
!!            dimID   ! (1) x, (2) y, (3) time
!!      
!!      !local only
!!      integer :: ncok      !ncdf return status
!!      
!!      wfile = trim(ncfile//".nc")
!!
!!         ! create netCDF dataset: enter define mode
!!         ncok = NF90_CREATE(path = wfile, cmode = nf90_noclobber, ncid = ncid)
!!            if (ncok /= nf90_noerr) call stderr_nc('ncdf creating ', ncfile) 
!!      
!!            ! define dimensions: from name and lengthCopy of UM in CABLE repository to be compatible with CABLE-1.9.
!!            call init_dims(ncid, dimID, num_dims)
!!     
!!            ! define variables: from name, type, dims
!!            call def_vars(num_vars, ncid,  nf90_float, dimID, var_name, varID )
!!       
!!            ncok = NF90_DEF_VAR(ncid, trim(var_name(1)), nf90_float, &
!!                            dimID(1), varID(1))
!!            
!!            ncok = NF90_DEF_VAR(ncid, trim(var_name(2)), nf90_float, &
!!                            dimID(2), varID(2))
!!            
!!            ncok = NF90_DEF_VAR(ncid, trim(var_name(3)), nf90_float, &
!!                            dimID, varID(3))
!!            
!!            
!!            ! define variable attributes
!!            call def_var_atts(wfile, ncid, varID )
!!               
!!            ncok = nf90_enddef(ncid) 
!!         
!!         ncok = nf90_close(ncid)            ! close: save new netCDF dataset
!!      
!!      call put_var(wfile, var_name(1), lat)
!!      call put_var(wfile, var_name(2), lon)
!!      call put_var(wfile, var_name(3), newdata)
!!
!!      return 
!!   end subroutine write_ncdf_file
!!
!!


