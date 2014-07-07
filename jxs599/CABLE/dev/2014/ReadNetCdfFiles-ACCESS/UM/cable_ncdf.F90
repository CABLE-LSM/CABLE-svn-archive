module cable_ncdf_module
contains
subroutine predef_grid
   use netcdf
   implicit none
 
   ! This is the name of the data file we will read. 
   character (len = *), parameter :: FILE_NAME = "LAI.nc"
 
   ! We are reading 3D data, a predefined grid. 
   integer, parameter :: nx = 192, ny = 145, nz=600
   integer :: data_in(nz,ny,nx)
 
   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, varid
 
   ! Loop indexes, and error handling.
   integer :: x, y, z 
 
   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
 
   ! Get the varid of the data variable, based on its name.
   call check( nf90_inq_varid(ncid, "", varid) )
 
   ! Read the data.
   call check( nf90_get_var(ncid, varid, data_in) )
 
   ! Check the data.
   do x = 1, NX
      do y = 1, NY
         if (data_in(y, x) /= (x - 1) * NY + (y - 1)) then
            print *, "data_in(", y, ", ", x, ") = ", data_in(y, x)
            stop "Stopped"
         end if
      end do
   end do
 
   ! Close the file, freeing all resources.
   call check( nf90_close(ncid) )
 
   print *,"*** SUCCESS reading example file ", FILE_NAME, "! "
   
   contains
   subroutine check(status)
     integer, intent ( in) :: status
     
     if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
     end if
   end subroutine check  

end subroutine predef_grid

end module cable_ncdf_module
