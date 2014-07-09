module cable_ncdf_module
contains
subroutine predef_grid
   use netcdf
   implicit none
 
   ! This is the name of the data file we will read. 
   character (len = *), parameter :: FILE_NAME = "LAI.nc"
   character (len = *), parameter :: VARNAME = "field1392"
 
   ! We are reading 3D data, a predefined grid. 
   integer, parameter :: nlon = 192, nlat = 145, npseudo=12 !(pseudo,lat,lon)-(nz,ny,nx)
   real:: data_in(nlon,nlat,npseudo)
   !real:: data_in(npseudo,nlat,nlon)
 
   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, varid
 
   ! Loop indexes, and error handling.
   integer :: x, y, z 
   
   integer :: tests
   integer :: latDimId, lonDimId, pseudoDimId
   integer :: IQnlat,IQnlon,IQnpseudo 
   integer :: numdims
   integer :: VarNumDimIds(3)

    
   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
   print *, "File opened on ncid ", ncid 
   ! Get the varid of the data variable, based on its name.
   call check( nf90_inq_varid(ncid, VARNAME, varid) )
   print *, "and varid is: ", varid 

   tests = nf90_inq_dimid(ncid, "latitude", latDimID)
   tests = nf90_inq_dimid(ncid, "longitude", lonDimID)
   tests = nf90_inq_dimid(ncid, "pseudo", pseudoDimID)
   print *, "latdimiD : ", latDimID
   print *, "lonDimId : ", lonDimID
   print *, "pseudoDimId : ", pseudoDimID


   tests = nf90_inquire_dimension(ncid, latDimID, len = IQnlat) 
   tests = nf90_inquire_dimension(ncid, lonDimID, len = IQnlon) 
   tests = nf90_inquire_dimension(ncid, pseudoDimID, len = IQnpseudo) 
   print *, "latdimS  : ", IQnlat 
   print *, "lonDimS  : ", IQnlon 
   print *, "pseudoDimS  : ", IQnpseudo 

   tests = nf90_inquire_variable(ncid, VarId, ndims = numDims)
   print *, "and ToTdim : ", numDims 

   tests = nf90_inquire_variable(ncid, VarId, dimids = VarNumDimIds(:numDims))
   print *, "and shape: ", VarNumDimIds 
                                              
   ! Read the data.
   call check( nf90_get_var(ncid, varid, data_in) )
 
   ! Check the data.
   do x = 1, nlon
      do y = 1, nlat
         do z = 1, nPseudo
         !if (data_in(z, y, x) /= (x - 1) * NY + (y - 1)) then
            print *, "data_in(", z, ", ", y, ", ", x, ") = ", data_in(x, y, z)
            !stop "Stopped"
         !end if
         end do
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
