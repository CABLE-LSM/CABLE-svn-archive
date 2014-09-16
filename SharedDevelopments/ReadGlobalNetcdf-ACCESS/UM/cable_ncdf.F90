!module cable_ncdf_module
!contains
subroutine predef_grid(latitude, longitude, node_gl, rows, row_length, mp,       &
                     npseudo,mydata)
   use netcdf
   implicit none
 
   ! we are passing these as they are declared i8 and r8
   integer*8 :: node_gl, rows, row_length, mp, npseudo
   REAL*8, dimension(mp) :: latitude, longitude 
   real*8, dimension(mp, npseudo) :: mydata 

   ! This is the name of the data file we will read. 
   character (len = *), parameter :: fbase_name = "LAI"
   character (len = *), parameter :: fext_name = ".nc"
   character (len =10) :: file_name = ".nc"
   character (len = *), parameter :: VARNAME = "field1392"
 
   ! We are reading 3D data, a predefined grid. 
   integer, parameter :: nlon = 192, nlat = 145  !(pseudo,lat,lon)-(nz,ny,nx)
   real:: lat_in(nlat,npseudo)
   real:: lon_in(nlon,npseudo)
   real:: data_in(nlon,nlat,npseudo)
   integer, dimension(mp) :: x_in, y_in
 
   ! This will be the netCDF ID for the file and data variable.
   integer :: ncid, varid, latid, lonid
 
   ! Loop indexes, and error handling.
   integer :: i,x, y, z 
   
   integer :: tests
   !integer :: latDimId, lonDimId, pseudoDimId
   !integer :: IQnlat,IQnlon,IQnpseudo 
   !integer :: numdims
   !integer :: VarNumDimIds(3)
   character(len=3) :: chnode
   character(len=6) :: frmat
   
   if(node_gl > 99) then
      frmat = "(I3.1)"
   elseif( node_gl > 9 .AND. node_gl < 100 ) then    
      frmat = "(I2.1)"
   else     
      frmat = "(I1.1)"
   endif   
   write(chnode,frmat) node_gl

   file_name = trim(fbase_name)//trim(chnode)//trim(fext_name) 

   ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
   ! the file.
   call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
   !print *, "File opened on ncid ", ncid 
   
   ! Get the varid of the data variable, based on its name.
   call check( nf90_inq_varid(ncid, VARNAME, varid) )
   call check( nf90_inq_varid(ncid, 'latitude', latid) )
   call check( nf90_inq_varid(ncid, 'longitude', lonid) )
   !print *, "and varid is: ", varid 

   !tests = nf90_inq_dimid(ncid, "latitude", latDimID)
   !tests = nf90_inq_dimid(ncid, "longitude", lonDimID)
   !tests = nf90_inq_dimid(ncid, "pseudo", pseudoDimID)
   !print *, "latdimiD : ", latDimID
   !print *, "lonDimId : ", lonDimID
   !print *, "pseudoDimId : ", pseudoDimID

   !tests = nf90_inquire_dimension(ncid, latDimID, len = IQnlat) 
   !tests = nf90_inquire_dimension(ncid, lonDimID, len = IQnlon) 
   !tests = nf90_inquire_dimension(ncid, pseudoDimID, len = IQnpseudo) 
   !print *, "latdimS  : ", IQnlat 
   !print *, "lonDimS  : ", IQnlon 
   !print *, "pseudoDimS  : ", IQnpseudo 

   !tests = nf90_inquire_variable(ncid, VarId, ndims = numDims)
   !print *, "and ToTdim : ", numDims 

   !tests = nf90_inquire_variable(ncid, VarId, dimids = VarNumDimIds(:numDims))
   !print *, "and shape: ", VarNumDimIds 
   ! Read the data.
   call check( nf90_get_var(ncid, latid, lat_in) )
   call check( nf90_get_var(ncid, lonid, lon_in) )
   call check( nf90_get_var(ncid, varid, data_in) )
   !call check( nf90_get_var(ncid, varid, mydata, & 
   !            start=(/latitude(1),longitude(1),1/),count=(/rows,row_length,12/) ) )
 
   ! Check the data.
   !lon, lat, t
   !row_length, rows, t
   x_in = 1; y_in =1
   z = 1
   do x = 1, nlon
      do y = 1, nlat
         do i = 1, mp
            if( ( abs(lon_in(x,1) - longitude(i) ) < 0.001 ) .And. & 
                ( abs(lat_in(y,1) - latitude(i)  ) < 0.001 ) ) then 
               x_in(i) = x
               !print *, "x_in ", x_in(i)
               y_in(i) = y
               if( i > mp) print *, "Count has reached mp"
            endif
         end do
      end do
   end do

   do z = 1, nPseudo
      do i = 1, mp
         mydata(i,z) =  DBLE( data_in(x_in(i), y_in(i), z) )  
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

!end module cable_ncdf_module
