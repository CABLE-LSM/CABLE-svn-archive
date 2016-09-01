! jtk561 - subroutines to read the binary files, and write ALMA nc files
MODULE fileio_module
CONTAINS
!!!!!!!!!!!!!!! binary file reader!!!!!!!!!!!!!!!!!!
SUBROUTINE read_binary_file(filename,data_out)
USE GlobalDefs
INTEGER :: iunit, iostatus
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL, DIMENSION(:), INTENT(OUT) :: data_out
! the *4 in recl is required due to a machine dependant thing (the file size
! divide by dx*dy gives exactly 4
! not initializing iunit seems to create weird problems.....
iunit = 10
OPEN(iunit, file=TRIM(filename),form='unformatted',status='old',action='read',access='direct',recl=np*4,iostat=iostatus)
IF (iostatus .NE. 0) THEN
   PRINT *, 'Problem opening ', TRIM(filename),' ,iostatus is: ', iostatus
   CLOSE(iunit)
   STOP
END IF
READ(unit=iunit,iostat=iostatus,rec=1) data_out
IF (iostatus .NE. 0) THEN
   PRINT *, 'Problem reading ', TRIM(filename),' ,iostatus is: ', iostatus
   CLOSE(iunit)
   STOP
END IF
CLOSE(iunit)
END SUBROUTINE read_binary_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_all_bom_data(datein,all_data)
USE GlobalDefs
USE DateFunctions
!INTEGER, INTENT(IN) :: year, month, day
type(dmydate), INTENT(IN) :: datein 
REAL, DIMENSION(:,:), INTENT(OUT) :: all_data 
REAL, DIMENSION(:), ALLOCATABLE :: data_out
CHARACTER(LEN=2) :: day_str, mon_str
CHARACTER(LEN=4) :: year_str
CHARACTER(LEN=200) :: filein_tmax, filein_tmin, filein_rain, filein_rad, &
                      filein_vph09, filein_vph15

! pad zeros
IF (datein%Day < 10) THEN
      WRITE(day_str,'(A,I1.2)') '0', datein%Day
ELSE      
      WRITE(day_str,'(I2)') datein%Day
ENDIF
IF (datein%Month < 10) THEN
      WRITE(mon_str,'(A,I1.2)') '0', datein%Month
ELSE
      WRITE(mon_str,'(I2)') datein%Month
ENDIF
WRITE(year_str,'(I4)') datein%Year
! make filename and get data
filein_rad= TRIM(base_dir)//'rad/'//year_str//mon_str//day_str//'_rad.flt'
!print *, filein_rad
ALLOCATE(data_out(np))
CALL read_binary_file(filein_rad,data_out)
all_data(:,1) = data_out
DEALLOCATE(data_out)

filein_rain= TRIM(base_dir)//'rain/'//year_str//mon_str//day_str//'_rain.flt'
ALLOCATE(data_out(np))
CALL read_binary_file(filein_rain,data_out)
all_data(:,2) = data_out
DEALLOCATE(data_out)

filein_tmax= TRIM(base_dir)//'tmax/'//year_str//mon_str//day_str//'_tmax.flt'
ALLOCATE(data_out(np))
CALL read_binary_file(filein_tmax,data_out)
all_data(:,3) = data_out
DEALLOCATE(data_out)

filein_tmin= TRIM(base_dir)//'tmin/'//year_str//mon_str//day_str//'_tmin.flt'
ALLOCATE(data_out(np))
CALL read_binary_file(filein_tmin,data_out)
all_data(:,4) = data_out
DEALLOCATE(data_out)

filein_vph09= TRIM(base_dir)//'vph09/bom-vph09-day-'//year_str//mon_str//day_str//'-'//year_str//mon_str//day_str//'.flt'
ALLOCATE(data_out(np))
CALL read_binary_file(filein_vph09,data_out)
all_data(:,5) = data_out
DEALLOCATE(data_out)

filein_vph15=TRIM(base_dir)//'vph15/bom-vph15-day-'//year_str//mon_str//day_str//'-'//year_str//mon_str//day_str//'.flt'
ALLOCATE(data_out(np))
CALL read_binary_file(filein_vph15,data_out)
all_data(:,6) = data_out
DEALLOCATE(data_out)
END SUBROUTINE get_all_bom_data

SUBROUTINE write_nc_file(dir_write,file_prefix,datein,data_in,lat1d,lon1d)
USE netcdf
USE DateFunctions
USE GlobalDefs
!INTEGER, INTENT(IN) :: year, month, day
type(dmydate), INTENT(IN) :: datein
CHARACTER(LEN=*), INTENT(IN) :: dir_write, file_prefix
REAL, DIMENSION(:,:,:), INTENT(IN) :: data_in
!REAL, DIMENSION(:,:), INTENT(IN) :: data_in
REAL, DIMENSION(:), INTENT(IN) :: lat1d, lon1d
integer :: ncid, lat_dimid, lon_dimid, time_dimid, lat_varid, lon_varid, time_varid, dimids(3)
integer :: Temp_varid, Speed_varid, Hum_varid, LW_varid, SW_varid, &
           Press_varid, Precip_varid, Coszen_varid
CHARACTER(LEN=200) :: FILE_NAME
CHARACTER(LEN=2) :: day_str, mon_str
CHARACTER(LEN=4) :: year_str
REAL,ALLOCATABLE, DIMENSION(:,:) :: Temp_out, Speed_out, Hum_out, LW_out, &
                                    SW_out, Press_out, Precip_out, Coszen_out 
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Temp2d, Speed2d, Hum2d, LW2d, &
                                    SW2d, Press2d, Precip2d, Coszen2d
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Temp2d_reverse, Speed2d_reverse, Hum2d_reverse, LW2d_reverse, &
                                    SW2d_reverse, Press2d_reverse, Precip2d_reverse, Coszen2d_reverse
INTEGER :: ii,i1,i2,i3
CHARACTER (LEN = *), PARAMETER :: UNITS = "units", LONGNAME="long_name",MISS="missing_value"
CHARACTER (LEN = *), PARAMETER :: Temp_UNITS = "K", Speed_UNITS = "m/s", &
                                  LW_UNITS = "W/m2", SW_UNITS = "W/m2", &
                                  Precip_UNITS = "mm/s", Press_UNITS = "hPa", &
                                  Coszen_UNITS = "-", Hum_UNITS = "kg/kg" , &
                                  Lat_UNITS = "degrees_north", Lon_UNITS ="degrees_east"

CHARACTER (LEN = *), PARAMETER :: Temp_LONGNAME = "Near surface air temperature", &
                                  Speed_LONGNAME = "Scalar windspeed", &
                                  LW_LONGNAME = "Surface incident longwave radiation", &
                                  SW_LONGNAME = "Surface incident shortwave radiation", &
                                  Precip_LONGNAME = "Rainfall Rate", &
                                  Press_LONGNAME = "Surface air pressure", &
                                  Coszen_LONGNAME = "Cosine of zenith angle", &
                                  Hum_LONGNAME = "Near surface specific humidity", &
                                  Lat_LONGNAME = "Latitude", Lon_LONGNAME = "Longitude"
REAL, PARAMETER :: miss_val = -999.

! pad zeros
IF (datein%Day < 10) THEN
      WRITE(day_str,'(A,I1.2)') '0', datein%Day
ELSE
      WRITE(day_str,'(I2)') datein%Day
ENDIF
IF (datein%Month < 10) THEN
      WRITE(mon_str,'(A,I1.2)') '0', datein%Month
ELSE
      WRITE(mon_str,'(I2)') datein%Month
ENDIF
WRITE(year_str,'(I4)') datein%Year
! make filename and get data
FILE_NAME = TRIM(dir_write//file_prefix//'_'//year_str//'-'//mon_str//'-'//day_str//'.nc')
!print *, FILE_NAME
call check(nf90_create(FILE_NAME, nf90_clobber, ncid) )
call check( nf90_def_dim(ncid, 'latitude', nrows, lat_dimid) )
call check( nf90_def_dim(ncid, 'longitude', ncols, lon_dimid) )
call check( nf90_def_dim(ncid, 'time'     , ntime,  time_dimid) )
dimids = (/lon_dimid, lat_dimid, time_dimid/)
call check( nf90_def_var(ncid, 'latitude', NF90_REAL, lat_dimid, lat_varid) )
call check( nf90_put_att(ncid, lat_varid, UNITS, Lat_UNITS) )
call check( nf90_put_att(ncid, lat_varid, LONGNAME, Lat_LONGNAME) )

call check( nf90_def_var(ncid, 'longitude', NF90_REAL, lon_dimid, lon_varid) )
call check( nf90_put_att(ncid, lon_varid, UNITS, Lon_UNITS) )
call check( nf90_put_att(ncid, lon_varid, LONGNAME, Lon_LONGNAME) )

call check( nf90_def_var(ncid, 'Tair', NF90_REAL, dimids, Temp_varid) )
call check( nf90_put_att(ncid, Temp_varid, UNITS, Temp_UNITS) )
call check( nf90_put_att(ncid, Temp_varid, LONGNAME, Temp_LONGNAME) )
call check( nf90_put_att(ncid, Temp_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'Wind', NF90_REAL, dimids, Speed_varid) )
call check( nf90_put_att(ncid, Speed_varid, UNITS, Speed_UNITS) )
call check( nf90_put_att(ncid, Speed_varid, LONGNAME, Speed_LONGNAME) )
call check( nf90_put_att(ncid, Speed_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'LWdown', NF90_REAL, dimids, LW_varid) )
call check( nf90_put_att(ncid, LW_varid, UNITS, LW_UNITS) )
call check( nf90_put_att(ncid, LW_varid, LONGNAME, LW_LONGNAME) )
call check( nf90_put_att(ncid, LW_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'SWdown', NF90_REAL, dimids, SW_varid) )
call check( nf90_put_att(ncid, SW_varid, UNITS, SW_UNITS) )
call check( nf90_put_att(ncid, SW_varid, LONGNAME, SW_LONGNAME) )
call check( nf90_put_att(ncid, SW_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'Rainf', NF90_REAL, dimids, Precip_varid) )
call check( nf90_put_att(ncid, Precip_varid, UNITS, Precip_UNITS) )
call check( nf90_put_att(ncid, Precip_varid, LONGNAME, Precip_LONGNAME) )
call check( nf90_put_att(ncid, Precip_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'PSurf', NF90_REAL, dimids, Press_varid) )
call check( nf90_put_att(ncid, Press_varid, UNITS, Press_UNITS) )
call check( nf90_put_att(ncid, Press_varid, LONGNAME, Press_LONGNAME) )
call check( nf90_put_att(ncid, Press_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'Coszen', NF90_REAL, dimids, Coszen_varid) )
call check( nf90_put_att(ncid, Coszen_varid, UNITS, Coszen_UNITS) )
call check( nf90_put_att(ncid, Coszen_varid, LONGNAME, Coszen_LONGNAME) )
call check( nf90_put_att(ncid, Coszen_varid, MISS, miss_val))

call check( nf90_def_var(ncid, 'Qair', NF90_REAL, dimids, Hum_varid) )
call check( nf90_put_att(ncid, Hum_varid, UNITS, Hum_UNITS) )
call check( nf90_put_att(ncid, Hum_varid, LONGNAME, Hum_LONGNAME) )
call check( nf90_put_att(ncid, Hum_varid, MISS, miss_val))

call check( nf90_enddef(ncid) )
call check( nf90_put_var(ncid, lat_varid, lat1d) )
call check( nf90_put_var(ncid, lon_varid, lon1d) )

!ShortWave_out = hMM(:,1,:)
!LongWave_out = hMM(:,2,:)
!Precip_out = hMM(:,3,:)
!Wind_out = hMM(:,4,:)
!Temp_out = hMM(:,5,:)
!Hum_out = hMM(:,6,:)
!Press_out = hMM(:,7,:)
!coszen_out = hMM(:,8,:)

ALLOCATE(SW_out(np,ntime))
SW_out =data_in(:,1,:)
ALLOCATE(SW2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(SW_out(:,ii),SW2d(ii,:,:))
END DO
ALLOCATE(SW2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(SW2d,SW2d_reverse)
DEALLOCATE(SW2d)
CALL check( nf90_put_var(ncid,SW_varid,SW2d_reverse))
DEALLOCATE(SW2d_reverse)

ALLOCATE(LW_out(np,ntime))
LW_out =data_in(:,2,:)
ALLOCATE(LW2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(LW_out(:,ii),LW2d(ii,:,:))
END DO
ALLOCATE(LW2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(LW2d,LW2d_reverse)
DEALLOCATE(LW2d)
CALL check( nf90_put_var(ncid,LW_varid,LW2d_reverse))
DEALLOCATE(LW2d_reverse)

ALLOCATE(Precip_out(np,ntime))
Precip_out =data_in(:,3,:)   ! convert to mm/s
WHERE (Precip_out .GT. 0.)
        Precip_out = Precip_out / (24.0*3600.0)
END WHERE     

ALLOCATE(Precip2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(Precip_out(:,ii),Precip2d(ii,:,:))
END DO
ALLOCATE(Precip2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(Precip2d,Precip2d_reverse)
DEALLOCATE(Precip2d)
CALL check( nf90_put_var(ncid,Precip_varid,Precip2d_reverse))
DEALLOCATE(Precip2d_reverse)

ALLOCATE(Speed_out(np,ntime))
Speed_out =data_in(:,4,:)
ALLOCATE(Speed2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(Speed_out(:,ii),Speed2d(ii,:,:))
END DO
ALLOCATE(Speed2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(Speed2d,Speed2d_reverse)
DEALLOCATE(Speed2d)
CALL check( nf90_put_var(ncid,Speed_varid,Speed2d_reverse))
DEALLOCATE(Speed2d_reverse)

ALLOCATE(Temp_out(np,ntime))
Temp_out =data_in(:,5,:) + 273.15 ! convert to Kelvin
ALLOCATE(Temp2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(Temp_out(:,ii),Temp2d(ii,:,:))
END DO
ALLOCATE(Temp2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(Temp2d,Temp2d_reverse)
DEALLOCATE(Temp2d)
CALL check( nf90_put_var(ncid,Temp_varid,Temp2d_reverse))
DEALLOCATE(Temp2d_reverse)

ALLOCATE(Hum_out(np,ntime))
Hum_out =data_in(:,6,:)
ALLOCATE(Hum2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(Hum_out(:,ii),Hum2d(ii,:,:))
END DO
ALLOCATE(Hum2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(Hum2d,Hum2d_reverse)
DEALLOCATE(Hum2d)
CALL check( nf90_put_var(ncid,Hum_varid,Hum2d_reverse))
DEALLOCATE(Hum2d_reverse)

ALLOCATE(Press_out(np,ntime))
Press_out =data_in(:,7,:)
ALLOCATE(Press2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(Press_out(:,ii),Press2d(ii,:,:))
END DO
ALLOCATE(Press2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(Press2d,Press2d_reverse)
DEALLOCATE(Press2d)
CALL check( nf90_put_var(ncid,Press_varid,Press2d_reverse))
DEALLOCATE(Press2d_reverse)

ALLOCATE(Coszen_out(np,ntime))
Coszen_out =data_in(:,8,:)
ALLOCATE(Coszen2d(ntime,nrows,ncols))
DO ii=1,ntime
    CALL make_2d(Coszen_out(:,ii),Coszen2d(ii,:,:))
END DO
ALLOCATE(Coszen2d_reverse(ncols,nrows,ntime))
CALL reverse_dims(Coszen2d,Coszen2d_reverse)
DEALLOCATE(Coszen2d)
CALL check( nf90_put_var(ncid,Coszen_varid,Coszen2d_reverse))
DEALLOCATE(Coszen2d_reverse)

call check( nf90_close(ncid) )

END SUBROUTINE write_nc_file

SUBROUTINE check(status)
USE netcdf
INTEGER, INTENT (IN) :: status
IF(status /= nf90_noerr) THEN 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
END IF
END SUBROUTINE check  

SUBROUTINE make_2d(vector_1d,array_2d)
USE GlobalDefs
REAL, DIMENSION(:), INTENT(IN) :: vector_1d
REAL, DIMENSION(:,:), INTENT(OUT) :: array_2d 
INTEGER :: ind1(nrows), ind2(nrows)
INTEGER :: kk
ind1(1) = ncols
ind2(1) = 1
DO kk = 1,nrows-1
  ind1(kk+1) = ind1(kk) + ncols
  ind2(kk+1) = ind2(kk) + ncols
END DO
DO kk=1,nrows
    array_2d(kk,:) = vector_1d(ind2(kk):ind1(kk))
END DO
END SUBROUTINE make_2d

SUBROUTINE make_1d(array_2d,vector_1d)
USE GlobalDefs
REAL, DIMENSION(:), INTENT(OUT) :: vector_1d
REAL, DIMENSION(:,:), INTENT(IN) :: array_2d
INTEGER :: ind1(nrows), ind2(nrows)
INTEGER :: kk
ind1(1) = ncols
ind2(1) = 1
DO kk = 1,nrows-1
  ind1(kk+1) = ind1(kk) + ncols
  ind2(kk+1) = ind2(kk) + ncols
END DO
DO kk=1,nrows
    vector_1d(ind2(kk):ind1(kk)) = array_2d(kk,:)
END DO
END SUBROUTINE make_1d

SUBROUTINE reverse_dims(a_in,a_out)
USE GlobalDefs
REAL, DIMENSION(:,:,:), INTENT(IN) :: a_in
REAL, DIMENSION(:,:,:), INTENT(OUT) :: a_out
INTEGER :: i1, i2, i3
DO i1 = 1,ncols
   DO i2 = 1,nrows
       DO i3 = 1,ntime
           a_out(i1,i2,i3) = a_in(i3,i2,i1)
       END DO
   END DO
END DO
END SUBROUTINE reverse_dims

END MODULE fileio_module
