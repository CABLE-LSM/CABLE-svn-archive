!module cable_ncdf_module
!contains
subroutine predef_grid(latitude, longitude, node_gl, rows, row_length, mp,       &
                     npseudo,mydata)
   use netcdf
   implicit none
 
   ! we are passing these as they are declared i8 and r8
   integer*8 :: node_gl, rows, row_length, mp, npseudo
   REAL*8, dimension(mp) :: latitude, longitude 
   real, dimension(mp, npseudo) :: mydata 

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
               y_in(i) = y
               if( i > mp) print *, "Count has reached mp"
            endif
         end do
      end do
   end do

   do z = 1, nPseudo
      do i = 1, mp
         mydata(i,z) =   data_in(x_in(i), y_in(i), z)   
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
!program sum                                          !a: name of program
!  implicit none
!  
!    ! temporary variables
!    integer, parameter :: file_id = 15,file_daily_id=20,file_shift_daily_id =25
!    character(len=100) :: junk_head,label, input_file, temp_shift_days
!    
!    integer  :: shift_days = 20
!    integer, parameter :: days_of_year=365
!    integer, parameter :: months_of_year=12
!    integer :: i_day
!
!    ! define arrary for reading the monthly lai data
!    real :: monthly_lai_array(2,12)  ! colomum major function pay attehtion
!    
!    ! variables lai data before and after interpolation and shift 
!    real :: monthly_lai(months_of_year)
!    real :: daily_lai(days_of_year)
!    
!    !initialize the daily lai array
!    DATA daily_lai /days_of_year*-9999/
!    
!    ! switch whether the using default lai input and shifting days or input as argument
!    IF (command_argument_count() >1) THEN
!      print *, "using the argument variables"
!      call get_command_argument(1,input_file)
!      call get_command_argument(2,temp_shift_days)
!      print *, "argument one ****   ", input_file
!      print *, "argument two *** ",    temp_shift_days
!      ! convert the string into iterger
!      open(file_id,file = input_file)
!      read(temp_shift_days,'(I2)') shift_days
!    ELSE ! using the defualt input file 
!      open(file_id, file='north_lai.txt')  ! it will the lai from the north part of earth
!    END IF
!
!   ! read the monhtly lai data
!   ! read the junk head and print
!     read(file_id,'(A)'), junk_head
!     print *,'The head of the file is :', junk_head
!   ! read the data into monthly lai arrary
!     read(file_id,*),monthly_lai_array
!
!   ! extract monthly lai for interpolation 
!     monthly_lai(:) = monthly_lai_array(2,:)
!     close(file_id)  ! do not forget the close the file
!  
!    
!   ! montly lai data is ready, it is time interplocation into daily lai
!   ! write(*,6),monthly_lai
!   print *, "interpolate the lai and write out"    
!   call interp_linear_lai(monthly_lai,daily_lai) ! call interpolation function 
!   ! write out the interpolated daily out
!   open(file_daily_id,file='daily_lai.txt')
!   do i_day=1, days_of_year
!       write(file_daily_id,*) daily_lai(i_day)
!   end do

!   print *, "shifting the daily lai defined shifted days"
   
!   call advanced_shift_daily_lai(shift_days,daily_lai,monthly_lai) 
   ! write out the shifted daily lai for comparison
!   open(file_shift_daily_id,file="advanced_shift_daily_lai.txt")
!    
!   do i_day=1, days_of_year
!       write(file_shift_daily_id,*) daily_lai(i_day)
!   end do 
!
!   print *, "successfully interpolate and shifted the lai data"
!   ! write(*,6),daily_lai
!  
!5   format(f16.4)
!6   format(f16.4)
!
!end  program sum                                     !h: end of program 
 
! subroutine for the interpolate the monthly lai data into daily lai
  
subroutine interp_linear_lai(monthly_lai,daily_lai)
    implicit none
    integer,parameter :: days_of_year=365
    integer,parameter :: months_of_year=12
    real,dimension(months_of_year) :: monthly_lai
    real,dimension(days_of_year) :: daily_lai
    integer :: month_date(months_of_year) = (/ 15,46,75,106,136,167,197,228,259,289,320,350 /) 
    
    integer :: i_day, i_month, first_date_id, second_date_id, delta_days
    real  :: first_date_lai, second_date_lai
    real  :: delta_lai
    
    ! interpolation start
    do i_month=1, months_of_year-1  ! do not go to last month 
       first_date_id   = month_date(i_month)
       second_date_id  = month_date(i_month+1)
       first_date_lai  = monthly_lai(i_month)
       second_date_lai = monthly_lai(i_month+1)
       delta_days      = second_date_id -first_date_id
       delta_lai       = (second_date_lai-first_date_lai)/delta_days
       
       do i_day=1, delta_days+1     ! start date loop
          daily_lai(first_date_id+i_day-1) = first_date_lai + (i_day-1)*delta_lai
       end do ! end with the date loop
       
    end do ! end with month loop 
    ! deal with special case day 1 to 15 and 350 to 365
    daily_lai(1:14)    =  monthly_lai(1) ! 
    daily_lai(351:365) =  monthly_lai(12)!


    ! end with interpolation
    return
end subroutine interp_linear_lai   

! define the shift subrountine
subroutine advanced_shift_daily_lai(shift_days,daily_lai,monthly_lai)
    integer :: shift_days
    integer,parameter :: days_of_year=365
    integer,parameter :: months_of_year=12
    real :: daily_lai(days_of_year), monthly_lai(months_of_year)

    integer :: pos_max(1),pos_min(1)
    real :: diff_monthly_lai(months_of_year-1)
    
    ! calculate the difference of the monthly LAI in order to provide condition for shift condition
    diff_monthly_lai = monthly_lai(2:months_of_year)-monthly_lai(1:(months_of_year-1))
    diff_monthly_lai = abs(diff_monthly_lai) 

    ! locate the position of the max and min lai vaule 
    pos_max = MAXLOC(daily_lai)
    pos_min = MINLOC(daily_lai) 

    IF(maxval(diff_monthly_lai) > 1.3) THEN 
      print *, "strange value, no shift"
    ELSE ! normal case shift the lai 
       IF(pos_max(1) > 90 .AND. pos_max(1) < 300)  THEN  ! north part   
         daily_lai(1:(pos_max(1)-shift_days))              = daily_lai((shift_days+1):pos_max(1))

         daily_lai( (pos_max(1)-shift_days+1):pos_max(1)) = daily_lai(pos_max(1)) ! repeat the max value shift times
       ELSE 
         daily_lai(pos_min(1):(days_of_year-shift_days))   = daily_lai((pos_min(1)+shift_days):days_of_year)  !(south part)

         daily_lai( (days_of_year-shift_days+1):days_of_year) = daily_lai(days_of_year)        
        END IF ! end with north south loop
    END IF ! end with normal if
     
end subroutine advanced_shift_daily_lai

