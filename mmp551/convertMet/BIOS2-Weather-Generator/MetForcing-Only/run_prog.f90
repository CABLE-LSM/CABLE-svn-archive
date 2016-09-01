PROGRAM run_prog
! jtk561: Main program to read input data, run weather generation, and write to Alma
! netcdf for any one particular day
! IMPORTANT NOTES:
! The IO part of this code is also hard coded to read the Bom data at /g/data1/fj2/AWAP/, and
! assumes a particular directory structure. If this code is to be used using
! inputs from other sources, changes will need to be made to GlobalDefs.f90 and
! fileio_module.f90, and this program. The rest of the code should not need any
! changes
USE TypeDef
USE GlobalDefs
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
USE SubDiurnalMetModule
USE GetSubDiurnalMetModule 
USE fileio_module
IMPLICIT NONE
INTEGER :: in_day, in_month, in_year, i, ii, jj
type(dmydate) :: Today, NextDay, PrevDay
REAL, DIMENSION(:,:), ALLOCATABLE :: MM,MMprev,MMnext ! contains all met variables for day,prev,next day
REAL, DIMENSION(:,:), ALLOCATABLE :: VV ! 2D fields needed by code
CHARACTER(LEN=200) :: tmp_str
REAL, DIMENSION(:), ALLOCATABLE :: lat1d, lon1d
REAL, DIMENSION(:,:), ALLOCATABLE :: lat2d, lat2d1
REAL, DIMENSION(:), ALLOCATABLE :: Latreshape
REAL, DIMENSION(:), ALLOCATABLE :: TTime
REAL, DIMENSION(:,:,:), ALLOCATABLE :: hMM
CHARACTER(LEN=100) :: dir_write="./", file_prefix

! Read input arguments and check if correct number given
i = iargc()
IF (i .LT. 3) THEN
     WRITE (*,*) "Wrong usage of program!"
     WRITE (*,*) "Correct Use:"
     WRITE (*,*) "./run_prog year month day"
     STOP
END IF
! read and covert to integer
CALL getarg(1, tmp_str)
READ (tmp_str, *) in_year
CALL getarg(2, tmp_str)
READ(tmp_str, *) in_month
CALL getarg(3, tmp_str)
READ(tmp_str, *) in_day

! Make Today date
Today%Day = in_day
Today%Month = in_month
Today%Year = in_year

! Calc Next day and previous dates
NextDay = AddDay(Today,1)
PrevDay = SubDay(Today,1)

! extract BoM data for Today, NextDay and PrevDay, as needed by the weather
! generator in the right order
ALLOCATE(MM(np,6))
ALLOCATE(MMprev(np,6))
ALLOCATE(MMnext(np,6))
CALL get_all_bom_data(Today,MM)
CALL get_all_bom_data(PrevDay,MMprev)
CALL get_all_bom_data(NextDay,MMnext)

! define the Ttime and VV arays as required by the weather generator
ALLOCATE(VV(np,72)) ! the orginal code checks for 72, need to change later
! create lat array
ALLOCATE(lat1d(nrows))
ALLOCATE(lon1d(ncols))
DO ii = 1,nrows
lat1d(ii) = yllcorner + (ii*cellsize)
END DO

lat1d = lat1d(nrows:1:-1)
!PRINT *, SHAPE(lat1d)
!PRINT *, lat1d

DO jj = 1,ncols
lon1d(jj) = xllcorner + (jj*cellsize)
END DO
ALLOCATE(lat2d(nrows,ncols))
!ALLOCATE(lon2d(nrows,ncols))
DO jj = 1,ncols
lat2d(:,jj) = lat1d
END DO
!DEALLOCATE(lat1d)
!DO ii =1,nrows
!lon2d(ii,:) = lon1d
!END DO
!DEALLOCATE(lon1d)

ALLOCATE(Latreshape(np))
!Latreshape = RESHAPE(lat2d,(/np/))
CALL make_1d(lat2d,Latreshape)
!DEALLOCATE(lat2d)

ALLOCATE(lat2d1(nrows,ncols))
CALL make_2d(Latreshape,lat2d1)

!PRINT *, lat2d - lat2d1

!PRINT *, SHAPE(lat2d)
!PRINT *, SHAPE(lat2d1)

!PRINT *, Latreshape(1:841)
!PRINT *, 'next'
!PRINT *, Latreshape(842:(842+840))

VV(:,3) = Latreshape
DEALLOCATE(Latreshape)

! define the Ttime array, only TTDay, TTMonth and TTYear seem to be used
ALLOCATE(TTime(5))
TTime(2) = Today%Day
TTime(3) = Today%Month
TTime(4) = Today%Year
!  Time1y     => TargetTTime(01)
!  TTDay      => TargetTTime(02)
!  TTMonth    => TargetTTime(03)
!  TTYear     => TargetTTime(04)
!  TTEndMth   => TargetTTime(05)
!real(sp),pointer:: Time1y       ! [y]           ! 01 run time at end of step
!real(sp),pointer:: TTDay        ! [d]           ! 02 current date: day
!real(sp),pointer:: TTMonth      ! [mth]         ! 03 current date: month
!real(sp),pointer:: TTYear       ! [y]           ! 04 current date: year
!real(sp),pointer:: TTEndMth     ! [-]           ! 05 (0,1) = last day of month(N,Y)

!!       hMM(np,nhMM,ntime) = hourly met variables
ALLOCATE(hMM(np,8,ntime))
! now CALL the weather generator
CALL GetSubdiurnalMet(Ttime,MM,MMprev,MMnext,VV,hMM) 
DEALLOCATE(TTime)
DEALLOCATE(MM)
DEALLOCATE(MMprev)
DEALLOCATE(MMnext)
DEALLOCATE(VV)

file_prefix = TRIM('BIOS2AWAP')
! dir_write is output dir, datein is used to make filename
CALL write_nc_file(TRIM(dir_write),TRIM(file_prefix),Today,hMM,lat1d,lon1d)
!  hFsd      => TargethMM(:,01,:)
!  hFld      => TargethMM(:,02,:)
!  hPrecip       => TargethMM(:,03,:)
!  hUa       => TargethMM(:,04,:)
!  hTc       => TargethMM(:,05,:)
!  hqv       => TargethMM(:,06,:)
!  hpmb      => TargethMM(:,07,:)
!  hcoszen       => TargethMM(:,08,:)

DEALLOCATE(hMM)
DEALLOCATE(lat1d)
DEALLOCATE(lon1d)

END PROGRAM run_prog
