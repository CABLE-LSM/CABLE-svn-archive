! Gab Abramowitz, University of New South Wales
! gabsun@gmail.com
! Writes a NetCDF binary file from text met forcing, following ALMA
! convention as well as being Grads readable. This version to read
! fluxnet format files from multiple input files. 

!**** NEE NOT ALMA UNITS - impractically small numbers

! **MAKE SURE** all parameters are set:
!   timeunits, timestpunits, timestepsize, latitude, longitude,
!   elevation 

! Outputs (LH,SH,Gfl etc) with missing values will result in no 
! variable in netcdf file******

! Use Makefile to build and run.

! Fluxnet flags refer to:  
!  "0" original value,
!  "1" missing in original,
!  "2" rejected from original,
!  "3" filled by redundant,
!  "4" last-second removes by data provider

PROGRAM fluxnet_to_netcdf 
  USE netcdf
  IMPLICIT NONE

!!!!!!!!==================USER ADJUST VARIABLES=====================!!!!!!!!!!!!

  ! Input/output details
  REAL, PARAMETER    :: latitude = 23.167  ! site latitude   ! BP
  REAL, PARAMETER    :: longitude = 112.533 ! site longitude   ! BP
  INTEGER, PARAMETER :: timestepsize=1800  ! timestep size in seconds
                                           ! of OUTPUT file
  REAL, PARAMETER :: elevation = 300.0 ! site elevation, to create pressure
                                       ! where missing   ! BP
  ! To use less than all the data in the first input file, specify the
  ! starting date in the following variable (only sensitive to hours):
  CHARACTER(LEN=20), PARAMETER :: start_time='2003-01-01 00:01:00' ! BP
  INTEGER, PARAMETER :: nfiles = 1   ! # met forcing/flux files ! BP
  INTEGER, PARAMETER :: za = 1       ! # vertical atmospheric layers
  INTEGER, PARAMETER :: zs = 1       ! # vertical soil layers
  CHARACTER(LEN=43), DIMENSION(nfiles) :: fileinmet ! set below
  CHARACTER(LEN=50),DIMENSION(nfiles) :: fileinflux ! set below
  CHARACTER(LEN=50),DIMENSION(nfiles) :: fileinlai  ! set below  ! BP
  DATA fileinmet(1) / 'dinghushanmet2003c.met' /    ! BP
  DATA fileinflux(1) / 'dinghushanmet2003c.flx'/     ! BP
  DATA fileinlai(1) / 'dhslai03.dat' /      ! BP
!  DATA fileinlai / nfiles*'noLAI.txt'/     ! BP
  CHARACTER(LEN=*), PARAMETER :: sitename='Dinghushan'
  CHARACTER(LEN=*), PARAMETER :: datalength='One years' ! BP
  CHARACTER(LEN=*), PARAMETER :: fileout = 'metDH.nc'   ! BP
  CHARACTER(LEN=*), PARAMETER :: contact='BP, bernard.pak@csiro.au'
  ! "compresstype": ='half' will only record timesteps in the fluxnet data that
  !  arerecorded twice, and use 0 for those recorded once in the PAR and SWdown
  !  fields.
  ! "compresstype": ='avgd' will halve the number of data by averaging 
  !  consecutive timesteps (e.g. half-hourly to hourly data).
  ! "compresstype": ='dont' will use every timestep in the dataset.
  CHARACTER(LEN=4), PARAMETER :: compresstype='dont'
  LOGICAL, PARAMETER :: fixedlai=.FALSE. ! Are we using fixed or varying LAI?
  LOGICAL, PARAMETER :: flag = .TRUE. ! record variable gapfilling flags
  ! If site is recognised, fix problems with data?
  LOGICAL, PARAMETER :: fixdata = .TRUE.
  TYPE par_type
     INTEGER :: isoil(1,1) = 8    ! BP
     INTEGER :: iveg(1,1) = 1    ! BP
     REAL :: frac4(1,1) = 0.0     ! BP
     REAL :: hc(1,1) = 18.0        ! BP
     REAL :: za(1,1) = 38.0       ! should be height of flux tower. BP
  END TYPE par_type
  TYPE(par_type) :: par
    ! Output units for all variables :
  TYPE units_type        
     CHARACTER(LEN=5) :: SWdown='W/m^2'
     CHARACTER(LEN=5) :: LWdown='W/m^2'
     CHARACTER(LEN=5) :: Qle='W/m^2'
     CHARACTER(LEN=5) :: Qh='W/m^2'
     CHARACTER(LEN=5) :: Qg='W/m^2'
     CHARACTER(LEN=5) :: SWnet='W/m^2'
     CHARACTER(LEN=5) :: LWnet='W/m^2'
     CHARACTER(LEN=3) :: CO2air='ppm'
     CHARACTER(LEN=4) :: Rainf='mm/s'
     CHARACTER(LEN=4) :: Snowf='mm/s'
     CHARACTER(LEN=8) :: Tair='K'
     CHARACTER(LEN=2) :: PSurf='Pa'
     CHARACTER(LEN=3) :: Wind='m/s'
     CHARACTER(LEN=5) :: Qair='kg/kg'
     CHARACTER(LEN=5) :: Rnet='W/m^2'    ! not alma variable
     CHARACTER(LEN=1) :: SoilTemp='K'
     CHARACTER(LEN=9) :: NEE='umol/m2/s' ! not alma units
     CHARACTER(LEN=1) :: LAI='-'         ! not alma variable
  END TYPE units_type
  TYPE(units_type) :: units
!!!!!!!!!!!!!----------------END USER ADJUSTED VARIABLES---------------!!!!!!!!!!!!!!!
  
  INTEGER,PARAMETER  :: mtin = nfiles*17520 ! # met forcing timesteps TOTAL in INPUT files
  INTEGER :: mtmid ! used if some of first file is to be ignored
  INTEGER :: mtout ! # met tsteps total in output file, set by 'compresstype'
  INTEGER :: tdrop ! 
  INTEGER :: xv = 1 ! grid x-coord #
  INTEGER :: yv = 1 ! grid y-coord #
  INTEGER,PARAMETER :: lineskip=3 ! Skip how many lines at top of text file?
                              ! 2 for Gab's met file, 3 for normal fluxnet file
                              ! but 2 for all flux files (use lineskip-1) BP
  CHARACTER(LEN=33) :: timeunits    ! units for time variable
  CHARACTER(LEN=35) :: timestpunits ! units for timestep variable
  CHARACTER(LEN=20) :: timeorigin   ! time origin 

  ! Define met type:
   TYPE met_input_type                      
     REAL, ALLOCATABLE,DIMENSION(:,:,:,:) :: CO2air   ! CO2 concentration
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: SWdown  !downward short-wave radiation 
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: LWdown  ! downward long-wave radiation 
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: Rainf  ! rainfall 
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: Snowf  ! snowfall 
     REAL, ALLOCATABLE,DIMENSION(:,:,:,:) :: Tair  ! surface air temperature 
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: PSurf  ! surface air pressure 
     REAL, ALLOCATABLE,DIMENSION(:,:,:,:) :: Wind   ! surface wind speed 
     REAL, ALLOCATABLE,DIMENSION(:,:,:,:) :: Qair   ! surface specific humidity 
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: LAI ! leaf area index timeseries
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:) :: CO2air_flag   ! 
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: SWdown_flag   ! equivalent
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: LWdown_flag   !  gapfilling
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: Rainf_flag   !     flags
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: Snowf_flag   !       |
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:) :: Tair_flag   !      |
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: PSurf_flag   !      \ /
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:) :: Wind_flag    !     |
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:) :: Qair_flag    ! 
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: LAI_flag  ! 
  END TYPE met_input_type  
  REAL :: templai

  ! Define output type:
  TYPE fluxstate_output_type
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: Qle ! latent heat flux
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: Qh  ! sensible heat flux
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: Qg  ! ground heat flux
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: Rnet ! net radiation
     REAL, ALLOCATABLE,DIMENSION(:,:,:) :: NEE  ! net ecosystem exchange
     REAL, ALLOCATABLE,DIMENSION(:,:,:,:) :: SoilTemp ! soil temp
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: Qle_flag ! gapfilling flag for Qle
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: Qh_flag  ! 
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: Qg_flag  ! 
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: Rnet_flag !
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: NEE_flag  ! 
     INTEGER, ALLOCATABLE,DIMENSION(:,:,:,:) :: SoilTemp_flag ! 
  END TYPE fluxstate_output_type
 
  ! 'finclude' determines which "fluxes" are going to be included in the netcdf file,
  LOGICAL, DIMENSION(6) :: finclude ! 1:Qle(LH),2:Qh(SH),3:Qg,4:Rnet,5:NEE,6:Stemp
  TYPE (met_input_type),TARGET       :: metin 
  TYPE (fluxstate_output_type),TARGET :: flstin
  TYPE (met_input_type),TARGET       :: met2
  TYPE (fluxstate_output_type),TARGET :: flst2
  TYPE (met_input_type), POINTER       :: met 
  TYPE (fluxstate_output_type), POINTER :: flst
  REAL(KIND(0.0D0)), ALLOCATABLE, DIMENSION(:),TARGET :: timein
  REAL(KIND(0.0D0)), ALLOCATABLE, DIMENSION(:),TARGET :: time2
  REAL(KIND(0.0D0)),POINTER,DIMENSION(:) :: time
  REAL,ALLOCATABLE,DIMENSION(:),TARGET :: timestpin
  REAL,ALLOCATABLE,DIMENSION(:),TARGET :: timestp2
  REAL,POINTER,DIMENSION(:) :: timestp
  REAL :: skipr ! skip column in read
  INTEGER :: skipi ! skip column in read
  INTEGER :: tscounter ! timestep counter for output file
  REAL :: sph ! temporary specific humidity

  ! I/o declarations
  LOGICAL			:: eof	 ! end of file?  
  INTEGER			:: i,j,k,l,m,n,o,bpcount ! do counters  BP
  INTEGER			:: ios	 ! I/O status
  LOGICAL			:: is_open
  INTEGER, PARAMETER:: unitin = 12, unitout = 11 ! input unitin (open & close here) 

  ! Declarations for netcdf writing 
  INTEGER :: start_out1,count_out1,start_out3(3),count_out3(3),start_out4(4), &
       count_out4a(4),start_out2(2),count_out2(2),count_out4s(4)
  INTEGER :: xID,yID,zID,zsID,tID  ! Dimension IDs
  INTEGER :: SWdownID, LWdownID, RainfID, SnowfID, TairID, PSurfID, WindID, QairID, &
       CO2airID, nav_lonID, nav_latID, xvID, yvID, levelID, timeID, timestpID, &
       QleID, QhID, QgID, RnetID, NEEID, SoilTempID, LAIID, isoilID, ivegID, &
       frac4ID,hcID, zaID! Variable IDs
  INTEGER :: SWdownfID, LWdownfID, RainffID, SnowffID, TairfID, &
       PSurffID, WindfID, QairfID, CO2airfID, &
       QlefID, QhfID, QgfID, RnetfID, NEEfID, SoilTempfID, LAIfID ! Variable flagIDs
  INTEGER :: ncid, status, out1(1)=1
  INTEGER :: ctr ! counter for file reading
  CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
  CHARACTER(LEN=3)  :: smonth ! starting month of dataset
  INTEGER :: shod,sdoy,syear,smonthi,sdoytmp ! starting time variables
  INTEGER :: tstoendyr ! timesteps to end of first dataset year
  INTEGER :: tsyrs     ! timesteps to data fix points from whole years
  INTEGER :: tstofix   ! timesteps to data fix period
  INTEGER :: startfix, endfix ! timestep numbers of suspect data range

  INTEGER :: day1, jday   ! BP
  REAL :: twoyrLAI(730), lai1, dayprecip   ! BP

  ! Write site name:
  WRITE(*,*) 'Site: ', sitename, '     Dataset length: ', datalength

  ! Create all starting time variables:
  timeunits='seconds since '//start_time
  timestpunits='timesteps since '//start_time
  ! Determine the number of timesteps to be dropped from first file
  ! (use internal file to convert from character to integer):
  READ(timeunits(15:18),*) syear    ! integer year
  READ(timeunits(20:21),*) smonthi  ! integer month
  READ(timeunits(23:24),*) sdoytmp  ! integer day of that month
  READ(timeunits(26:27),*) shod     ! starting hour of day 
  SELECT CASE(smonthi)
  CASE(1) 
     smonth='JAN'
     sdoy=sdoytmp
  CASE(2) 
     smonth='FEB'
     sdoy=sdoytmp+31
  CASE(3) 
     smonth='MAR'
     sdoy=sdoytmp+59
  CASE(4)
     smonth='APR'
     sdoy=sdoytmp+90
  CASE(5)
     smonth='MAY'
     sdoy=sdoytmp+120
  CASE(6)
     smonth='JUN'
     sdoy=sdoytmp+151
  CASE(7)
     smonth='JUL'
     sdoy=sdoytmp+181
  CASE(8)
     smonth='AUG'
     sdoy=sdoytmp+212
  CASE(9)
     smonth='SEP'
     sdoy=sdoytmp+243
  CASE(10)
     smonth='OCT'
     sdoy=sdoytmp+273
  CASE(11)
     smonth='NOV'
     sdoy=sdoytmp+304
  CASE(12) 
     smonth='DEC'
     sdoy=sdoytmp+334
  CASE DEFAULT
     WRITE(*,*)'Could not interpret month (',timeunits(20:21),&
          ') in start_time variable.'
     STOP
  END SELECT
  timeorigin=start_time(1:5)//smonth//start_time(8:20)

  ! Check that starting years match fluxnet file read in:
  IF(timeunits(17:18)/=fileinmet(1)(16:17)) THEN     ! filenames are shorter. BP
     WRITE(*,*) 'Year of first met file (',fileinmet(1)(16:17), &
          ') does not match starting year in start_time variable (', &
          timeunits(17:18), ').'
     STOP
  END IF
!  IF(timeunits(17:18)/=fileinmet(1)(28:29)) THEN
!     WRITE(*,*) 'Year of first met file (',fileinmet(1)(28:29), &
!          ') does not match starting year in start_time variable (', &
!          timeunits(17:18), ').'
!     STOP
!  END IF
  ! Set number of timesteps to be dropped from first dataset
  ! (timestep # to start reading - 1)
  tdrop=INT(2.0*(24*(sdoy-1)+shod))
  WRITE(*,*) 'Dropping first',tdrop,'timesteps from first met file ', &
       '(begins: ',timeorigin,').'  
  ! Check timestepsize and compresstype compatibility:
  IF(timestepsize/=1800.0.AND.compresstype=='dont') THEN
     WRITE(*,*) 'No data compression but timestep size not native fluxnet!'
     STOP
  ELSE IF(timestepsize==1800.0.AND.compresstype/='dont') THEN
     WRITE(*,*) 'Data being compressed but timestep size unchanged!'
     STOP
  END IF

  ! The number of timesteps after the first portion (no observations) of the
  ! first data file is removed:
  mtmid = mtin-tdrop

  ALLOCATE(timein(mtmid), timestpin(mtmid))

  ! Allocate initial met/flux variable sizes:
  ALLOCATE(metin%SWdown(1,1,mtmid),metin%LWdown(1,1,mtmid),metin%Rainf(1,1,mtmid),& 
       metin%Snowf(1,1,mtmid), metin%PSurf(1,1,mtmid),metin%CO2air(1,1,za,mtmid), &
       metin%Tair(1,1,za,mtmid),metin%Wind(1,1,za,mtmid), metin%Qair(1,1,za,mtmid), &
       flstin%Qle(1,1,mtmid),flstin%Qh(1,1,mtmid),flstin%Qg(1,1,mtmid), &
       flstin%NEE(1,1,mtmid),flstin%Rnet(1,1,mtmid),flstin%SoilTemp(1,1,zs,mtmid), &
       metin%LAI(1,1,mtmid))
  ! Allocate initial met/flux variable flag sizes:
  ALLOCATE(metin%SWdown_flag(1,1,mtmid),metin%LWdown_flag(1,1,mtmid),&
       metin%Rainf_flag(1,1,mtmid),metin%Snowf_flag(1,1,mtmid), &
       metin%PSurf_flag(1,1,mtmid),metin%CO2air_flag(1,1,za,mtmid), &
       metin%Tair_flag(1,1,za,mtmid),metin%Wind_flag(1,1,za,mtmid), &
       metin%Qair_flag(1,1,za,mtmid), flstin%Qle_flag(1,1,mtmid), &
       flstin%Qh_flag(1,1,mtmid),flstin%Qg_flag(1,1,mtmid), &
       flstin%NEE_flag(1,1,mtmid),flstin%Rnet_flag(1,1,mtmid), &
       flstin%SoilTemp_flag(1,1,zs,mtmid), metin%LAI_flag(1,1,mtmid))

  !------------------------- BEGIN READING --------------------------------------------

  tscounter=1 ! initialise
  
  DO k=1,nfiles ! for each met forcing/flux file
     ! Open met file:
     INQUIRE(UNIT = unitin, OPENED = is_open)
     IF (.NOT. is_open) THEN
        OPEN (UNIT=unitin, FILE=fileinmet(k), STATUS='OLD', ACTION='READ', IOSTAT = ios)
        IF (ios /= 0) THEN
           WRITE(*,*) 'Error opening file ', fileinmet(k)
           STOP
        END IF
     END IF
     ! Read from file and write to 'met' one timestep at a time:
     ! starting with met file:
     DO j=1, lineskip
        READ (unitin, *, IOSTAT = ios)
     END DO
     IF (ios /= 0) WRITE(*,*) 'Error reading file ', fileinmet(k)

     DO i = tscounter, tscounter + mtin/nfiles-1
        ! If we have not yet reached the first timestep to read
        IF(i<=tdrop) THEN
           ctr=1 ! keep rewriting over the first entry in arrays
        ELSE
           ctr=i-tdrop
        END IF
        
        READ (unitin, *, IOSTAT = ios) & 
             skipi, skipr, metin%SWdown(1,1,ctr), metin%SWdown_flag(1,1,ctr), &
             skipr, skipi, metin%Tair(1,1,1,ctr), metin%Tair_flag(1,1,1,ctr), &
             flstin%SoilTemp(1,1,1,ctr), flstin%SoilTemp_flag(1,1,1,ctr), &
             metin%Qair(1,1,1,ctr), metin%Qair_flag(1,1,1,ctr), skipr, skipi, &
             metin%CO2air(1,1,1,ctr), metin%CO2air_flag(1,1,1,ctr), &
             flstin%Rnet(1,1,ctr),flstin%Rnet_flag(1,1,ctr), metin%Rainf(1,1,ctr), &
             metin%Rainf_flag(1,1,ctr), skipr, skipi, metin%Wind(1,1,1,ctr), &
             metin%Wind_flag(1,1,1,ctr), metin%PSurf(1,1,ctr), metin%PSurf_flag(1,1,ctr)
        timestpin(ctr)=ctr
        timein(ctr)=timestpin(ctr)*timestepsize
        eof = i==mtin.AND.ios==0
        IF (eof) THEN
           CLOSE(unitin)
        ELSE
           IF (ios /= 0) WRITE(*,*) 'Error reading met file', fileinmet(k)
        END IF
!      IF(metin%Rainf(1,1,ctr)==-9999.0) print*, i
     END DO
     CLOSE(unitin)

     ! Then flux file:
     INQUIRE(UNIT = unitin, OPENED = is_open)
     IF (.NOT. is_open) THEN
        OPEN (UNIT=unitin, FILE=fileinflux(k), STATUS='OLD', ACTION='READ', IOSTAT = ios)
        IF (ios /= 0) WRITE(*,*) 'Error opening file ', fileinflux(k)
     END IF
     DO j=1, lineskip -1           ! BP
        READ (unitin, *, IOSTAT = ios)
     END DO
     IF (ios /= 0) WRITE(*,*) 'Error reading file ', fileinflux(k)

     DO i =  tscounter, tscounter + mtin/nfiles-1       
        ! If we have not yet reached the first timestep to read
        IF(i<=tdrop) THEN
           ctr=1 ! keep rewriting over the first entry in arrays
        ELSE
           ctr=i-tdrop
        END IF
            
        READ (unitin, 12, IOSTAT = ios) & 
             skipi, skipr, flstin%NEE(1,1,ctr),flstin%NEE_flag(1,1,ctr), &
             flstin%Qle(1,1,ctr),flstin%Qle_flag(1,1,ctr),flstin%Qh(1,1,ctr), &
             flstin%Qh_flag(1,1,ctr),flstin%Qg(1,1,ctr),flstin%Qg_flag(1,1,ctr)
12      FORMAT(I4,F6.2,F10.3,I4,F9.2,I4,F9.2,I4,F9.2,I4)
        eof = i==mtin.AND.ios==0
        IF (eof) THEN
           CLOSE(unitin)
        ELSE
           IF (ios /= 0) THEN
              WRITE(*,*) 'Error reading flux file', fileinflux(k)
              STOP
           END IF
        END IF
     END DO
     CLOSE(unitin)

     ! Then lai file:
     IF(.NOT.fixedlai) THEN
        INQUIRE(UNIT = unitin, OPENED = is_open)
        IF (.NOT. is_open) THEN
           OPEN (UNIT=unitin, FILE=fileinlai(k), STATUS='OLD', &
                    & ACTION='READ', IOSTAT = ios)
           IF (ios /= 0) WRITE(*,*) 'Error opening file ', fileinlai(k)
        END IF
        WRITE(*,*) 'Start reading file ', fileinlai(k)
        READ (unitin, *, IOSTAT = ios) ! skip header line  ! BP
        DO i=1,365
           READ (unitin, *, IOSTAT = ios) skipi, templai   ! BP
!           READ (unitin, *, IOSTAT = ios) skipr, templai   ! bondville 97 only
           DO l = tscounter+(i-1)*86400/timestepsize, &
                & tscounter-1+i*86400/timestepsize
              metin%LAI(1,1,l)=templai
           END DO
           eof = i==mtin.AND.ios==0
           IF (eof) THEN
              CLOSE(unitin)
           ELSE
              IF (ios /= 0) WRITE(*,*) 'Error reading file 2', fileinlai(k)
           END IF
        END DO
     END IF
     CLOSE(unitin)
     ! increment counter
     tscounter=tscounter+mtin/nfiles
!     print *, 'l = ', l, ' tscounter = ', tscounter
!     IF(l /= tscounter-1) STOP 'tscounter and l are different in LAI section.'
     IF(l /= tscounter) STOP 'tscounter and l are different in LAI section.'  ! BP
  END DO

!---------------------------END READING FROM FILE-------------------------------
  
  
! -----------------------------SITE DEPENDENT FIXES----------------------------
  ! Site specific adjustments:
  IF(fixdata) THEN ! if we want to fix known problems with the data:
     SELECT CASE(sitename)
     CASE('Loobos','loobos') 
        ! Make sure period to be changed is in selected period for dataset:
        IF(syear<1998.OR.(syear==1998.AND.sdoy<=122)) THEN 
           ! ie changes are there to be made
           WRITE(*,*) '  Changes made to Loobos data:'
           ! Find relative timestep:
           IF(syear==1996) THEN
              ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
              tstoendyr = (365-sdoy)*24*2+(24-shod)*2
              ! Find # timesteps from whole years b/w start and fix:
              tsyrs = MAX(1997-syear,0)*365*48
              ! Find # timesteps to fix in fix year (1998)
              tstofix = 23*48
              ! Set index for beginning of fix period
              startfix=tstoendyr+tsyrs+tstofix
              ! Set index for end of fix period
              endfix=startfix+33*48
              ! First fix (replace with one year earlier):
              metin%SWdown(:,:,startfix:endfix) = &
                   metin%SWdown(:,:,startfix-48*365:endfix-48*365)
              metin%SWdown_flag(:,:,startfix:endfix) = 3
              WRITE(*,*) '    SWdown replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'        
              ! Find # timesteps to fix in fix year (1998)
              tstofix = 17*48
              ! Set index for beginning of fix period
              startfix=tstoendyr+tsyrs+tstofix
              ! Set index for end of fix period
              endfix=startfix+39*48
              ! Second fix (replace with one year earlier):
              metin%Tair(:,:,:,startfix:endfix) = &
                   metin%Tair(:,:,:,startfix-48*365:endfix-48*365)
              metin%Tair_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%PSurf(:,:,startfix:endfix) = &
                   metin%PSurf(:,:,startfix-48*365:endfix-48*365)
              metin%PSurf_flag(:,:,startfix:endfix) = 3
              WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%Qair(:,:,:,startfix:endfix) = &
                   metin%Qair(:,:,:,startfix-48*365:endfix-48*365)
              metin%Qair_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Qair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%CO2air(:,:,:,startfix:endfix) = &
                   metin%CO2air(:,:,:,startfix-48*365:endfix-48*365)
              metin%CO2air_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    CO2air replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%Wind(:,:,:,startfix:endfix) = &
                   metin%Wind(:,:,:,startfix-48*365:endfix-48*365)
              metin%Wind_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Wind replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              ! Find # timesteps to fix in fix year (1998)
              tstofix = 88*48
              ! Set index for beginning of fix period
              startfix=tstoendyr+tsyrs+tstofix
              ! Set index for end of fix period
              endfix=startfix+34*48
              ! Second fix (replace with one year earlier:
              metin%Wind(:,:,:,startfix:endfix) = &
                   metin%Wind(:,:,:,startfix-48*365:endfix-48*365)
              metin%Wind_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Wind replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           ELSE 
              WRITE(*,*) '  No replacement of data possible!'           
           END IF
        END IF
     CASE('Little Washita','Washita','washita','Little washita','little washita')
        WRITE(*,*) '  Changes made to Little Washita data:'
        ! Fix odd problem with relative humididty:
        WHERE(metin%Qair==-0.6)
           metin%Qair=60.0
        END WHERE
        WHERE(metin%Qair==-0.7) 
           metin%Qair=70.0
        END WHERE     
        WRITE(*,*) '    Relative humidity < 0 (fraction) converted to %.'
        !----- Make adjustment for temperature problems on 3/4/5 Oct 98:-----
        IF(syear==1996) THEN
           ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
           tstoendyr = (365-sdoy)*24*2+(24-shod)*2
           ! Find # timesteps from whole years b/w start and fix:
           tsyrs = 365*48
           ! Find # timesteps to fix in fix year (1998)
           tstofix = 13247
           ! Set index for beginning of fix period
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+26
           ! Apply changes (replace with specifically chosen period):
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,tstoendyr+12999:tstoendyr+12999+26)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           ! Set index for beginning of fix period
           tstofix = 13293
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+26
           ! Apply changes (replace with specifically chosen period):
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,tstoendyr+3834:tstoendyr+3834+26)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           ! Set index for beginning of fix period
           tstofix = 13342
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+26
           ! Apply changes (replace with specifically chosen period):
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,tstoendyr+31479:tstoendyr+31479+26)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           !------------------Now adjust for problems 15/12/96 - 2/1/97-------
           IF(sdoy<349) THEN
              ! Find # timesteps to fix in fix year (1996)
              tstofix = tstoendyr-17*48
              ! Set index for beginning of fix period
              startfix=tstofix
              ! Set index for end of fix period
              endfix=startfix+19*48
              ! Fix (replace with one year later):      
              metin%SWdown(:,:,startfix:endfix) = &
                   metin%SWdown(:,:,startfix+48*365:endfix+48*365)
              metin%SWdown_flag(:,:,startfix:endfix) = 3
              WRITE(*,*) '    SWdown replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%Tair(:,:,:,startfix:endfix) = &
                   metin%Tair(:,:,:,startfix+48*365:endfix+48*365)
              metin%Tair_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%Wind(:,:,:,startfix:endfix) = &
                   metin%Wind(:,:,:,startfix+48*365:endfix+48*365)
              metin%Wind_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Wind replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              ! Now all PSurf for period leading up to 1997:
              ! Set index for beginning of fix period
              startfix=1
              ! Set index for end of fix period (same as above)
              endfix=tstofix+19*48
              ! Fix (replace with one year later):      
              metin%PSurf(:,:,startfix:endfix) = &
                   metin%PSurf(:,:,startfix+48*365:endfix+48*365)
              metin%PSurf_flag(:,:,startfix:endfix) = 3
              WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           END IF
           !------------- Now pressure spike on 19/2/98------------------
           ! Find # timesteps to fix in fix year (1998)
           tstofix = 49*48+26
           ! Set index for beginning of fix period
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+2
           ! Apply changes (replace with specifically chosen period):
           metin%PSurf(:,:,startfix:endfix) = &
                metin%PSurf(:,:,startfix+2:endfix+2)
           metin%PSurf_flag(:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
        ELSE IF(syear==1997) THEN
           IF(sdoy/=1) THEN
              WRITE(*,*) 'Cannot make corrections for bad data!'
              STOP
           END IF
           ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
           tstoendyr = (365-sdoy)*24*2+(24-shod)*2
           ! Find # timesteps from whole years b/w start and fix:
           ! Find # timesteps to fix in fix year (1998)
           tstofix = 13247
           ! Set index for beginning of fix period
           startfix=tstoendyr+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+26
           ! Apply changes (replace with specifically chosen periods):
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,12999:12999+26)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           ! Set index for beginning of fix period
           tstofix = 13293
           startfix=tstoendyr+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+26
           ! Apply changes (replace with specifically chosen periods):
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,3834:3834+26)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           ! Set index for beginning of fix period
           tstofix = 13342
           startfix=tstoendyr+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+26
           ! Apply changes (replace with specifically chosen periods):
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,31479:31479+26)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           !------------- Now pressure spike on 19/2/98------------------
           ! Find # timesteps to fix in fix year (1998)
           tstofix = 49*48+26
           ! Set index for beginning of fix period
           startfix=tstoendyr+tstofix
           ! Set index for end of fix period (3 13hour drops in temp)
           endfix=startfix+2
           ! Apply changes (replace with specifically chosen period):
           metin%PSurf(:,:,startfix:endfix) = &
                metin%PSurf(:,:,startfix+2:endfix+2)
           metin%PSurf_flag(:,:,startfix:endfix) = 3 ! ie "filled by redundant"
           WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
        END IF
        !---------------------------------------------------------------------
     CASE('Norunda','norunda')
         IF(syear==1996) THEN
              ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
              tstoendyr = (365-sdoy)*24*2+(24-shod)*2
              ! Find # timesteps from whole years b/w start and fix:
              tsyrs = MAX(1997-syear,0)*365*48
              ! Find # timesteps to fix in fix year (1998)
              tstofix = 0
              ! Set index for beginning of fix period
              startfix=tstoendyr+tsyrs+tstofix
              ! Set index for end of fix period
              endfix=startfix+181*48
              ! First fix (replace with one year earlier):
              metin%SWdown(:,:,startfix:endfix) = &
                   metin%SWdown(:,:,startfix-48*365:endfix-48*365)
              metin%SWdown_flag(:,:,startfix:endfix) = 3
              WRITE(*,*) '    SWdown replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'        
              metin%Tair(:,:,:,startfix:endfix) = &
                   metin%Tair(:,:,:,startfix-48*365:endfix-48*365)
              metin%Tair_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Tair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%PSurf(:,:,startfix:endfix) = &
                   metin%PSurf(:,:,startfix-48*365:endfix-48*365)
              metin%PSurf_flag(:,:,startfix:endfix) = 3
              WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%Qair(:,:,:,startfix:endfix) = &
                   metin%Qair(:,:,:,startfix-48*365:endfix-48*365)
              metin%Qair_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Qair replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%CO2air(:,:,:,startfix:endfix) = &
                   metin%CO2air(:,:,:,startfix-48*365:endfix-48*365)
              metin%CO2air_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    CO2air replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
              metin%Wind(:,:,:,startfix:endfix) = &
                   metin%Wind(:,:,:,startfix-48*365:endfix-48*365)
              metin%Wind_flag(:,:,:,startfix:endfix) = 3
              WRITE(*,*) '    Wind replaced from timestep',startfix, &
                   'to timestep',endfix, ' - change flagged.'
           ELSE 
              WRITE(*,*) '  No replacement of data possible!'           
           END IF
     CASE('Harvard Forest','Harvard','harvard','harvard forest', &
          'Harvard forest')
        ! No fixes applied
        !---------------------------------------------------------------------
     CASE('Flak','Flakaliden','flak','flakaliden')
     CASE('Shidler','shidler')
     CASE('Castelporziano','castelporziano')
        IF(start_time=='1997-01-01 00:01:00') THEN
           WHERE(metin%Rainf_flag/=0)
              flstin%NEE_flag=5
              flstin%Qle_flag=5
              flstin%Qh_flag=5
           END WHERE
           metin%PSurf(1,1,1:17520)=metin%PSurf(1,1,17521:35040)
           metin%PSurf_flag(1,1,1:17520)=3
           WRITE(*,*) &
                '    PSurf replaced from timestep 1 to timestep 17520 - change flagged.'
        END IF
     CASE('Collelongo','collelongo')
     CASE('Sarrebourg/Hesse')
        IF(syear==1996) THEN
           ! Plug large hole in data 1/1/00 - 16/4/00:
           ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
           tstoendyr = (365-sdoy)*24*2+(24-shod)*2
           ! Find # timesteps from whole years b/w start and fix:
           tsyrs = 3*365*48
           ! Find # timesteps to fix in fix year (2000)
           tstofix = 0
           ! Set index for beginning of fix period
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period
           endfix=startfix+107*48
           ! First fix (replace with one year earlier):
           metin%SWdown(:,:,startfix:endfix) = &
                metin%SWdown(:,:,startfix-48*365:endfix-48*365)
           metin%SWdown_flag(:,:,startfix:endfix) = 3
           WRITE(*,*) '    SWdown replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.' 
           metin%Rainf(:,:,startfix:endfix) = &
                metin%Rainf(:,:,startfix-48*365:endfix-48*365)
           metin%Rainf_flag(:,:,startfix:endfix) = 3
           WRITE(*,*) '    Rainf replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.' 
           metin%Tair(:,:,:,startfix:endfix) = &
                metin%Tair(:,:,:,startfix-48*365:endfix-48*365)
           metin%Tair_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    Tair replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           metin%PSurf(:,:,startfix:endfix) = &
                metin%PSurf(:,:,startfix-48*365:endfix-48*365)
           metin%PSurf_flag(:,:,startfix:endfix) = 3
           WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           metin%Qair(:,:,:,startfix:endfix) = &
                metin%Qair(:,:,:,startfix-48*365:endfix-48*365)
           metin%Qair_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    Qair replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           metin%CO2air(:,:,:,startfix:endfix) = &
                metin%CO2air(:,:,:,startfix-48*365:endfix-48*365)
           metin%CO2air_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    CO2air replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           metin%Wind(:,:,:,startfix:endfix) = &
                metin%Wind(:,:,:,startfix-48*365:endfix-48*365)
           metin%Wind_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    Wind replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           ! NOW fix pressure problem for 1997:
           tsyrs = 0
           ! Find # timesteps to fix in fix year (2000)
           tstofix = -48
           ! Set index for beginning of fix period
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period
           endfix=startfix+366*48
           ! First fix (replace with two years later):
           metin%PSurf(:,:,startfix:endfix) = &
                metin%PSurf(:,:,startfix+48*365*2:endfix+48*365*2)
           metin%PSurf_flag(:,:,startfix:endfix) = 3
           WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
        END IF
     CASE('Bordeaux','bordeaux')
        WHERE(metin%Rainf==-9999.0)
           metin%Rainf=0.5
           metin%Rainf_flag=1
        END WHERE
        WHERE(metin%SWdown==-9999.0)
           metin%SWdown=150.0
           metin%SWdown_flag=1
        END WHERE
        WHERE(metin%Tair==-9999.0)
           metin%Tair=15.0
           metin%Tair_flag=1
        END WHERE
        WHERE(metin%Qair==-9999.0)
           metin%Qair=50.0
           metin%Qair_flag=1
        END WHERE
        WHERE(metin%Wind==-9999.0)
           metin%Wind=2.0
           metin%Wind_flag=1
        END WHERE
     CASE('Soroe','soroe')
        WHERE(metin%Qair<0.0)
           metin%Qair=50.0
           metin%Qair_flag=1
        END WHERE
        IF(syear==1996) THEN
           ! Plug large hole in SWdown data 13/2/99 - 21/2/99:
           ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
           tstoendyr = (365-sdoy)*24*2+(24-shod)*2
           ! Find # timesteps from whole years b/w start and fix:
           tsyrs = 2*365*48
           ! Find # timesteps to fix in fix year (2000)
           tstofix = 43*48
           ! Set index for beginning of fix period
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period
           endfix=startfix+9*48
           ! First fix (replace with one year earlier):
           metin%SWdown(:,:,startfix:endfix) = &
                metin%SWdown(:,:,startfix-48*365:endfix-48*365)
           metin%SWdown_flag(:,:,startfix:endfix) = 3
           WRITE(*,*) '    SWdown replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.' 
        END IF

     CASE('Weiden Brunnen')
        IF(syear==1996) THEN
           ! Plug large hole in CO2 data 1/1/97 - 20/5/97:
           ! Find # timesteps to end of 1st yr (timestep size is 1800, fluxnet):
           tstoendyr = (365-sdoy)*24*2+(24-shod)*2
           ! Find # timesteps from whole years b/w start and fix:
           tsyrs = 0
           ! Find # timesteps to fix in fix year (2000)
           tstofix = 0
           ! Set index for beginning of fix period
           startfix=tstoendyr+tsyrs+tstofix
           ! Set index for end of fix period
           endfix=startfix+140*48
           ! First fix (replace with one year later):
           metin%CO2air(:,:,:,startfix:endfix) = &
                metin%CO2air(:,:,:,startfix+48*365:endfix+48*365)
           metin%CO2air_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    CO2air replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.' 
        END IF
     CASE('Tharandt','tharandt')
        WHERE(metin%Rainf_flag/=0)
           flstin%NEE_flag=5
           flstin%Qle_flag=5
           flstin%Qh_flag=5
        END WHERE
        IF(start_time=='1996-01-01 00:01:00') THEN
           ! Set index for beginning of fix period
           startfix=1
           ! FIX problem with PSurf 1/1/96 - 24/2/97:
           ! Set index for end of fix period
           endfix=startfix+365*48+56*48
           ! First fix (replace with two years later):
           metin%PSurf(:,:,startfix:endfix) = &
                metin%PSurf(:,:,startfix+48*365*2:endfix+48*365*2)
           metin%PSurf_flag(:,:,startfix:endfix) = 3
           WRITE(*,*) '    PSurf replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           ! FIX problem with Wind 1/1/96 - 28/3/96:
           ! Set index for end of fix period
           endfix=startfix+87*48
           ! repalce with one year later
           metin%Wind(:,:,:,startfix:endfix) = &
                metin%Wind(:,:,:,startfix+48*365:endfix+48*365)
           metin%Wind_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    Wind replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.'
           ! FIX problem with CO2air 1/1/96 - 10/10/96:
           ! Set index for end of fix period
           endfix=startfix+283*48
           ! replace with one year later:
           metin%CO2air(:,:,:,startfix:endfix) = &
                metin%CO2air(:,:,:,startfix+48*365:endfix+48*365)
           metin%CO2air_flag(:,:,:,startfix:endfix) = 3
           WRITE(*,*) '    CO2air replaced from timestep',startfix, &
                'to timestep',endfix, ' - change flagged.' 
        END IF
        CASE('Brasschaat','brasschaat')
           WHERE(metin%Rainf==-9999.0)
           metin%Rainf=100.0
           metin%Rainf_flag=1
        END WHERE
        WHERE(metin%SWdown==-9999.0)
           metin%SWdown=1000.0
           metin%SWdown_flag=1
        END WHERE
        WHERE(metin%Tair==-9999.0)
           metin%Tair=35.0
           metin%Tair_flag=1
        END WHERE
        WHERE(metin%Qair==-9999.0)
           metin%Qair=0.0
           metin%Qair_flag=1
        END WHERE
        WHERE(metin%Wind==-9999.0)
           metin%Wind=0.0
           metin%Wind_flag=1
        END WHERE
     CASE('Gunnarsholt','gunnarsholt')
        WHERE(metin%Rainf==-9999.0)
           metin%Rainf=100.0
           metin%Rainf_flag=1
        END WHERE
        WHERE(metin%SWdown==-9999.0)
           metin%SWdown=1000.0
           metin%SWdown_flag=1
        END WHERE
        WHERE(metin%Tair==-9999.0)
           metin%Tair=35.0
           metin%Tair_flag=1
        END WHERE
        WHERE(metin%Qair==-9999.0)
           metin%Qair=0.0
           metin%Qair_flag=1
        END WHERE
        WHERE(metin%Wind==-9999.0)
           metin%Wind=0.0
           metin%Wind_flag=1
        END WHERE
     CASE('Blodgett','blodget','Blodgett Forest','blodgett forest')
        WHERE(metin%Rainf==-9999.0)
           metin%Rainf=100.0
           metin%Rainf_flag=1
        END WHERE
        WHERE(metin%SWdown==-9999.0)
           metin%SWdown=1000.0
           metin%SWdown_flag=1
        END WHERE
        WHERE(metin%Tair==-9999.0)
           metin%Tair=35.0
           metin%Tair_flag=1
        END WHERE
        WHERE(metin%Qair==-9999.0)
           metin%Qair=0.0
           metin%Qair_flag=1
        END WHERE
        WHERE(metin%Wind==-9999.0)
           metin%Wind=0.0
           metin%Wind_flag=1
        END WHERE
     CASE('Bondville','bondville')
        IF(start_time=='1997-01-01 00:01:00') THEN
!           metin%PSurf(1,1,1:17520) = metin%PSurf(1,1,17521:35040)
           OPEN(UNIT=33,FILE='BV98n_hh.met',FORM='formatted')
           READ (33,*)
           READ (33,*)
           DO bpcount = 1, 17520
             READ (33,*) skipi, skipr, skipr, skipi, skipr, skipi, &
                  & skipr, skipi, skipr, skipi, skipr, skipi, skipr, skipi, &
                  & skipr, skipi, skipr, skipi, skipr, skipi, skipr, skipi,&
                  & skipr, skipi, metin%PSurf(1,1,bpcount)
           END DO
           CLOSE(33)
           metin%PSurf_flag(1,1,1:17520) = 1
           WRITE(*,*) '    Year 1997 PSurf replaced. Change flagged' 
        END IF
     CASE DEFAULT
        WRITE(*,*) 'Site not known; no fixes to dataset.'
     END SELECT
  END IF
!------------------------END SITE DEPENDENT FIXES--------------------------------

! All site changes:

  WHERE(metin%Rainf_flag/=0)
     flstin%NEE_flag=5
     flstin%Qle_flag=5
     flstin%Qh_flag=5
  END WHERE


  !------------------------------UNITS CHANGES--------------------------------------
  ! Check for ridiculous values in essential variables (pre units change):
  IF(ANY(metin%Rainf<0.0.OR.metin%Rainf>300.0)) STOP "undefined values of Rainf!"
  IF(ANY(metin%SWdown<0.0.OR.metin%SWdown>1360.0)) STOP "undefined values of SWdown!"
  IF(ANY(metin%Tair<-100.0.OR.metin%Tair>80.0)) STOP "undefined values of Tair!"
  IF(ANY(metin%Qair<0.0.OR.metin%Qair>100.0)) STOP "undefined values of Qair!"
  IF(ANY(metin%Wind<0.0.OR.metin%Wind>100.0)) STOP "undefined values of Wind!"

  ! Any units corrections to data:
  metin%Tair=metin%Tair+273.15 ! C to K
  flstin%SoilTemp=flstin%SoilTemp+273.15 ! C to K

  IF(timestepsize==3600.0)  metin%Rainf=2*metin%Rainf ! from mm/timestep(halfhour) to mm/h
  metin%Rainf=metin%Rainf/timestepsize ! convert from mm/h to mm/s
  metin%Snowf=0
  metin%Snowf_flag=5
  DO i=1,mtmid ! synthesize longwave with Swinbank formula
     metin%LWdown(1,1,i)=0.0000094*0.0000000567*(metin%Tair(1,1,1,i)**6.0)
     metin%LWdown_flag(1,1,i) = 1 ! ie "missing from original"
  END DO
  ! Correct for missing CO2 values:
  IF(ANY(metin%CO2air==-9999.0.OR.metin%CO2air==9999.0)) &
       WRITE(*,*) 'Some undef CO2Air: using 350 ppm for these...'
  WHERE(metin%CO2air==-9999.0.OR.metin%CO2air==9999.0)
     metin%CO2air=350.0
     metin%CO2air_flag = 1 ! ie "missing from original"
  END WHERE
  ! Correct for missing PSurf values:
  IF(ANY(metin%PSurf==-9999.0.OR.metin%PSurf==9999.0.OR.metin%PSurf==0.0)) &
       WRITE(*,*) 'Some undef PSurf: using PSurf based on elevation and temp for these...'
  WHERE(metin%PSurf==-9999.0.OR.metin%PSurf==9999.0.OR.metin%PSurf==0.0)
     metin%PSurf=101.325*(metin%Tair(1,:,:,:)/(metin%Tair(1,:,:,:)+0.0065*elevation)) &
          **(9.80665/287.04/0.0065)
     metin%PSurf_flag = 1 ! ie "missing from original"
  END WHERE
  metin%PSurf=metin%PSurf*1000.0 ! Convert from kPa to Pa

  ! Convert relative to specific humididty:
  DO o = 1, SIZE(metin%Qair)
     CALL rh_sh(metin%Qair(1,1,1,o), metin%Tair(1,1,1,o), &
          metin%PSurf(1,1,o)/100.0,sph) ! PSurf: in Pa; conv to hPa
     metin%Qair(1,1,1,o) = sph
  END DO

  ! Again, check no ridiculous values (after units change):
  IF(ANY(metin%Rainf<0.0.OR.metin%Rainf>300.0/timestepsize)) &
    & STOP "undefined values of Rainf!"
  IF(ANY(metin%CO2air<0.0.OR.metin%CO2air>1000.0)) &
    & STOP "undefined values of CO2air!"
  IF(ANY(metin%LWdown<0.0.OR.metin%LWdown>600.0)) &
    & STOP "undefined values of LWdown!"
  IF(ANY(metin%SWdown<0.0.OR.metin%SWdown>1360.0)) &
    & STOP "undefined values of SWdown!"
  IF(ANY(metin%Tair<200.0.OR.metin%Tair>400.0)) STOP "undefined values of Tair!"
  IF(ANY(metin%Qair<0.0.OR.metin%Qair>0.04)) STOP "undefined values of Qair!"
  IF(ANY(metin%PSurf<30000.0.OR.metin%PSurf>120000.0)) &
    & STOP "undefined values of PSurf!"
  IF(ANY(metin%Wind<0.0.OR.metin%Wind>100.0)) STOP "undefined values of Wind!"

  !---------------------------------END UNITS CHANGES-----------------------------------

  WRITE(*,*) 'Met variables created;'

  finclude=.TRUE. ! initialise

  ! Decide which output varibles are available:
  WHERE(flstin%Qle<-800.0.OR.flstin%Qle>1000.0)
    flstin%Qle=0.0
    flstin%Qle_flag=1
  END WHERE
  WHERE(flstin%Qh<-800.0.OR.flstin%Qh>1000.0)
     flstin%Qh=0.0
    flstin%Qh_flag=1
  END WHERE
  WHERE(flstin%Qg<-800.0.OR.flstin%Qg>1000.0)
    flstin%Qg=0.0
    flstin%Qg_flag=1
  END WHERE
  WHERE(flstin%Rnet<-800.0.OR.flstin%Rnet>1100.0)
     flstin%Rnet=0.0
     flstin%Rnet_flag=1
  END WHERE
  WHERE(flstin%NEE<-100.0.OR.flstin%NEE>100.0)
     flstin%NEE=0.0
     flstin%NEE_flag=1
  END WHERE
  WHERE(flstin%SoilTemp<-203.0.OR.flstin%SoilTemp>363.0)
     flstin%SoilTemp=0.0
     flstin%SoilTemp_flag=1
  END WHERE

!----------------------------------DATA COMPRESSION------------------------
  ! Begin compressing data if required:
  IF(compresstype=='half') THEN
     IF(MOD(mtmid,2)==1) THEN
        WRITE(*,*) 'Cannot halve dataset with an odd number of timesteps.', &
             'Try changing starting time.'
        STOP
     END IF
     mtout=mtmid/2
     ! allocate memory for new half-size variables:
     ALLOCATE(met2%SWdown(1,1,mtout), met2%LWdown(1,1,mtout), &
          met2%Rainf(1,1,mtout), met2%Snowf(1,1,mtout), &
          met2%PSurf(1,1,mtout), met2%CO2air(1,1,za,mtout), &
          met2%Tair(1,1,za,mtout), met2%Wind(1,1,za,mtout), &
          met2%Qair(1,1,za,mtout), flst2%Qle(1,1,mtout), &
          flst2%Qh(1,1,mtout), flst2%Qg(1,1,mtout), &
          flst2%NEE(1,1,mtout), flst2%Rnet(1,1,mtout), &
          flst2%SoilTemp(1,1,zs,mtout),met2%LAI(1,1,mtout))
     ALLOCATE(time2(mtout),timestp2(mtout))
     ! Allocate memory for new half size met/flux variable flags:
     ALLOCATE(met2%SWdown_flag(1,1,mtout),met2%LWdown_flag(1,1,mtout), &
          met2%Rainf_flag(1,1,mtout),met2%Snowf_flag(1,1,mtout), &
          met2%PSurf_flag(1,1,mtout),met2%CO2air_flag(1,1,za,mtout), &
          met2%Tair_flag(1,1,za,mtout),met2%Wind_flag(1,1,za,mtout), &
          met2%Qair_flag(1,1,za,mtout), flst2%Qle_flag(1,1,mtout), &
          flst2%Qh_flag(1,1,mtout),flst2%Qg_flag(1,1,mtout), &
          flst2%NEE_flag(1,1,mtout),flst2%Rnet_flag(1,1,mtout), &
          flst2%SoilTemp_flag(1,1,zs,mtout), met2%LAI_flag(1,1,mtout))

     DO m=1,mtout !-1
        IF(metin%SWdown(1,1,m*2)==metin%SWdown(1,1,m*2+1)) THEN
           met2%SWdown(1,1,m)=metin%SWdown(1,1,m*2)
           met2%SWdown_flag(1,1,m)=metin%SWdown_flag(1,1,m*2)
        ELSE
           met2%SWdown(1,1,m)=0.0
        END IF
        met2%LWdown(1,1,m)=metin%LWdown(1,1,m*2)
        met2%Rainf(1,1,m)=metin%Rainf(1,1,m*2)
        met2%Snowf(1,1,m)=metin%Snowf(1,1,m*2)
        met2%PSurf(1,1,m)=metin%PSurf(1,1,m*2)
        IF(.NOT.fixedlai) met2%LAI(1,1,m)=metin%LAI(1,1,m*2)
        met2%CO2air(1,1,:,m)=metin%CO2air(1,1,:,m*2)
        met2%Tair(1,1,:,m)=metin%Tair(1,1,:,m*2)
        met2%Wind(1,1,:,m)=metin%Wind(1,1,:,m*2)
        met2%Qair(1,1,:,m)=metin%Qair(1,1,:,m*2)
        flst2%Qle(1,1,m)=flstin%Qle(1,1,m*2)
        flst2%Qh(1,1,m)=flstin%Qh(1,1,m*2)
        flst2%Qg(1,1,m)=flstin%Qg(1,1,m*2)
        flst2%NEE(1,1,m)=flstin%NEE(1,1,m*2)
        flst2%Rnet(1,1,m)=flstin%Rnet(1,1,m*2)
        flst2%SoilTemp(1,1,:,m)=flstin%SoilTemp(1,1,:,m*2)
        time2(m)=timein(m)
        timestp2(m)=timestpin(m)
        met2%LWdown_flag(1,1,m)=metin%LWdown_flag(1,1,m*2)
        met2%Rainf_flag(1,1,m)=metin%Rainf_flag(1,1,m*2)
        met2%Snowf_flag(1,1,m)=metin%Snowf_flag(1,1,m*2)
        met2%PSurf_flag(1,1,m)=metin%PSurf_flag(1,1,m*2)
        IF(.NOT.fixedlai) met2%LAI_flag(1,1,m)=metin%LAI_flag(1,1,m*2)
        met2%CO2air_flag(1,1,:,m)=metin%CO2air_flag(1,1,:,m*2)
        met2%Tair_flag(1,1,:,m)=metin%Tair_flag(1,1,:,m*2)
        met2%Wind_flag(1,1,:,m)=metin%Wind_flag(1,1,:,m*2)
        met2%Qair_flag(1,1,:,m)=metin%Qair_flag(1,1,:,m*2)
        flst2%Qle_flag(1,1,m)=flstin%Qle_flag(1,1,m*2)
        flst2%Qh_flag(1,1,m)=flstin%Qh_flag(1,1,m*2)
        flst2%Qg_flag(1,1,m)=flstin%Qg_flag(1,1,m*2)
        flst2%NEE_flag(1,1,m)=flstin%NEE_flag(1,1,m*2)
        flst2%Rnet_flag(1,1,m)=flstin%Rnet_flag(1,1,m*2)
        flst2%SoilTemp_flag(1,1,:,m)=flstin%SoilTemp_flag(1,1,:,m*2)
     END DO
     flst=>flst2
     met=>met2
     time=>time2
     timestp=>timestp2
  ELSE IF(compresstype=='avgd') THEN
     IF(MOD(mtmid,2)==1) THEN
        WRITE(*,*) 'Cannot halve dataset with an odd number of timesteps.', &
             'Try changing starting time.'
        STOP
     END IF
      mtout=mtmid/2
     ! allocate memory for new half-size variables:
     ALLOCATE(met2%SWdown(1,1,mtout),met2%LWdown(1,1,mtout), &
          met2%Rainf(1,1,mtout),met2%Snowf(1,1,mtout), &
          met2%PSurf(1,1,mtout),met2%CO2air(1,1,za,mtout), &
          met2%Tair(1,1,za,mtout),met2%Wind(1,1,za,mtout), &
          met2%Qair(1,1,za,mtout),flst2%Qle(1,1,mtout), &
          flst2%Qh(1,1,mtout),flst2%Qg(1,1,mtout), &
          flst2%NEE(1,1,mtout),flst2%Rnet(1,1,mtout), &
          flst2%SoilTemp(1,1,zs,mtout),met2%LAI(1,1,mtout))
     ALLOCATE(time2(mtout),timestp2(mtout))
     ! Allocate memory for new half size met/flux variable flags:
     ALLOCATE(met2%SWdown_flag(1,1,mtout),met2%LWdown_flag(1,1,mtout), &
          met2%Rainf_flag(1,1,mtout),met2%Snowf_flag(1,1,mtout), &
          met2%PSurf_flag(1,1,mtout),met2%CO2air_flag(1,1,za,mtout), &
          met2%Tair_flag(1,1,za,mtout),met2%Wind_flag(1,1,za,mtout), &
          met2%Qair_flag(1,1,za,mtout), flst2%Qle_flag(1,1,mtout), &
          flst2%Qh_flag(1,1,mtout),flst2%Qg_flag(1,1,mtout), &
          flst2%NEE_flag(1,1,mtout),flst2%Rnet_flag(1,1,mtout), &
          flst2%SoilTemp_flag(1,1,zs,mtout), met2%LAI_flag(1,1,mtout))
     DO m=1,mtout !-1
        met2%SWdown(1,1,m)=(metin%SWdown(1,1,m*2)+metin%SWdown(1,1,m*2-1))/2.0
        met2%LWdown(1,1,m)=(metin%LWdown(1,1,m*2)+metin%LWdown(1,1,m*2-1))/2.0
        ! fluxnet units for precip (mm/tstepsize) already 
        ! converted to mm/h if timestepsize=3600
        met2%Rainf(1,1,m)=(metin%Rainf(1,1,m*2)+metin%Rainf(1,1,m*2-1))/2.0
        met2%Snowf(1,1,m)=(metin%Snowf(1,1,m*2)+metin%Snowf(1,1,m*2-1))/2.0
        met2%PSurf(1,1,m)=(metin%PSurf(1,1,m*2)+metin%PSurf(1,1,m*2-1))/2.0
        IF(.NOT.fixedlai) met2%LAI(1,1,m)=(metin%LAI(1,1,m*2) &
                                         & +metin%LAI(1,1,m*2-1))/2.0
        met2%CO2air(1,1,:,m)=(metin%CO2air(1,1,:,m*2) &
                                         & +metin%CO2air(1,1,:,m*2-1))/2.0
        met2%Tair(1,1,:,m)=(metin%Tair(1,1,:,m*2)+metin%Tair(1,1,:,m*2-1))/2.0
        met2%Wind(1,1,:,m)=(metin%Wind(1,1,:,m*2)+metin%Wind(1,1,:,m*2-1))/2.0
        met2%Qair(1,1,:,m)=(metin%Qair(1,1,:,m*2)+metin%Qair(1,1,:,m*2-1))/2.0
        flst2%Qle(1,1,m)=(flstin%Qle(1,1,m*2)+flstin%Qle(1,1,m*2-1))/2.0
        flst2%Qh(1,1,m)=(flstin%Qh(1,1,m*2)+flstin%Qh(1,1,m*2-1))/2.0
        flst2%Qg(1,1,m)=(flstin%Qg(1,1,m*2)+flstin%Qg(1,1,m*2-1))/2.0
        flst2%NEE(1,1,m)=(flstin%NEE(1,1,m*2)+flstin%NEE(1,1,m*2-1))/2.0
        flst2%Rnet(1,1,m)=(flstin%Rnet(1,1,m*2)+flstin%Rnet(1,1,m*2-1))/2.0
        flst2%SoilTemp(1,1,:,m)=(flstin%SoilTemp(1,1,:,m*2) &
                                         & +flstin%SoilTemp(1,1,:,m*2))/2.0
        time2(m)=timein(m)
        timestp2(m)=timestpin(m)
        ! gapfilling flag data:
        met2%SWdown_flag(1,1,m)=MAX(metin%SWdown_flag(1,1,m*2), &
                                         & metin%SWdown_flag(1,1,m*2-1))
        met2%LWdown_flag(1,1,m)=MAX(metin%LWdown_flag(1,1,m*2), &
                                         & metin%LWdown_flag(1,1,m*2-1))
        met2%Rainf_flag(1,1,m)=MAX(metin%Rainf_flag(1,1,m*2), &
                                         & metin%Rainf_flag(1,1,m*2-1))
        met2%Snowf_flag(1,1,m)=MAX(metin%Snowf_flag(1,1,m*2), &
                                         & metin%Snowf_flag(1,1,m*2-1))
        met2%PSurf_flag(1,1,m)=MAX(metin%PSurf_flag(1,1,m*2), &
                                         & metin%PSurf_flag(1,1,m*2-1))
        IF(.NOT.fixedlai) met2%LAI_flag(1,1,m) &
                    & =MAX(metin%LAI_flag(1,1,m*2),metin%LAI_flag(1,1,m*2-1))
        met2%CO2air_flag(1,1,:,m)=MAX(metin%CO2air_flag(1,1,:,m*2), &
                                         & metin%CO2air_flag(1,1,:,m*2-1))
        met2%Tair_flag(1,1,:,m)=MAX(metin%Tair_flag(1,1,:,m*2), &
                                         & metin%Tair_flag(1,1,:,m*2-1))
        met2%Wind_flag(1,1,:,m)=MAX(metin%Wind_flag(1,1,:,m*2), &
                                         & metin%Wind_flag(1,1,:,m*2-1))
        met2%Qair_flag(1,1,:,m)=MAX(metin%Qair_flag(1,1,:,m*2), &
                                         & metin%Qair_flag(1,1,:,m*2-1))
        flst2%Qle_flag(1,1,m)=MAX(flstin%Qle_flag(1,1,m*2), &
                                         & flstin%Qle_flag(1,1,m*2-1))
        flst2%Qh_flag(1,1,m)=MAX(flstin%Qh_flag(1,1,m*2), &
                                         & flstin%Qh_flag(1,1,m*2-1))
        flst2%Qg_flag(1,1,m)=MAX(flstin%Qg_flag(1,1,m*2), &
                                         & flstin%Qg_flag(1,1,m*2-1))
        flst2%NEE_flag(1,1,m)=MAX(flstin%NEE_flag(1,1,m*2), &
                                         & flstin%NEE_flag(1,1,m*2-1))
        flst2%Rnet_flag(1,1,m)=MAX(flstin%Rnet_flag(1,1,m*2), &
                                         & flstin%Rnet_flag(1,1,m*2-1))
        flst2%SoilTemp_flag(1,1,:,m)=MAX(flstin%SoilTemp_flag(1,1,:,m*2), &
                                         & flstin%SoilTemp_flag(1,1,:,m*2))
     END DO
     flst=>flst2
     met=>met2
     time=>time2
     timestp=>timestp2
  ELSE
     flst=>flstin
     met=>metin
     time=>timein
     timestp=>timestpin
     mtout=mtmid
  END IF

  ! Set netcdf counters:
  start_out1=1
  count_out1=mtout
  start_out2=1
  count_out2=1
  start_out3=1
  count_out3=1
  count_out3(3)=mtout
  start_out4=1
  count_out4a=1
  count_out4a(3)=za
  count_out4a(4)=mtout
  count_out4s=1
  count_out4s(3)=zs
  count_out4s(4)=mtout

  !--------------------- BEGIN WRITING TO NETCDF FILE ---------------------------
  ! Create NetCDF file:
  status = NF90_create(fileout, NF90_CLOBBER, ncid)
  IF (status /= NF90_noerr) CALL handle_err(status)

  ! Put the file in define mode:
  status = NF90_redef(ncid)                   

  ! Define dimensions:
  status = NF90_def_dim(ncid, 'x', 1, xID)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_def_dim(ncid, 'y', 1, yID) 
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_def_dim(ncid, 'z', za, zID)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_def_dim(ncid, 'z_down', zs, zsID)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_def_dim(ncid, "time", NF90_unlimited, tID)
  IF (status /= NF90_noerr) CALL handle_err(status)

  ! Define dimension variables:
  status = NF90_def_var(ncid,'nav_lon',NF90_float,(/yID,xID/),nav_lonID)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,nav_lonID,'long_name', 'Longitude')
  status = NF90_PUT_ATT(ncid,nav_lonID,'units','degrees_east')
  status = NF90_def_var(ncid,'nav_lat' ,NF90_float,(/yID,xID/),nav_latID)
  status = NF90_PUT_ATT(ncid,nav_latID, "long_name", "Latitude")
  status = NF90_PUT_ATT(ncid,nav_latID, "units","degrees_north" )
  status = NF90_def_var(ncid,'x',NF90_int,(/xID/),xvID)
  status = NF90_def_var(ncid,'y',NF90_int,(/yID/),yvID)
  status = NF90_def_var(ncid,'elevation',NF90_float,(/zID/),levelID)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_def_var(ncid,'time',NF90_double,(/tID/),timeID)
  status = NF90_PUT_ATT(ncid,timeID,'long_name', 'Time')
  status = NF90_PUT_ATT(ncid,timeID,'units',timeunits)
  status = NF90_PUT_ATT(ncid,timeID,'calendar', 'gregorian')
  status = NF90_PUT_ATT(ncid,timeID,'time_origin', timeorigin)
  status = NF90_def_var(ncid,'timestp',NF90_float,(/tID/),timestpID)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,timestpID,'long_name', 'Timesteps')
  status = NF90_PUT_ATT(ncid,timestpID,'units',timestpunits)
  status = NF90_PUT_ATT(ncid,timestpID,'tstep_sec',REAL(timestepsize))

  ! Define met variables:
  status = NF90_def_var(ncid, "SWdown", NF90_float,  &
       (/ xID, yID, tID /), SWdownID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,SWdownID,"units",units%SWdown )
  status = NF90_PUT_ATT(ncid,SWdownID,"long_name", &
                       & 'Surface incident shortwave radiation' )
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_def_var(ncid, "LWdown", NF90_float,  &
       (/ xID, yID, tID /), LWdownID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,LWdownID,"units",units%LWdown)
  status = NF90_PUT_ATT(ncid,LWdownID, "long_name", &
                       & "Surface incident longwave radiation")
  status = NF90_PUT_ATT(ncid,LWdownID, "source", &
                       & "Downward longwave from Swinbank formula")
  status = NF90_def_var(ncid, "Rainf", NF90_float, &
       (/ xID, yID, tID /), RainfID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,RainfID, "units",units%Rainf )
  status = NF90_PUT_ATT(ncid,RainfID, "long_name",'Rainfall rate')
  status = NF90_def_var(ncid, "Snowf", NF90_float, &
       (/ xID, yID, tID /), SnowfID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,SnowfID, "units",units%Snowf )
  status = NF90_PUT_ATT(ncid,SnowfID, "long_name",'Snowfall rate')
  status = NF90_def_var(ncid, "Tair", NF90_float,  &
       (/ xID, yID, zID, tID /), TairID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,TairID, "units",units%Tair )
  status = NF90_PUT_ATT(ncid,TairID, "long_name",'Near surface air temperature')
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_def_var(ncid, "PSurf", NF90_float, &
       (/ xID, yID, tID /), PSurfID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,PSurfID, "units",units%PSurf)
  status = NF90_PUT_ATT(ncid,PSurfID, "long_name",'Surface pressure')
  status = NF90_def_var(ncid, "Wind", NF90_float, &
       (/ xID, yID, zID, tID /), WindID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,WindID,"units",units%Wind)
  status = NF90_PUT_ATT(ncid,WindID,"long_name", &
                       & 'Near surface module of the wind')
  status = NF90_def_var(ncid, "Qair", NF90_float,  &
       (/ xID, yID, zID, tID /), QairID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,QairID, "units",units%Qair)
  status = NF90_PUT_ATT(ncid,QairID, "long_name", &
                       & 'Near surface specific humidity')
  status = NF90_def_var(ncid, "CO2air", NF90_float, &
       (/ xID, yID, zID, tID /), CO2airID)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,CO2airID, "units",units%CO2air)
  status = NF90_PUT_ATT(ncid,CO2airID, "long_name", &
                       & 'Near surface CO2 concentration')
  IF(.NOT.fixedlai) THEN
     status = NF90_def_var(ncid, "LAI", NF90_float, &
          (/ xID, yID, tID /), LAIID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,LAIID, "units",units%LAI)
     status = NF90_PUT_ATT(ncid,LAIID, "long_name",'Average Leaf Area Index')
     status = NF90_PUT_ATT(ncid,LaiID, "source", &
            & "data from Junhua Yan and Guoyi Zhou")    ! BP
  END IF

  ! Define flux/state variables:
  IF(finclude(1)) THEN
     status = NF90_def_var(ncid, "Qle", NF90_float, &
          (/ xID, yID, tID /), QleID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,QleID, "units",units%Qle)
     status = NF90_PUT_ATT(ncid,QleID, "long_name",'Latent heat flux')
     WRITE(*,*) 'Qle created;'
  END IF
  IF(finclude(2)) THEN
     status = NF90_def_var(ncid, "Qh", NF90_float, &
          (/ xID, yID, tID /), QhID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,QhID, "units",units%Qh )
     status = NF90_PUT_ATT(ncid,QhID, "long_name",'Sensible heat flux')
     WRITE(*,*) 'Qh created;'
  END IF
  IF(finclude(3)) THEN
     status = NF90_def_var(ncid, "Qg", NF90_float,  &
          (/ xID, yID, tID /), QgID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,QgID, "units",units%Qg)
     status = NF90_PUT_ATT(ncid,QgID, "long_name",'Ground heat flux')
     WRITE(*,*) 'Qg created;'
  END IF
  IF(finclude(4)) THEN
     status = NF90_def_var(ncid, "Rnet", NF90_float, &
          (/ xID, yID, tID /), RnetID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,RnetID, "units",units%Rnet)
     status = NF90_PUT_ATT(ncid,RnetID, "long_name",'Net absorbed radiation')
     WRITE(*,*) 'Rnet created;'
  END IF
  IF(finclude(5)) THEN
     status = NF90_def_var(ncid, "NEE", NF90_float, &
          (/ xID, yID, tID /), NEEID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,NEEID, "units",units%NEE)
     status = NF90_PUT_ATT(ncid,NEEID, "long_name",'Net Ecosystem Exchange')
     WRITE(*,*) 'NEE created;'
  END IF
  IF(finclude(6)) THEN
     status = NF90_def_var(ncid, "SoilTemp", NF90_float, &
          (/ xID, yID, zsID, tID /), SoilTempID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_ATT(ncid,SoilTempID, "units",units%SoilTemp)
     status = NF90_PUT_ATT(ncid,SoilTempID, "long_name", &
                       & 'Average layer soil temperature')
     WRITE(*,*) 'SoilTemp created;'
  END IF

  ! Define parameters:
  ! isoil
  status=NF90_DEF_VAR(ncid,'isoil',NF90_FLOAT,(/xID,yID/),isoilID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,isoilID,"long_name",&
       "Soil type from Zobler")
  status = NF90_PUT_ATT(ncid,isoilID,"units","-")
  ! iveg
  status=NF90_DEF_VAR(ncid,'iveg',NF90_FLOAT,(/xID,yID/),ivegID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,ivegID,"long_name",&
       "Vegetation type from Potter et al")
  status = NF90_PUT_ATT(ncid,ivegID,"units","-")
  ! frac4 
  status=NF90_DEF_VAR(ncid,'frac4',NF90_FLOAT,(/xID,yID/),frac4ID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,frac4ID,"long_name",&
       "Fraction of vegetation which is c4")
  status = NF90_PUT_ATT(ncid,frac4ID,"units","-")
  ! hc 
  status=NF90_DEF_VAR(ncid,'hc',NF90_FLOAT,(/xID,yID/),hcID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,hcID,"long_name",&
       "Vegetation height")
  status = NF90_PUT_ATT(ncid,hcID,"units","m")
  ! za 
  status=NF90_DEF_VAR(ncid,'za',NF90_FLOAT,(/xID,yID/),zaID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,zaID,"long_name",&
       "Reference height")
  status = NF90_PUT_ATT(ncid,zaID,"units","m")

  ! Write variable flags:
  IF(flag) THEN
     WRITE(*,*) 'Writing variable flags;'
     ! Define met variable flags:
     status = NF90_def_var(ncid, "SWdown_flag", NF90_int, &
                       & (/ xID, yID, tID /), SWdownfID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "LWdown_flag", NF90_int, &
                       & (/ xID, yID, tID /), LWdownfID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "Rainf_flag", NF90_int, &
                       & (/ xID, yID, tID /), RainffID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "Snowf_flag", NF90_int, &
                       & (/ xID, yID, tID /), SnowffID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "Tair_flag", NF90_int,  &
          (/ xID, yID, zID, tID /), TairfID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "PSurf_flag", NF90_int, &
                       & (/ xID, yID, tID /), PSurffID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "Wind_flag", NF90_int,  &
          (/ xID, yID, zID, tID /), WindfID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "Qair_flag", NF90_int, &
          (/ xID, yID, zID, tID /), QairfID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_def_var(ncid, "CO2air_flag", NF90_int, &
          (/ xID, yID, zID, tID /), CO2airfID)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(.NOT.fixedlai) THEN
        status = NF90_def_var(ncid, "LAI_flag", NF90_int, &
             (/ xID, yID, tID /), LAIfID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     ! Define flux/state variable flags:
     IF(finclude(1)) THEN
        status = NF90_def_var(ncid,"Qle_flag",NF90_int, &
                             & (/ xID, yID, tID /),QlefID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     IF(finclude(2)) THEN
        status = NF90_def_var(ncid,"Qh_flag",NF90_int,(/ xID, yID, tID /),QhfID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     IF(finclude(3)) THEN
        status = NF90_def_var(ncid,"Qg_flag", NF90_int, &
                             & (/ xID, yID, tID /),QgfID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     IF(finclude(4)) THEN
        status = NF90_def_var(ncid,"Rnet_flag",NF90_int, &
                             & (/ xID, yID, tID /),RnetfID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     IF(finclude(5)) THEN
        status = NF90_def_var(ncid,"NEE_flag",NF90_int, &
                             & (/ xID, yID, tID /),NEEfID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     IF(finclude(6)) THEN
        status = NF90_def_var(ncid, "SoilTemp_flag", NF90_int, &
             (/ xID, yID, zsID, tID /), SoilTempfID)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
  END IF

  ! Write time of production to netcdf file:
  CALL DATE_AND_TIME(todaydate, nowtime)
  todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
  nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Production", &
       TRIM(todaydate)//' at '//TRIM(nowtime))
  ! Other global details:
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Site_name",sitename)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Length_of_dataset",datalength)
  IF (status /= NF90_noerr) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Contact",contact)

  ! End define mode:
  status = NF90_enddef(ncid) 
  IF (status /= NF90_noerr) CALL handle_err(status)

  ! Write variables:
  status = NF90_PUT_VAR(ncid, nav_lonID, longitude)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, nav_latID,  latitude)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, xvID, xv)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, yvID,  yv)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, levelID, elevation)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, timeID,  time)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, timestpID,  timestp)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, SWdownID,  met%SWdown)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, LWdownID,  met%LWdown)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, RainfID,  met%Rainf)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, SnowfID,  met%Snowf)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, TairID,  met%Tair)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, PSurfID, met%PSurf)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, WindID,  met%Wind)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, QairID,  met%Qair)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, CO2airID,  met%CO2air)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  IF(.NOT.fixedlai) THEN
     status = NF90_PUT_VAR(ncid, LAIID, met%LAI)
     IF(status /= NF90_NoErr) CALL handle_err(status)
  END IF
  IF(finclude(1)) status = NF90_PUT_VAR(ncid, QleID, flst%Qle)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  IF(finclude(2)) status = NF90_PUT_VAR(ncid, QhID,  flst%Qh)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  IF(finclude(3)) status = NF90_PUT_VAR(ncid, QgID,  flst%Qg)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  IF(finclude(4)) status = NF90_PUT_VAR(ncid, RnetID, flst%Rnet)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  IF(finclude(5)) status = NF90_PUT_VAR(ncid, NEEID, flst%NEE)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  IF(finclude(6)) status = NF90_PUT_VAR(ncid, SoilTempID,  flst%SoilTemp)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, isoilID,par%isoil)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, ivegID,par%iveg)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, frac4ID,par%frac4)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, hcID,par%hc)
  IF(status /= NF90_NoErr) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, zaID,par%za)
  IF(status /= NF90_NoErr) CALL handle_err(status)

  IF(flag) THEN
     ! write flag vriables:
     status = NF90_PUT_VAR(ncid, SWdownfID, met%SWdown_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, LWdownfID,  met%LWdown_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, RainffID,  met%Rainf_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, SnowffID,  met%Snowf_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, TairfID,  met%Tair_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, PSurffID,  met%PSurf_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, WindfID, met%Wind_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, QairfID,  met%Qair_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     status = NF90_PUT_VAR(ncid, CO2airfID,  met%CO2air_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(.NOT.fixedlai) THEN
        status = NF90_PUT_VAR(ncid, LAIfID,  met%LAI_flag)
        IF(status /= NF90_NoErr) CALL handle_err(status)
     END IF
     IF(finclude(1)) status = NF90_PUT_VAR(ncid, QlefID, flst%Qle_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(finclude(2)) status = NF90_PUT_VAR(ncid, QhfID,flst%Qh_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(finclude(3)) status = NF90_PUT_VAR(ncid, QgfID, flst%Qg_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(finclude(4)) status = NF90_PUT_VAR(ncid, RnetfID, flst%Rnet_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(finclude(5)) status = NF90_PUT_VAR(ncid, NEEfID, flst%NEE_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
     IF(finclude(6)) status = NF90_PUT_VAR(ncid, SoilTempfID,flst%SoilTemp_flag)
     IF(status /= NF90_NoErr) CALL handle_err(status)
  END IF

  ! Close NetCDF file:
  status = NF90_close(ncid) 
  IF (status /= NF90_noerr) CALL handle_err(status)

  !------------------------------ END OUTPUT ----------------------------------

  WRITE(*,*) 'CHECK netcdf header for lat, lon, elevation and site details.'

CONTAINS
  SUBROUTINE handle_err(status)
    INTEGER, INTENT(IN) :: status
    IF(status /= NF90_noerr) THEN
       PRINT*, TRIM(NF90_strerror(status))
       STOP "Stopped"
    END IF
  END SUBROUTINE handle_err

  SUBROUTINE rh_sh (rh, tl, pl,        & ! I
       sh)          ! O
    !
    !    Purpose - To convert relative humidity to specific humidity
    !
    !
    ! !Input Parameters:
    !   pl - pressure level (hPa)
    !   rh - relative humidity value [0-100]
    !   tl - temperature at pressure level (K)
    !
    ! !Output Parameters:
    !   sh - specific humidity value (g/g)
    !
    !  Internal Routines Called:
    !    savtap
    !
    ! !Revision History:
    !  Version 2.0   August 15, 1997
    !                Lisa Coleman, SAIC, l.h.coleman@larc.nasa.gov
    !                Shalini Gupta, SAIC, s.gupta@larc.nasa.gov
    !                Fred Rose, AS&M, f.g.rose@larc.nasa.gov
    !                Initial Version of Code
    !

    REAL, INTENT (IN)  :: pl
    REAL, INTENT (IN)  :: rh
    REAL, INTENT (OUT) :: sh
    REAL, INTENT (IN)  :: tl
    REAL :: es

    REAL :: ws

    es = satvap (tl)
    ws = 0.622 * es / (pl - es)
    sh = (rh/100.0) * ws
  END SUBROUTINE rh_sh
  !----------------------------------------------------------------------------

  FUNCTION satvap (tl) RESULT (F_Result)

    !    Purpose - To compute the saturation value used in converting between
    !              specific humidity and relative humidity
    !
    !
    ! !Input Parameters:
    !   tl - temperature at pressure level being converted
    !
    ! !Output Parameters:
    !   F_Result - saturation value
    !
    ! !Revision History:
    !  Version 2.0   August 15, 1997
    !                Lisa Coleman, SAIC, l.h.coleman@larc.nasa.gov
    !                Shalini Gupta, SAIC, s.gupta@larc.nasa.gov
    !                Fred Rose, AS&M, f.g.rose@larc.nasa.gov
    !                Initial Version of Code

    REAL :: eilog
    REAL :: ewlog, ewlog2, ewlog3, ewlog4
    REAL :: F_Result
    REAL :: temp, tl
    REAL :: toot, toto, tsot
    temp = tl - 273.155
    IF (temp < -20.) THEN
       !     *** ice saturation
       toot = 273.16 / tl
       toto = 1. / toot
       eilog =   -9.09718 * (toot-1)                                      &
            -  3.56654 * (LOG (toot) / LOG (10.))                      &
            +  .876793 * (1-toto)                                      &
            +  (LOG (6.1071) / LOG (10.))
       F_Result = 10.**eilog
    ELSE
       tsot = 373.16 / tl
       ewlog = -7.90298 * (tsot-1) + 5.02808 * (LOG (tsot) / LOG (10.))
       ewlog2 =   ewlog                                                   &
            - 1.3816e-07 * (10**(11.344 * (1 - (1/tsot))) - 1)
       ewlog3 =    ewlog2                                                 &
            + .0081328 * (10**(-3.49149 * (tsot-1)) - 1)
       ewlog4 = ewlog3 + (LOG (1013.246) / LOG (10.))
       F_Result = 10.**ewlog4
    END IF
  END FUNCTION satvap


END PROGRAM fluxnet_to_netcdf
