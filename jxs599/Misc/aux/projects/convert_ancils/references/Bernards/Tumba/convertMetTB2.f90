! Writes a NetCDF binary file from text met forcing, essentially following ALMA
! convention (some units changes) as well as being Grads readable.
! Modified by BP to write site specific LAI to the NetCDF file (Sep 2007).

PROGRAM texttonc
  USE netcdf
  IMPLICIT NONE
  ! Input/output details
  !======================================================================
  INTEGER, PARAMETER  :: mt = 56184   ! # met forcing timesteps      ! BP
  CHARACTER(LEN=*), PARAMETER :: filein = 'TumbaDataii.csv'          ! BP
  CHARACTER(LEN=*), PARAMETER :: fileout = 'Tumbarumba2001-07.nc'    ! BP
  CHARACTER(LEN=33), PARAMETER :: timeunits= &
                                  & 'seconds since 2001-01-01 00:00:01' ! BP
  CHARACTER(LEN=33), PARAMETER :: timestpunits= &
                                  & 'timesteps since 2001-02-21 23:59:00'
  INTEGER, PARAMETER :: timestepsize=3600 ! timestep size in seconds
  REAL :: latitude = -35.6557 ! site latitude
  REAL :: longitude = 148.152 ! site longitude
  REAL :: elevation = 1200.0  ! site elevation in metres
  CHARACTER(LEN=*), PARAMETER :: sitename='Tumbarumba, NSW Australia'
  CHARACTER(LEN=*), PARAMETER :: datalength='Six and a half years'
  CHARACTER(LEN=*), PARAMETER :: contact='Gab Abramowitz, gabsun@gmail.com, +BP'
  ! Units to write to netcdf file:
  TYPE units_type
     CHARACTER(LEN=5) :: SWdown='W/m^2'
     CHARACTER(LEN=5) :: LWdown='W/m^2'
     CHARACTER(LEN=5) :: Qle='W/m^2'
     CHARACTER(LEN=5) :: Qh='W/m^2'
     CHARACTER(LEN=5) :: Qg='W/m^2'
     CHARACTER(LEN=5) :: SWnet='W/m^2'
     CHARACTER(LEN=5) :: LWnet='W/m^2'
     CHARACTER(LEN=3) :: CO2air='ppm'
     CHARACTER(LEN=4) :: Rainf='mm/s' ! changed from mm/h to mm/s below
     CHARACTER(LEN=4) :: Snowf='mm/s' ! changed from mm/h to mm/s below
     CHARACTER(LEN=1) :: Tair='K'   ! changed from C to K below
     CHARACTER(LEN=2) :: PSurf='mb'
     CHARACTER(LEN=3) :: Wind='m/s'
     CHARACTER(LEN=5) :: Qair='kg/kg' ! changed from % to kg/kg below
     CHARACTER(LEN=10):: NEE='umol/m^2/s'    ! BP
     CHARACTER(LEN=5) :: Rnet='W/m^2'
     CHARACTER(LEN=4) :: LAI='unit'   ! BP
  END TYPE units_type
  ! Specify any known CABLE parameters (others will be read through default_pars): 
  TYPE par_type
     INTEGER :: isoil(1,1) = 2
     INTEGER :: iveg(1,1) = 1
     REAL :: albsoil(1,1) = 0.1
     REAL :: clay(1,1) = 0.24
     REAL :: frac4(1,1) = 0.0
     REAL :: hc(1,1) = 40.0
     REAL :: sand(1,1) = 0.48
     REAL :: silt(1,1) = 0.28
     REAL :: za(1,1) = 70.0
  END TYPE par_type
  TYPE(par_type) :: par
  !=============================================================
  INTEGER, PARAMETER  :: za = 1       ! # vertical atmospheric layers
  INTEGER, PARAMETER  :: zs = 2       ! # vertical soil layers
  INTEGER :: lineskip=3 ! Skip how many lines in text file   ! BP
  INTEGER :: xv = 1 ! grid x-coord #
  INTEGER :: yv = 1 ! grid y-coord #
  TYPE(units_type) :: units
  ! Define met type
  TYPE met_input_type ! All timestep met data
     REAL,   DIMENSION(1,1,za,mt) :: CO2air   ! CO2 concentration
     REAL,   DIMENSION(1,1,mt) :: SWdown  ! downward short-wave radiation (W/m2)
     REAL,   DIMENSION(1,1,mt) :: LWdown  ! downward long-wave radiation (W/m2)
     REAL,   DIMENSION(1,1,mt) :: Rainf  ! rainfall 
     REAL,   DIMENSION(1,1,mt) :: Snowf  ! snowfall 
     REAL,   DIMENSION(1,1,za,mt) :: Tair  ! surface air temperature 
     REAL,   DIMENSION(1,1,mt) :: PSurf  ! surface air pressure 
     REAL,   DIMENSION(1,1,za,mt) :: Wind   ! surface wind speed 
     REAL,   DIMENSION(1,1,za,mt) :: Qair   ! surface relative humidity
     REAL,   DIMENSION(1,1,mt) :: LAI       ! BP
  END TYPE met_input_type
  TYPE fluxstate_output_type
     REAL, DIMENSION(1,1,mt) :: Qle ! latent heat
     REAL, DIMENSION(1,1,mt) :: Qh  ! sensible heat
     REAL, DIMENSION(1,1,mt) :: NEE ! net ecosystem exchange
     REAL, DIMENSION(1,1,mt) :: Qg  ! ground heat flux
     REAL, DIMENSION(1,1,zs,mt) :: SoilTemp ! soil temperature
     REAL, DIMENSION(1,1,mt) :: Rnet ! net radiation
     REAL, DIMENSION(1,1,mt) :: SWnet ! net shortwave radiation
  END TYPE fluxstate_output_type
  TYPE (met_input_type)       :: met 
  TYPE (fluxstate_output_type) :: flst
  REAL(KIND(0.0D0)) :: time(mt) 
  REAL :: timestp(mt)
  REAL :: skip ! skip column in read
  ! I/o declarations
  LOGICAL :: eof	 ! end of file?  
  INTEGER :: i,j,o ! do counters
  INTEGER :: ios	 ! I/O status
  LOGICAL :: is_open
  INTEGER :: unitin = 12 ! input unitin (open & close here)
  ! Declarations for netcdf writing 
  INTEGER :: xID,yID,zID,zsID,tID  ! Dimension IDs
  INTEGER :: SWdownID, LWdownID, RainfID, SnowfID, TairID, PSurfID, WindID, & 
       QairID, LaiID, &        ! BP
       nav_lonID, nav_latID, xvID, yvID, levelID, timeID, timestpID, &
       QleID, QhID, NEEID, QgID, SoilTempID, RnetID, SWnetID, ivegID,isoilID,&
       albsoilID,clayID,frac4ID,hcID,sandID,siltID,zaID! Variable IDs
  INTEGER :: ncid, status
  REAL :: sph ! temporary specific humidity
  CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file

  REAL :: tempLAI(365)      ! BP
  INTEGER :: bp1, bp2, bp3, bp4, bp5    ! BP

  ! Open input text file:
  INQUIRE(UNIT = unitin, OPENED = is_open)
  IF (.NOT. is_open) THEN
     OPEN (UNIT=unitin, FILE=filein, STATUS='OLD', ACTION='READ', IOSTAT = ios)
     IF (ios /= 0) WRITE(*,*) 'Error opening file ', filein
  END IF

  ! Read from file and write to 'met' one timestep at a time:
  DO j=1, lineskip
     READ (unitin, *, IOSTAT = ios)
  END DO
  IF (ios /= 0) WRITE(*,*) 'Error reading file ', filein

  DO i = 1, mt                    
     READ (unitin, *, IOSTAT = ios) bp1, bp2, bp3, bp4, bp5,  &
          met%SWdown(1,1,i), skip, flst%Rnet(1,1,i), met%Rainf(1,1,i), &
          met%Tair(1,1,1,i), met%Wind(1,1,1,i), met%Qair(1,1,1,i), &
          met%PSurf(1,1,i), met%LWdown(1,1,i), flst%Qh(1,1,i), &
          flst%Qle(1,1,i), flst%NEE(1,1,i)
     timestp(i)=REAL(i)
!     IF (met%LWdown(1,1,i) < 0.0) THEN ! BP to fix negative data
!       PRINT *, i, met%LWdown(1,1,i), 'LWdown' ! BP for 2004Nov2
!       met%LWdown(1,1,i) = 0.0
!     END IF
     IF (flst%Rnet(1,1,i) < -3500.0) THEN ! BP to fix very negative data
       PRINT *, i, flst%Rnet(1,1,i), 'Rnet' ! BP for 2004Nov12 noon
       flst%Rnet(1,1,i) = 0.0
     END IF
     time(i)=timestp(i)*REAL(timestepsize)+4489200.0
     eof = i==mt.AND.ios==0
     IF (eof) THEN
        CLOSE(unitin)
     ELSE
        IF (ios /= 0) WRITE(*,*) 'Error reading file ', filein
     END IF
     IF (bp1==29 .and. bp2==2 .and. bp4==0) THEN   ! BP to get rid of leap day
!       print *, 'getting to a leap day.'
!       print *, bp1, bp2, bp3, bp4, bp5, met%SWdown(1,1,i), flst%Rnet(1,1,i), met%Rainf(1,1,i), flst%NEE(1,1,i)
!       print *, '29-2-04, i = ', i
       DO j = 1, 24
         READ (unitin, *, IOSTAT = ios)
       END DO
     END IF
!     IF (bp3==4 .and. bp2==3 .and. bp1==1 .and. bp4==1) &
!       & print *, '1-3-04, i = ', i
  END DO
  CLOSE(unitin)

  ! BP added this section to read in Tumbarumba daily LAI and interpolate to hourly
  OPEN (UNIT=unitin, FILE='tumb_lai.txt', FORM='formatted')
!  OPEN (UNIT=unitin, FILE='tumb02_std_input', FORM='formatted')
!  DO j=1, 110            ! this loop is for file tumb02_std_input
!     READ (unitin, *)
!  END DO
  DO j = 1, 365
     READ (unitin, *)  tempLAI(j)   ! read daily value from tumb_lai.txt
!     READ (unitin, *) i, tempLAI(j)   ! read daily value
  END DO
  DO j = 53, 365     ! because data starts from 22 Feb
    DO i = (j-53)*24+1, (j-52)*24    ! replicate for every hour in that day
      met%LAI(1,1,i) = tempLAI(j)
    END DO
  END DO
!  print *, 'i value after first year assignment: ', i
  DO j = 1, 365
    DO i = (j-1)*24+1, j*24
      DO bp1 = 0, 4
        met%LAI(1,1,7512+bp1*8760+i) = tempLAI(j)
      END DO
    END DO
  END DO
!  print *, 'i value after sixth year assignment: ', i
  DO j = 1, 203      ! because data ends at 22 Jul 23:00
    DO i = (j-1)*24+1, j*24
      met%LAI(1,1,7512+5*8760+i) = tempLAI(j)
    END DO
  END DO
!  print *, 'Final value = ', 7512+5*8760+i-1
  CLOSE(unitin)
  ! BP end of addition

  ! Any units corrections to data:
  met%Snowf=0
  met%Tair=met%Tair+273.15
  WHERE(met%SWdown<0.0)
     met%SWdown=0.0
  END WHERE
  DO o = 1, SIZE(met%Qair)
     CALL rh_sh(met%Qair(1,1,1,o), met%Tair(1,1,1,o), &
          met%PSurf(1,1,o),sph)
!     IF (sph < 0.0) THEN  ! BP to fix negative data for 2001Sep5 at 9am
!       PRINT *, o, sph, met%Qair(1,1,1,o), 'Qair', met%Qair(1,1,1,o-1)
!       CALL rh_sh(met%Qair(1,1,1,o+1), met%Tair(1,1,1,o+1), &
!          met%PSurf(1,1,o+1),sph)
!       PRINT *, o+1, sph, met%Qair(1,1,1,o+1), 'next Qair'
!       sph=met%Qair(1,1,1,o-1)
!     END IF
     met%Qair(1,1,1,o) = sph
  END DO
  met%Rainf = met%Rainf/REAL(timestepsize)

  ! Create NetCDF file:
  status = NF90_CREATE(fileout, NF90_CLOBBER, ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Put the file in define mode:
  status = NF90_REDEF(ncid)                   

  ! Define netcdf file dimensions:
  status = NF90_DEF_DIM(ncid, 'x', 1, xID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid, 'y', 1, yID) 
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid, 'z', za, zID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid, 'z_down', zs, zsID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_DIM(ncid, "t", NF90_UNLIMITED, tID)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Define dimension variables:
  status = NF90_DEF_VAR(ncid,'nav_lon',NF90_FLOAT,(/yID,xID/),nav_lonID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,nav_lonID,'long_name','Longitude')
  status = NF90_PUT_ATT(ncid,nav_lonID,'units','degrees_east')
  status = NF90_DEF_VAR(ncid,'nav_lat' ,NF90_float,(/yID,xID/),nav_latID)
  status = NF90_PUT_ATT(ncid,nav_latID, "long_name","Latitude")
  status = NF90_PUT_ATT(ncid,nav_latID, "units","degrees_north" )
  status = NF90_DEF_VAR(ncid,'x',NF90_INT,(/xID/),xvID)
  status = NF90_DEF_VAR(ncid,'y',NF90_INT,(/yID/),yvID)
  status = NF90_DEF_VAR(ncid,'elevation',NF90_FLOAT,(/zID/),levelID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_DEF_VAR(ncid,'time',NF90_DOUBLE,(/tID/),timeID)
  status = NF90_PUT_ATT(ncid,timeID,'long_name','Time')
  status = NF90_PUT_ATT(ncid,timeID,'units',timeunits)
  status = NF90_DEF_VAR(ncid,'timestp',NF90_FLOAT,(/tID/),timestpID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,timestpID,'long_name', 'Timesteps')
  status = NF90_PUT_ATT(ncid,timestpID,'units',timestpunits)

  ! Define met variables:
  status = NF90_DEF_VAR(ncid, "SWdown", NF90_FLOAT,  &
       (/ xID, yID, tID /), SWdownID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,SWdownID,"units",units%SWdown )
  status = NF90_PUT_ATT(ncid,SWdownID,"long_name",'Surface incident shortwave radiation' )
  status = NF90_DEF_VAR(ncid, "LWdown", NF90_FLOAT, &
       (/ xID, yID, tID /), LWdownID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,LWdownID,"units",units%LWdown)
  status = NF90_PUT_ATT(ncid,LWdownID, "long_name","Surface incident longwave radiation")
  status = NF90_PUT_ATT(ncid,LWdownID, "source","Cloudiness-adjusted Brutsaert synthesis")
  status = NF90_DEF_VAR(ncid, "Rainf", NF90_FLOAT, &
       (/ xID, yID, tID /), RainfID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,RainfID, "units",units%Rainf )
  status = NF90_PUT_ATT(ncid,RainfID, "long_name",'Rainfall rate')
  status = NF90_DEF_VAR(ncid, "Snowf", NF90_FLOAT, &
       (/ xID, yID, tID /), SnowfID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,SnowfID, "units",units%Snowf )
  status = NF90_PUT_ATT(ncid,SnowfID, "long_name",'Snowfall rate')
  status = NF90_DEF_VAR(ncid, "Tair", NF90_FLOAT, &
       (/ xID, yID, zID, tID /), TairID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,TairID, "units",units%Tair )
  status = NF90_PUT_ATT(ncid,TairID, "long_name",'Near surface air temperature')
  status = NF90_DEF_VAR(ncid, "PSurf", NF90_FLOAT,  &
       (/ xID, yID, tID /), PSurfID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,PSurfID, "units",units%PSurf )
  status = NF90_PUT_ATT(ncid,PSurfID, "long_name",'Surface pressure')
  status = NF90_DEF_VAR(ncid, "Wind", NF90_FLOAT, &
       (/ xID, yID, zID, tID /), WindID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,WindID,"units",units%Wind )
  status = NF90_PUT_ATT(ncid,WindID,"long_name",'Near surface module of the wind')
  status = NF90_DEF_VAR(ncid, "Qair", NF90_FLOAT, &
       (/ xID, yID, zID, tID /), QairID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,QairID, "units",units%Qair )
  status = NF90_PUT_ATT(ncid,QairID, "long_name",'Near surface specIFic humidity')
  status = NF90_DEF_VAR(ncid, "LAI", NF90_FLOAT, &
       (/ xID, yID, tID /), LaiID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,LaiID,"units",units%LAI)
  status = NF90_PUT_ATT(ncid,LaiID,"long_name","Leaf Area Index")
  status = NF90_PUT_ATT(ncid,LaiID,"source","provided by Ray Leuning")

  ! Define flux/state variable
  status = NF90_DEF_VAR(ncid, "Qle", NF90_FLOAT, &
       (/ xID, yID, tID /), QleID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,QleID, "units",units%Qle )
  status = NF90_PUT_ATT(ncid,QleID, "long_name",'Latent heat flux')
  status = NF90_DEF_VAR(ncid, "Qh", NF90_FLOAT,  &
       (/ xID, yID, tID /), QhID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,QhID, "units",units%Qh )
  status = NF90_PUT_ATT(ncid,QhID, "long_name",'Sensible heat flux')
  status = NF90_DEF_VAR(ncid, "NEE", NF90_FLOAT,  &
       (/ xID, yID, tID /), NEEID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,NEEID, "units",units%NEE )
  status = NF90_PUT_ATT(ncid,NEEID, "long_name",'Net Ecosystem Exchange')
  status = NF90_DEF_VAR(ncid, "Rnet", NF90_FLOAT, &
       (/ xID, yID, tID /), RnetID)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,RnetID, "units",units%Rnet )
  status = NF90_PUT_ATT(ncid,RnetID, "long_name",'Net absorbed radiation')

  ! Define parameters:
  ! Firstly soil and veg type:
  status=NF90_DEF_VAR(ncid,'isoil',NF90_FLOAT,(/xID,yID/),isoilID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,isoilID,"long_name",&
       "Soil type from Zobler")
  status = NF90_PUT_ATT(ncid,isoilID,"units","-")
  status=NF90_DEF_VAR(ncid,'iveg',NF90_FLOAT,(/xID,yID/),ivegID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,ivegID,"long_name",&
       "Vegetation type from Potter et al")
  status = NF90_PUT_ATT(ncid,ivegID,"units","-")
  ! albsoil :
  status=NF90_DEF_VAR(ncid,'albsoil',NF90_FLOAT,(/xID,yID/),albsoilID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,albsoilID,"long_name",&
       "Soil albedo")
  status = NF90_PUT_ATT(ncid,albsoilID,"units","-")
  ! clay (fraction of soil which is clay):
  status=NF90_DEF_VAR(ncid,'clay',NF90_FLOAT,(/xID,yID/),clayID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,clayID,"long_name",&
       "Fraction of soil which is clay")
  status = NF90_PUT_ATT(ncid,clayID,"units","-")
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
  ! sand (fraction of soil which is sand):
  status=NF90_DEF_VAR(ncid,'sand',NF90_FLOAT,(/xID,yID/),sandID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,sandID,"long_name",&
       "Fraction of soil which is sand")
  status = NF90_PUT_ATT(ncid,sandID,"units","-")
  ! silt (fraction of soil which is silt):
  status=NF90_DEF_VAR(ncid,'silt',NF90_FLOAT,(/xID,yID/),siltID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,siltID,"long_name",&
       "Fraction of soil which is silt")
  status = NF90_PUT_ATT(ncid,siltID,"units","-")
  ! za 
  status=NF90_DEF_VAR(ncid,'za',NF90_FLOAT,(/xID,yID/),zaID)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,zaID,"long_name",&
       "Reference height")
  status = NF90_PUT_ATT(ncid,zaID,"units","m")

  ! Write time of production to netcdf file:
  CALL DATE_AND_TIME(todaydate, nowtime)
  todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
  nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Production", &
       TRIM(todaydate)//' at '//TRIM(nowtime))
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Contact",contact)
  ! Other global details:
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Site_name",sitename)
  IF (status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_ATT(ncid,NF90_GLOBAL,"Length_of_dataset",datalength)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! End define mode:
  status = NF90_enddef(ncid) 
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Write variables:
  status = NF90_PUT_VAR(ncid, nav_lonID, longitude)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, nav_latID, latitude)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, xvID, xv)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, yvID, yv)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, levelID, elevation)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, timeID, time)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, timestpID, timestp)
  IF(status /= NF90_NOERR) CALL handle_err(status)

  status = NF90_PUT_VAR(ncid, SWdownID, met%SWdown)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, LWdownID, met%LWdown)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, RainfID, met%Rainf)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, SnowfID, met%Snowf)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, TairID, met%Tair)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, PSurfID, met%PSurf)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, WindID, met%Wind)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, QairID, met%Qair)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, LaiID, met%LAI)
  IF(status /= NF90_NOERR) CALL handle_err(status)

  status = NF90_PUT_VAR(ncid, QleID, flst%Qle)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, QhID, flst%Qh)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, NEEID, flst%NEE)
  IF(status /= NF90_NOERR) CALL handle_err(status)

  ! Non-ALMA:
  status = NF90_PUT_VAR(ncid, RnetID, flst%Rnet)
  IF(status /= NF90_NOERR) CALL handle_err(status)

  ! Parameters:
  status = NF90_PUT_VAR(ncid, isoilID, par%isoil)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, ivegID, par%iveg)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, albsoilID, par%albsoil)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, clayID, par%clay)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, frac4ID, par%frac4)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, hcID, par%hc)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, sandID, par%sand)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, siltID, par%silt)
  IF(status /= NF90_NOERR) CALL handle_err(status)
  status = NF90_PUT_VAR(ncid, zaID, par%za)
  IF(status /= NF90_NOERR) CALL handle_err(status)

  ! Close NetCDF file:
  status = NF90_close(ncid) 
  IF (status /= NF90_NOERR) CALL handle_err(status)

CONTAINS
  SUBROUTINE handle_err(status)
    INTEGER, INTENT(IN) :: status
    IF(status /= NF90_NOERR) THEN
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
    !  Internal Routines called:
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
  !------------------------------------------------------------------------------------

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


END PROGRAM texttonc
