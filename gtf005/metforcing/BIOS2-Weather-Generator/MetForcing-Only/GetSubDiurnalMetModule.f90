MODULE GetSubDiurnalMetModule
CONTAINS 
Subroutine GetSubdiurnalMet(Ttime,MM,MMprev,MMnext,VV,hMM)

USE TypeDef
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
USE SubDiurnalMetModule
USE GlobalDefs
implicit none

! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time from host program
real(sp),intent(in)::  MM(:,:)      ! MM(np,nmm) = daily met variables
real(sp),intent(in)::  MMprev(:,:)      ! MM(np,nmm) = daily met variables
real(sp),intent(in)::  MMnext(:,:)      ! MM(np,nmm) = daily met variables
real(sp),intent(in)::  VV(:,:)      ! VV(np,nmm) = spatially variable params
real(sp),intent(inout):: hMM(:,:,:)     ! hMM(np,nhmm,ntime) = hourly met variables
! Local variables
integer(i4b) :: itry, it, imm, ic, iunit, k
integer(i4b) :: MetDay, MetMonth, MetYear, doy
integer(i4b) :: ObsDay, ObsMonth, ObsYear
type(dmydate):: TTDate              ! current date from host program
type(dmydate):: MetDate
real(sp):: rntime 
real(sp),     parameter:: SolarConst = 1370*SecDay/1e6   ! Solar constant [MJ/m2/day]
real(sp), dimension(np) :: SolarMJday   ! 24hr-total incident solar radiation [MJ/day]
real(sp)                :: DecRad       ! Declination in radians
real(sp), dimension(np) :: LatRad       ! Latitude in radians
real(sp), dimension(np) :: PrecipDay    ! 24hr-total precipitation [m/day] 
real(sp), dimension(np) :: PrecipDayNext    ! 24hr-total precipitation [m/day] 
real(sp), dimension(np) :: WindDay      ! 24hr-av wind [m/s] 
real(sp), dimension(np) :: TempMinDay   ! Daily minimum air temp [degC]
real(sp), dimension(np) :: TempMaxDay   ! Daily maximum air temp [degC]
real(sp), dimension(np) :: TempMinDayNext
real(sp), dimension(np) :: TempMaxDayPrev
real(sp), dimension(np) :: WindDark     ! Wind [m/s] during night hrs
real(sp), dimension(np) :: WindLite     ! Wind [m/s] during daylight hrs
real(sp), dimension(np) :: VapPmbDay    ! 24hr-av water vapour pressure [mb]
real(sp), dimension(np) :: VapPmb0900    ! 0900 water vapour pressure [mb]
real(sp), dimension(np) :: VapPmb1500    ! 1500 water vapour pressure [mb]
real(sp), dimension(np) :: VapPmb1500Prev    ! 1500 (prev day) water vapour pressure [mb]
real(sp), dimension(np) :: VapPmb0900Next    ! 0900(next day) water vapour pressure [mb]
real(sp), dimension(np) :: PmbDay       ! 24hr-av pressure [mb]
real(sp), dimension(np) :: DayLength       ! day length (dawn:dusk)    [hr]
real(sp), dimension(np) :: SolarNorm       ! (daily solar irradiance)/(solar const)
real(sp), dimension(np) :: SolarFracObs    ! (obs daily solar)/(solar geometry value)
real(sp), dimension(np) :: SolarTotChk     ! check: day-integrated solar [MJ/m2/day] 
real(sp), dimension(np) :: WindAvgChk      ! day-average wind (check)
real(sp), dimension(np) :: RelHumTempMin   ! rel humidity at minimum air temp (check)
real(sp), dimension(np) :: PhiSd     ! downward solar irradiance [W/m2]
real(sp), dimension(np) :: coszen     ! cos(theta)
real(sp), dimension(np) :: PhiLd     ! downward longwave irradiance [W/m2]

real(sp), dimension(np) :: Wind      ! wind   [m/s]
real(sp), dimension(np) :: Temp      ! temp   [degC]
real(sp), dimension(np) :: VapPmb    ! vapour pressure [mb]
real(sp), dimension(np) :: Pmb       ! pressure [mb]

real(sp), dimension(np) :: TimeSunsetPrev
real(sp), dimension(np) :: TimeSunrise
real(sp), dimension(np) :: TimeMaxTemp
real(sp), dimension(np) :: TimeSunset
real(sp), dimension(np) :: TempSunsetPrev
real(sp), dimension(np) :: TempSunset
real(sp), dimension(np) :: TempNightRate
real(sp), dimension(np) :: TempNightRatePrev 
real(sp), dimension(np) :: TempRangeDay
real(sp), dimension(np) :: TempRangeAft

real(sp), dimension(np,ntime) :: PhiSdAll     ! downward solar irradiance [W/m2]
real(sp), dimension(np,ntime) :: PhiLdAll     ! downward longwave irradiance [W/m2]
real(sp), dimension(np,ntime) :: WindAll      ! wind   [m/s]
real(sp), dimension(np,ntime) :: TempAll      ! temp   [degC]
real(sp), dimension(np,ntime) :: rh      ! relative humidity  
integer(i4b) :: ip, itime
                                    
!-------------------------------------------------------------------------------

! * set output arrays to null obs flag (values will be overwritten)
hMM   =  -999.0 

! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime, TargetMM=MM,TargethMM=hMM,TargetVV=VV)
! * Construct date from TTime, in date format and as char string for filenames
TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
doy = YearDay(TTDate)
rntime = real(ntime)

!-------------------------------------------------------------------------------

SolarMJDay     = SolarMJ
PrecipDay      = Precip
PrecipDayNext  = MMNext(:,2)
WindDay        = 1.0
VapPmbDay      = vph09
VapPmb0900     = vph09
VapPmb1500     = vph15
PmbDay         = 1000.0
TempMinDay     = TempMin
TempMinDayNext = MMnext(:,4)
TempMaxDay     = TempMax
TempMaxDayPrev = MMprev(:,3)
VapPmb1500Prev = MMprev(:,6)
VapPmb0900Next = MMNext(:,5)



CALL DailyConstants (np,doy,LatDeg,WindDay,TempMinDay,TempMaxDay,TempMinDayNext, &
                     TempMaxDayPrev,WindDark,WindLite,SolarNorm,DecRad,LatRad, &
					 DayLength,TimeSunsetPrev,TimeSunrise,TimeMaxTemp,TimeSunset, &
					 TempSunsetPrev,TempSunset,TempNightRate,TempNightRatePrev, &
					 TempRangeDay,TempRangeAft)

SolarTotChk = 0.0
WindAvgChk  = 0.0

do itime = 1,ntime

  CALL SubDiurnalMet (np,itime,ntime,delT,SolarMJday,TempMinDay,PrecipDay,PrecipDayNext, &
                      VapPmbDay,VapPmb0900,VapPmb1500, &
                      VapPmb1500prev, VapPmb0900Next,PmbDay,WindDark,WindLite, &
                      SolarNorm,DecRad,LatRad,TimeSunsetPrev,TimeSunrise,TimeMaxTemp, &
					  TimeSunset,TempSunsetPrev,TempSunset,TempNightRate,TempNightRatePrev, &
					  TempRangeDay,TempRangeAft,PhiSd,PhiLd,hPrecip(:,itime),Wind,Temp,VapPmb,Pmb,coszen)

  SolarTotChk = SolarTotChk + PhiSd/1e6*SecDay  ! check: day-integrated solar [cnvrt to MJ/m2/day] 
  WindAvgChk  = WindAvgChk + Wind               ! check: day-integrated wind

! Collect hourly values into arrays for convenient output
  PhiSdAll(:,itime) = PhiSd  ! downward solar irradiance [W/m2]
  PhiLdAll(:,itime) = PhiLd  ! downward longwave irradiance [W/m2]
  WindAll(:,itime)  = Wind   ! wind   [m/s]
  TempAll(:,itime)  = Temp   ! temp   [degC]
  hqv(:,itime) = VapPmb/Pmb*rmw/rma ! specific humidity
  hcoszen(:,itime) = coszen
  hpmb(:,itime) = pmb
  rh(:,itime) = VapPmb/ESatf(Temp)
  !write(7,*) itime, VapPmb(1)/ESatf(Temp(1))
 
end do

! variables for hMM
  hFsd = PhiSdAll; hFld = PhiLdAll; hUa = WindAll; hTc = TempAll


SolarTotChk = SolarTotChk / rntime
RelHumTempMin = VapPmbDay / ESatf(TempMinDay)
WindAvgChk = WindAvgChk / rntime
! Observed daily solar irradiance as frac of value from solar geometry (should be < 1)
SolarFracObs = SolarMJDay / (SolarNorm*SolarConst + 1.0e-6)

END subroutine GetSubdiurnalMet

END MODULE GetSubDiurnalMetModule
