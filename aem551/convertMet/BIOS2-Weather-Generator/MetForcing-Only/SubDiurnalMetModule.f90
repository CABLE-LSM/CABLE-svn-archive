
MODULE SubDiurnalMetModule
CONTAINS
SUBROUTINE DailyConstants (np,YearDay,LatDeg,WindDay,TempMinDay,TempMaxDay,TempMinDayNext, &
                           TempMaxDayPrev,WindDark,WindLite,SolarNorm,DecRad,LatRad, &
						   DayLength,TimeSunsetPrev,TimeSunrise,TimeMaxTemp,TimeSunset, &
						   TempSunsetPrev,TempSunset,TempNightRate,TempNightRatePrev, &
						   TempRangeDay,TempRangeAft)

!-------------------------------------------------------------------------------

! Routine for calculating solar and met parameters that are constant over the
! the course of a 24hr day. These have been split off from version 01 of 
! SubDiurnalMet which no longer calculates a day's worth of subdiurnal met in
! one call. Subdiurnal met for ntime times is now calculated using ntime calls,
! using the day-constants calculated here.

!-------------------------------------------------------------------------------
USE TypeDef


implicit none
! Parameters
real(sp),parameter:: Pi         = 3.14159265        ! Pi
real(sp),parameter:: PiBy2      = 1.57079632        ! Pi/2
real(sp),parameter:: SecDay     = 86400.0           ! Seconds/day
real(sp),parameter:: SolarConst = 1370*SecDay/1e6   ! Solar constant [MJ/m2/day]

! Global variables
integer(i4b),intent(in)           :: np       ! number of spatial elements
integer(i4b),intent(in)           :: YearDay  ! Day of the year (1-366)
real(sp),dimension(np),intent(in) :: LatDeg   ! Latitude [degs, neg for SH]

real(sp),dimension(np),intent(in) :: WindDay
real(sp),dimension(np),intent(in) :: TempMinDay      ! Daily minimum air temp [degC]
real(sp),dimension(np),intent(in) :: TempMaxDay      ! Daily maximum air temp [degC]                                                  Tx
real(sp),dimension(np),intent(in) :: TempMinDayNext  ! Daily minimum air temp tomorrow [degC]
real(sp),dimension(np),intent(in) :: TempMaxDayPrev  ! Daily maximum air temp yesterday [degC] 

real(sp),dimension(np),intent(out) :: WindDark    ! Wind [m/s] during night hrs
real(sp),dimension(np),intent(out) :: WindLite    ! Wind [m/s] during daylight hrs

real(sp),dimension(np),intent(out) :: SolarNorm   ! (daily solar)/(solar const): geometry
real(sp),              intent(out) :: DecRad      ! Declination in radians
real(sp),dimension(np),intent(out) :: LatRad      ! Latitude in radians
real(sp),dimension(np),intent(out) :: DayLength      
real(sp),dimension(np),intent(out) :: TimeSunsetPrev
real(sp),dimension(np),intent(out) :: TimeSunrise
real(sp),dimension(np),intent(out) :: TimeMaxTemp
real(sp),dimension(np),intent(out) :: TimeSunset
real(sp),dimension(np),intent(out) :: TempSunsetPrev
real(sp),dimension(np),intent(out) :: TempSunset
real(sp),dimension(np),intent(out) :: TempNightRate
real(sp),dimension(np),intent(out) :: TempNightRatePrev 
real(sp),dimension(np),intent(out) :: TempRangeDay
real(sp),dimension(np),intent(out) :: TempRangeAft

! Local variables
real(sp) :: YearRad
real(sp),dimension(np) :: LatDeg1 
real(sp),dimension(np) :: TanTan 
real(sp),dimension(np) :: HDLRad 
real(sp),dimension(np) :: TimeSunriseNext
real(sp) :: RatioWindLiteDark 

! -------------------------
! Downward solar irradiance
! -------------------------
LatDeg1 = sign(min(abs(LatDeg),89.9), LatDeg)   ! avoid singularity at pole
LatRad  = LatDeg1*Pi/180.0                      ! latitude in radians
YearRad = 2.0*Pi*(YearDay-1)/365.0              ! day of year in radians
!YearRadPrev = 2.0*Pi*(YearDay-2)/365.0          ! previous day of year in radians
                                                ! (Ok for YD - 2 = -1)

! DecRad = Declination in radians (+23.5 deg on 22 June, -23.5 deg on 22 Dec):
DecRad = 0.006918 - 0.399912*cos(YearRad) + 0.070257*sin(YearRad)     &
         - 0.006758*cos(2.0*YearRad) + 0.000907*sin(2.0*YearRad)      &
         - 0.002697*cos(3.0*YearRad) + 0.001480*sin(3.0*YearRad)
                                        ! Paltridge and Platt eq [3.7]

! Daylength: HDLRad = Half Day Length in radians (dawn:noon = noon:dusk):
TanTan = -tan(LatRad)*tan(DecRad)
where (TanTan .le. -1.0) 
  HDLRad = Pi                           ! polar summer: sun never sets
elsewhere (TanTan .ge. 1.0)
  HDLRad = 0.0                          ! polar winter: sun never rises
elsewhere
  HDLRad = acos(TanTan)                 ! Paltridge and Platt eq [3.21]
end where                                  ! (HDLRad = their capital PI)
DayLength = 24.0*2.0*HDLRad / (2.0*Pi)  ! Daylength (dawn:dusk) in hours

! Daily solar irradiance without atmosphere, normalised by solar constant,
! with both energy fluxes in MJ/m2/day, calculated from solar geometry:
SolarNorm =  &                          ! Paltridge and Platt eq [3.22]
  (HDLRad*sin(LatRad)*sin(DecRad) + cos(LatRad)*cos(DecRad)*sin(HDLRad)) / Pi

! ----
! Wind
! ----

RatioWindLiteDark = 3.0                 ! (daytime wind) / (nighttime wind)
WindDark  = WindDay / ( (24.0-DayLength)/24.0 + RatioWindLiteDark*DayLength/24.0 )
WindLite  = WindDay / ( (1.0/RatioWindLiteDark)*(24.0-DayLength)/24.0 + DayLength/24.0 )

! -----------
! Temperature
! -----------
! These are parameters required for the calculation of temperature according to 
! Cesaraccio et al 2001 for sunrise-to-sunrise-next-day. Because we are calculating temp for
! midnight-to-midnight, we need to calculate midnight-to-sunrise temps using data for the 
! previous day (-24h), hence the extra parameters denoted by *'s below which are not
! mentioned per se in Cesaraccio. Cesaraccio symbology for incoming met data is included
! here as comments for completeness:
!                                                                 Sym in Cesaraccio et al 2001
! TempMinDay                                                        Tn
! TempMaxDay                                                        Tx
! TempMinDayNext                                                    Tp
! TempMaxDayPrev                                                    * Tx-24h
 
TimeSunrise = (acos(tan(LatRad)*tan(DecRad)))*12./Pi              ! Hn
TimeSunset  = TimeSunrise + DayLength                             ! Ho
TimeSunsetPrev = TimeSunset - 24.                                 ! * Ho-24h (a negative hour)  
TimeMaxTemp = TimeSunset - 4.                                     ! Hx
TimeSunriseNext = TimeSunrise + 24.                               ! Hp
TempSunset = TempMaxDay - & 
             (0.39 * (TempMaxDay - TempMinDayNext))               ! To
TempSunsetPrev = TempMaxDayPrev - &
                 (0.39 * (TempMaxDayPrev - TempMinDay))           ! * To-24h
TempRangeDay = TempMaxDay - TempMinDay                            ! alpha = Tx-Tn
TempRangeAft = TempMaxDay - TempSunset                            ! R = Tx-To
TempNightRate = (TempMinDayNext - TempSunset)/ &
                sqrt(TimeSunriseNext-TimeSunset)                  ! b = (Tp-To)/sqrt(Hp-Ho)
TempNightRatePrev = (TempMinDay - TempSunsetPrev)/ &
                sqrt(TimeSunrise-TimeSunsetPrev)                  ! * b-24h = (Tn-(To-24h))/sqrt(Hn-(Ho-24h))

END SUBROUTINE DailyConstants

!*******************************************************************************

SUBROUTINE SubDiurnalMet (np,itime,ntime,delT,SolarMJday,TempMinDay,PrecipDay,PrecipDayNext,&
                          VapPmbDay,VapPmb0900,VapPmb1500, &
                          VapPmb1500prev, VapPmb0900Next,PmbDay,WindDark,WindLite, &
                          SolarNorm,DecRad,LatRad,TimeSunsetPrev,TimeSunrise,TimeMaxTemp, &
						  TimeSunset,TempSunsetPrev,TempSunset,TempNightRate,TempNightRatePrev, &
						  TempRangeDay,TempRangeAft,PhiSd,PhiLd,Precip,Wind,Temp,VapPmb,Pmb,coszen)

!-------------------------------------------------------------------------------
! Routine for downscaling daily met data to ntime subdiurnal values per day,
! at the instants ((it/ntime)*2400, it=1,ntime).
!
! ALGORITHMS:
! * Downward solar irradiance:
!   * Fit a daily irradiance time series based on solar geometry, normalised to
!     observed daily total SolarMJDay. Reference: Paltridge and Platt (1976).
! * Downward longwave irradiance:
!   * Computed from downscaled air temperature using Swinbank (1963) formula.
! * Precipitation:
!   * Assume steady through 24 hours at average value (needs improvement?)
! * Wind:
!   * Set RatioWindLiteDark = (daytime wind) / (nighttime wind), a predetermined 
!     value, typically 3. Then calculate wind by day and wind by night (both 
!     steady, but different) to make the average work out right.
! * Temperature:
!   * Fit a sine wave from sunrise (TempMinDay) to peak at TempMaxDay, another
!     from TempMaxDay to sunset, TempMinDay, and a square-root function from 
!     sunset to sunrise the next day (Cesaraccio et al 2001).
! * Water Vapour Pressure, Air Pressure:
!   * Hold constant. Return rel humidity at TempMinDay as a diagnostic on obs. 
!
! HISTORY:
! * 04-sep-2003 (MRR): 01 Written and tested in Program SubDiurnalWeather
! * 30-nov-2007 (PRB): 02 Vectorise with deferred shape arrays, adding np dimension.
! * 11-dec-2007 (PRB): 03 Switch to Cesaraccio et al 2001 algorithm for temperature
! * 28/02/2012 (VH) :  04 cahnge precip from uniform distribution to evenly distributed over
! the periods 0600:0700; 0700:0800; 1500:1600; 1800:1900
!-------------------------------------------------------------------------------
USE TypeDef
USE DateFunctions
USE constants
implicit none
! Parameters
!real(sp),parameter:: Pi         = 3.14159265        ! Pi
real(sp),parameter:: PiBy2      = 1.57079632        ! Pi/2
!real(sp),parameter:: SecDay     = 86400.0           ! Seconds/day
real(sp),parameter:: SolarConst = 1370*SecDay/1e6   ! Solar constant [MJ/m2/day]
real(sp),parameter:: epsilon  = 0.736
! Global variables

!   Space & Time
integer(i4b),intent(in) :: np       ! number of spatial elements
integer(i4b),intent(in) :: itime    ! current timestep
integer(i4b),intent(in) :: ntime    ! number of subdiurnal steps in 24 hr
real(sp),intent(in) :: delT ! subdiurnal timestep (in s)
!   Daily Met
real(sp),dimension(np),intent(in) :: SolarMJday  ! 24hr-total incident solar radiation [MJ/day]
real(sp),dimension(np),intent(in) :: TempMinDay  ! Minimum temperature current day (°C)
real(sp),dimension(np),intent(in) :: PrecipDay   ! 24hr-total precipitation [m/day] (current day = 24h from 0900 on previous day)
real(sp),dimension(np),intent(in) :: PrecipDayNext   ! 24hr-total precipitation [m/day] (next day = 24 h from 0900 on current day)
real(sp),dimension(np),intent(in) :: VapPmbDay    ! 24hr-av water vapour pressure [mb]
real(sp),dimension(np),intent(in) :: VapPmb0900    ! 0900 water vapour pressure [mb]
real(sp),dimension(np),intent(in) :: VapPmb1500    ! 1500 water vapour pressure [mb]
real(sp),dimension(np),intent(in) :: VapPmb1500Prev    ! 1500 (prev day) water vapour pressure [mb]
real(sp),dimension(np),intent(in) :: VapPmb0900Next    ! 0900(next day) water vapour pressure [mb]
real(sp),dimension(np),intent(in) :: PmbDay       ! 24hr-av pressure [mb]
real(sp),dimension(np),intent(in) :: WindDark    ! Wind [m/s] during night hrs
real(sp),dimension(np),intent(in) :: WindLite    ! Wind [m/s] during daylight hrs
!   Solar and Temperature Params
real(sp),dimension(np),intent(in) :: SolarNorm   ! (daily solar)/(solar const): geometry
real(sp),              intent(in) :: DecRad      ! Declination in radians
real(sp),dimension(np),intent(in) :: LatRad      ! Latitude in radians
real(sp),dimension(np),intent(in) :: TimeSunsetPrev
real(sp),dimension(np),intent(in) :: TimeSunrise
real(sp),dimension(np),intent(in) :: TimeMaxTemp
real(sp),dimension(np),intent(in) :: TimeSunset
real(sp),dimension(np),intent(in) :: TempSunsetPrev
real(sp),dimension(np),intent(in) :: TempSunset
real(sp),dimension(np),intent(in) :: TempNightRate
real(sp),dimension(np),intent(in) :: TempNightRatePrev 
real(sp),dimension(np),intent(in) :: TempRangeDay
real(sp),dimension(np),intent(in) :: TempRangeAft
!   Hourly Met outgoing
real(sp),dimension(np),intent(out):: PhiSd     ! downward solar irradiance [W/m2]
real(sp),dimension(np),intent(out):: PhiLd     ! down longwave irradiance  [W/m2]
real(sp),dimension(np),intent(out):: Precip    ! precip [mm/h]
real(sp),dimension(np),intent(out):: Wind      ! wind   [m/s]
real(sp),dimension(np),intent(out):: Temp      ! temp   [degC]
real(sp),dimension(np),intent(out):: VapPmb    ! vapour pressure [mb]
real(sp),dimension(np),intent(out):: Pmb       ! pressure [mb]
real(sp),dimension(np),intent(out):: coszen      ! cos(theta)
! Local variables
real(sp) :: TimeNoon , test1, test2, adjust_fac(np)
real(sp) :: TimeRad
real(sp) :: rntime    ! Real version of ntime 
real(sp) :: ritime    ! Real version of current time
real(sp),dimension(np):: PhiLd_Swinbank     ! down longwave irradiance  [W/m2]


!-------------------------------------------------------------------------------

ritime = real(itime)*delT/3600.  ! Convert the current time to real
rntime = real(ntime)*delT/3600.  ! Convert ntime to real

! Instantaneous downward hemispheric solar irradiance PhiSd
TimeNoon = ritime/rntime - 0.5   ! Time in day frac (-0.5 to 0.5, zero at noon)
TimeRad  = 2.0*Pi*TimeNoon       ! Time in day frac (-Pi to Pi, zero at noon)
where (ritime >= TimeSunrise .and. ritime <= TimeSunset) ! Sun is up
  PhiSd = max((SolarMJDay/SolarNorm)  &   ! PhiSd [MJ/m2/day]
      * ( sin(DecRad)*sin(LatRad) + cos(DecRad)*cos(LatRad)*cos(TimeRad) ), 0.0)  ! fix to avoid negative PhiSD vh 13/05/08
                                        ! Paltridge and Platt eq [3.4]
  coszen = ( sin(DecRad)*sin(LatRad) + cos(DecRad)*cos(LatRad)*cos(TimeRad) )
elsewhere ! sun is down
  PhiSd = 0.0
  coszen = 0.0
end where

PhiSd    = PhiSd*1e6/SecDay       ! Convert PhiSd: [MJ/m2/day] to [W/m2]
!test1 = minval(PhiSD)
!test2 = maxval(PhiSD)

!write(99,'(6000e14.5)')  test1, test2

! -------------
! Precipitation
! -------------

!Precip = PrecipDay*1000./rntime  ! convert from m/d to mm/h
!Precip = (PrecipDay*1000.*9./24. + PrecipDayNext*1000.*15./24.)/rntime
!Precip = (PrecipDay*1000.*9./24. + PrecipDayNext*1000.*15./24.)/24.  ! hourly precip [mm/h]


if ((ritime >= 15. .and. ritime < 16.).or.(ritime >= 18. .and. ritime < 19.)) then
	Precip = PrecipDay*1000./2.
endif

! ----
! Wind
! ----

where (ritime >= TimeSunrise .and. ritime <= TimeSunset) ! Sun is up
  Wind = WindLite
elsewhere ! Sun is down
  Wind = WindDark
end where

! -----------
! Temperature
! -----------
! Calculate temperature according to Cesaraccio et al 2001, including midnight to 
! sunrise period using previous days info, and ignoring the period from after the 
! following midnight to sunrise the next day, normally calculated by Cesaraccio.

where (ritime > 0. .and. ritime <= TimeSunrise)
! Midnight to sunrise
  Temp = TempSunsetPrev + TempNightRatePrev * sqrt(ritime - TimeSunsetPrev)
elsewhere (ritime > TimeSunrise .and. ritime <= TimeMaxTemp) 
! Sunrise to time of maximum temperature
  Temp = TempMinDay + &
         TempRangeDay *sin (((ritime-TimeSunrise)/(TimeMaxTemp-TimeSunrise))*PiBy2)
elsewhere (ritime > TimeMaxTemp .and. ritime <= TimeSunset) 
! Time of maximum temperature to sunset
  Temp = TempSunset + &
         TempRangeAft * sin(PiBy2 + ((ritime-TimeMaxTemp)/4.*PiBy2))
elsewhere (ritime > TimeSunset .and. ritime <= 24.)
! Sunset to midnight
  Temp = TempSunset + TempNightRate * sqrt(ritime - TimeSunset)
end where
      
! -----------------------------------
! Water Vapour Pressure, Air Pressure
! -----------------------------------

VapPmb = VapPmbDay
Pmb    = PmbDay

if (ritime <= 9.) then
   ! before 9am
   VapPmb = VapPmb1500Prev + (VapPmb0900 - VapPmb1500Prev) * (9. + ritime)/18.
elseif (ritime > 9 .and. ritime <= 15.) then
! between 9am and 15:00
   VapPmb = VapPmb0900 + (VapPmb1500 - VapPmb0900) * (ritime - 9.)/(15.-9.)
elseif (ritime > 15.) then
! after 15:00
  VapPmb = VapPmb1500 + (VapPmb0900Next - VapPmb1500) * (ritime - 15.)/18.
end if

! ----------------------------
! Downward longwave irradiance
! ----------------------------

PhiLd_Swinbank = 335.97 * (((Temp + 273.16) / 293.0)**6)   ! [W/m2] (Swinbank 1963)

! -------------------------------
 ! Alternate longwave formulation
PhiLd = epsilon * SBoltz * (Temp + 273.16)**4       ! [W/m2] (Brutsaert)
where (PhiSd.gt.50.0)
	adjust_fac = ((1.17)**(SolarNorm))/1.17
elsewhere
	adjust_fac = 0.9
endwhere

PhiLd = PhiLd /adjust_fac * (1.0 + PhiSd/8000.)       ! adjustment (formulation from Gab Abramowitz)

where ((PhiLd.gt.500.00).or.(PhiLd.lt.100.00))
	PhiLd = PhiLd_Swinbank
endwhere

if (any((PhiLd.gt.750.00).or.(PhiLd.lt.100.00))) then
!write(*,*) 'PhiLD out of range'
endif


END SUBROUTINE SubDiurnalMet

END MODULE SubDiurnalMetModule
