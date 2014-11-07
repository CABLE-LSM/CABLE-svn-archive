! Program BIOS2
! Vanessa Haverd, Peter Briggs, Mike Raupach Aug 2012
!
! PROGRAM STRUCTURE
! -----------------
! MODULE TypeDef
!   * This module explicitly defines the sizes of variable types
! MODULE Constants
!   * Physical constants and parameters
!   USE TypeDef
! MODULE DateFunctions
!   * Date arithmetic with derived type dmydate
!   USE TypeDef
!   CONTAINS AddDay, SubDay, EQDates, NEDates, LTDates, GTDates, LEDates, GEDates,
!            JulianDay, GregDate, LeapDay, DaysInMonth, YearDay, DayDifference,
!            LegalDate, EndMonth, EndYear
! MODULE Utils
!   * Utility subroutines from library Milib
!   USE TypeDef
!   CONTAINS Qsatf, Esatf, Epsif, AstronDailySolar, AstronSZAtoTime, 
!            ComSkp, AnStep, AnMin
! MODULE PointerModule
!   * Defines and assigns pointers to specific variables in WaterDynStep
!   USE TypeDef
!   CONTAINS PointAll
! MODULE WaterDynModule
!   * This module contains all subroutines for WaterDyn dynamic and observation
!     models at a single time step.
!   * Everything in WaterDynModule is private unless otherwise declared
!   USE TypeDef, Utils
!   CONTAINS WaterDynStep (CONTAINS InitModel, InitVV), InitStores,
!            ReadControlFile, FindUniqueRegions, ReadMetObs, 
!            WaterDynFlux, WaterDynObs,
!            SpatialStats, TimeAverage, 
!            TimeSeriesOutput, CatchmentTSOutput, MapOutput
! PROGRAM WaterDynM
!   USE TypeDef, WaterDynModule
!
! GENERIC ARRAYS (declared in WaterDynStep and assigned with PointAll)
! --------------
!   TTime(ntt)   = time variables
!   XX(np,nxx)   = state variables
!   FF(np,nff)   = fluxes and dXXdt
!   DD(np,ndd)   = derived quantities
!   MM(np,nmm)   = met variables
!	hMM(np,nhMM,ntime) = hourly met variables
!   UU(nuu)      = spatially uniform params
!   VV(np,nvv)   = spatially variable params
!   ZZP(np,nzzP) = predicted point observables 
!   AAP(np,naaP) = actual and ancillary point observations 
!   ZZPh(np,nzzPh,ntime) = predicted point observables (hourly)
!   AAPh(np,naaPh,ntime) = actual and ancillary point observations (hourly)
!   ZZC(nc,nzzC) = predicted catchment observables
!   AAC(nc,naaC) = actual and ancillary catchment observations
!                  (ancillaries are variables needed by obs model, like obs time)
!
! FILE LOCATIONS
! --------------
!       program file            HomeDir\'WaterDyn'xxx.f90
!       control file            HomeDir\CtlFileName.ctl (copy in DataOutRootDir\DomName)
!       spatial params (VV)     DataInRootDir\DomName\'Params'\VVFileName
!       met data       (MM)     DataInRootDir\DomName\'Met'\MetFileName(1..4)
!		Remote senisng driver data (RR) DataInRootDir\DomName\'Met'\RRFileName(1..5)
!       NDVI    obs data (AA)   DataInRootDir\DomName\'Obs'\FileName
!       LST     obs data (AA)   DataInRootDir\DomName\'Obs'\FileName
!       LSTTime obs data (AA)   DataInRootDir\DomName\'Obs'\FileName
!       discharge data   (AA)   DataInRootDir\'DischargeData.yyyymm'\iCatchID'monrun'.bin
!       catchment properties    DataInRootDir\'DischargeData.yyyymm'\FileName
!       output (time series)    DataOutRootDir\DomName\Runyyy\OutputTS\'OutAv'x.csv (x=0,1,2,3)
!                               DataOutRootDir\DomName\Runyyy\OutputTS\'CatchAv'x.csv (x=0,1,2,3)
!       output (maps)           DataOutRootDir\DomName\Runyyy\OutputMap\QQname.xxx.Date.bin
!                                              (xxx = run,ann,mth,day; yyy = RunID)
!
! HISTORY
! -------
! WaterDyn01 (17-mar-05) to WaterDyn11M (16-jul-2006): see history in older versions
!
! WaterDyn12M (MRR, 26-feb to 03-mar-2007)
! * Go to 2-layer soil hydrology.
!
! WaterDyn13M (MRR, 07-mar-07)
! * Include daily runoff data (ADisCD) in ReadMetObs
! * Remove SpatialModeFlag: henceforth only use "combined-met" (SpatialModeFlag=1).
!
! WaterDyn14M (MRR, 13-mar-07)
! * ALLOCATE all arrays with np or nc dimension, to make code general for any domain
!   without need to recompile.  Arrays are allocated and initialised in 2 places: 
!   (1) InitMain: XX0, XX1, ZZP1, AAP1, ZZC1, AAC1
!   (2) InitModel: FF, MM, DD, VV, accumulator arrays with np dimension
! * (14-mar-07) 
!   * Model changes required for sensitivity tests (CMT program WaterDyn14S):
!     * Move declaration of XXSpTAv, XXSpTAvCount (and similar for FF, MM etc) from 
!       WaterDynStep to WaterDynModule.  These arrays can now be seen in host program.
!     * Improve accumulation and include AvType = (0,1,2,3) = (day,mth,ann,run) for 
!       each of domain-av TS, catchment-av TS, maps. Reinitialise all accumulator 
!       arrays after each output (for sensitivity tests).
!
! WaterDyn15M (MRR, 19-mar-07)
! * Introduce new met and AVHRR binary data files, in form (idate(3), vals(np))
!   where idate = (D,M,Y) as 3 4-byte integers. 
!
! WaterDyn16M (MRR, 19-mar-07)
! * Combine forward run (WaterDyn15M) and sensitivity run (WaterDyn15S, written by CMT)
!   into a single code, using SensTestFlag = (0,1,2,3) to switch between (M,S) modes
!   * If SensTestFlag=0, do forward run 
!   * If SensTestFlag=1, do sensitivity test on nominated parameters.
!     Note: diags, TS, map, PEST output are written only for reference case.
! * Include separate multipliers ZSoil1Mult, ZSoil1Mult for 2 soil layers
!
! WaterDyn17M (MRR, 20-apr-07)
! * Read LST and LSTTime from files specified in CTL file
! * To look at TS output for a single point, include flag ipDiag:
!     ipDiag = (0, >=1) for TS output = (domain-av, point ipDiag)
! * Fix problem in AllocLea code in WaterDynFlux:
!     OLD (wrong):   WRelAv   = (WRel1*FWTr1W + WRel1*FWTr2W) / FWTraW
!     NEW (right):   WRelAv   = (WRel1*FWTr1W + WRel2*FWTr2W) / FWTraW
! * Fix problem in MapOutput:
!     OLD (wrong):   OPEN(unit=170, file=GenFileName, form='unformatted')
!     NEW (right):   OPEN(unit=170, file=GenFileName, form='binary')
!     This problem put an extra 4-byte char at the start and end of all files 
!     written by MapOutput up to this point!!
! * (Mon 23-apr-2007) 
!   Fix errors in LST obs model but leave overpass time = noon in V17
!  
! WaterDyn18M (MRR, sun 29-apr-07)
! * In LST obs model, properly set overpass time and rerun for Murrumbidgee
! * (Wed 02-may-07)
!   * For EnKF compatibility, move InitMain (initialisation of XX0) from a CONTAINed
!     subroutine in host program to a subroutine in WaterDynModule, like WaterDynStep.
!     Rename InitMain to InitStores (its new sole function).
!   * Other small changes for EnKF, PEST compatibility:
!     * Don't open TS output files unless needed
!     * Include UUmin, UUmax, UUPESTFlag, UUKFFlag  in ctl file
! * (Wed 30-may-07)
!   * Split discharge data into 2 directories, 
!     DischargeData.200509 and DischargeData.200610.
!     Monthly and daily discharge data files in these directories need different READs.
!     Hardwire DischargeData.200509 read formats in V18.
! * (Mon 04-jun-07)
!   * In subroutine WaterDynObs, ensure that ZLST is always calculated whether
!     or not ALSTTime is available.  If no ALSTTime from data, use ALSTTime = 12.00.
!   * Redefine ZDisPM and ZDisCM, to:
!     *   ZDisPM = ZDisPD * DaysInMonth   for all days
!         ZDisCM = ZDisCD * DaysInMonth   for all days
!     * Thus ZDisPM, ZDisCM are both in m/mth (since ZDisPD, ZDisCD are in m/day).
!     * This means that a monthly AVERAGE of ZDisCM is directly comparable with
!       ADisCM (which is only defined on last day of month, and is ErrVal otherwise).
!     * Previously, ZDisPM, ZDisCM were set to SUM of ZDisPD, ZDisCD on last day
!       of month, and to ErrVal otherwise.
!
! WaterDyn19  (CMT, xx-jun-2007)
!
! WaterDyn20M (MRR, thurs 28-jun-07)
! * Include new discharge obs model:
!     FWRunC, FWLchC = sum(FWRun), sum(FWLchB) over catchment
!     ZRunCD(ic) = low-pass filtered FWRunC with time scale TZRunCat(ic)
!     ZLchCD(ic) = low-pass filtered FWLchC with time scale TZLchCat(ic)
!     ZDisCD(ic) = ZRunCD(ic) + ZLchCD(ic) = daily catchment outflow [m/day]
! * Shorten ZZP, AAP by removing ZDisPD, ZDisPM, so ZZP = (ZNDVI, ZLST, ZLSTTime)
! * Extend  ZZC, AAC to ZZC = (ZDisCD, ZDisCM, ZRunCD, ZLchCD)
!   (ZRunCD, ZLchCD included to provide accumulation access for main program for EnKF)
! * (Fri 29-jun-07)
!   * Include changes in CMT WaterDyn19P:
!     * include AAPdaymth, AACdaymth flags
!     * different DischargeData formats for 200509, 200610
!     * different ZDisCM models for PEST, EnKF (select with nEnsemble)
!     * PestOutput subroutine from WD19P
!
! WaterDyn21M (MRR, Tue 22-jan-08 to Thu 24-jan-08)
! * Include option for determining FracV: 
!   FracVegFlag = (0,1) => FracV from (external VegFPC, internal CLea)
! * Fix problem with reading of DischargeData in sensitivity mode. 
!   For all input data files:
!   * ensure REWIND is done immediately after OPEN
!   * apply action='read' in OPEN statement
! * Use 200610 discharge data throughout
! * Introduce maximum rLAI into UU parameters
! * Introduce TZRunDef, TZLchDef into UU parameters: these are default values of
!   TZRunCat(nc), TZLchCat(nc), used for param estimation on a single catchment.
! * (29-jan-08) Initialise ZRunCD=0, ZLchCD=0 in InitModel 
!
!CableDyn01M (VH Monday 11-feb-08 to Friday 15-Feb-08)
! * Created new variable hMM to hold hourly met data required for CABLE
! * Intoduced Mike's SubDiurnalMet to generate hourly met forcing
! * Introduced addtional parameters required for CABLE to UU and VV
! * ForwardModelFlag allows user to choose between Waterdyn and CABLE
! * Created subroutine "InitCable" to transfer information from Waerdyn variables to CABLE structures
! * Within WaterDynStep, now have the option to call CABLE 24 times and then aggregate results
! * Checks for water balance in canopy and soil after daily aggregation
! * Changed nmaes of Water fluxes so that "A" refers to teh A horizon, "B" refers to B horizon and 1:6 refere to 6 soil layers
! * Extract info required for LST prediction after call to CABLE
! * Elimiate call to WaterDynFlux within WaterDyn Obs

!CableDyn01M (VH Thursday 21-feb-08 to -Feb-08)
!*  define and allocate MMprev and MMnext
! * Introduce Peter's weather generator

!CableDyn02M (VH 15-May-08)
! Changed scope of VV so that it is a global variable
!modified sensitivity analysis to enable sensitivity tests of VV

!CableDyn03M (VH 20-May_08)
! enable monthly sensitivy analysis
! create integer nMonths (number of times end-of-month is reached)

!CableDyn03a (VH 23_May-08)
! Modified Subdiurnalmet to linearly interpolate between 0900 and 1500 vapour pressure data

!CableDyn04 (VH 29_May-08)
! Modified Subdiurnalmet to use next day's precip(*15/24) and current day's precip(*9/24)

!CableDyn04 (VH 4-Jun-08)
!Corrected Temperature downscaling algorithm (previous versions introduced -ve bias to houly air temp)
!
!Cabledyn05 (VH 13-JUN-08)
! B-Horizon soil properties set as elements of CABLE "soil" structure and passed to "call_soil_Ross".
! Depth of bottom soil layer modified to be consistent with depth of B-Horizon
!"jt" set to be consistent with depth of A Horizon
! soil layers adjusted for comparison with soil-moisture measurements
!  Noted that Whole-of-run discharge shows some dependence on soil layers: need to investigate this further.
!
!Cabledyn06P (VH 17-JUN-08)
! Included multipliers for Vcmax and Rootbeta
!
!CableDyn07a (VH 14-JUL-2008)
!  New obs type: AAPh, ZZPh, which have an 3rd dimension, specifying time of day. This allows multiple LSTs per day.
!
!Cabledyn08 (VH 25-JUL-2008)
! Add phiHh to AAPh and ZZPh
!
!CableDyn08 (VH 29-JUL-2008)
!add kth as to soil structure in CABLE (allows output of ccnsw (soil themal conductivity).
! Replaced soil-snow formulation of kth with Campbell eq (4.20), p.32, "Soil Physics with Basic"
!
!Cabledyn09 (VH 1-4-AUG-2008)
!Added Rneth, phiEh, NEEh to AAPh and ZZPh
! Added NEE, phiH, phiE and phiRnet to AAP and ZZP
!Added phiGh to AAPh and ZZPh
!
!Cabledyn11 (VH 16-OCT-1008)
! Insert soil-moisture variables from Cabledyn10
! Option to use site met data (UseLocalMetDataFlag)
!
!Cabledyn11 (VH&BS 22-DEC-2008)
! UseLAIFlag for direct LAI
! monthly LAI added to VV parameter list
!
!CableDyn12 (VH 24-FEB-2009)
!Upated to soil-litter05 (uses Lapack instead of IMSL; new upper BC calc; optional isotopes) 
! Model can be forced using local met data even when there is no gridded data available
!
!CableDyn12P (VH 1-APR-2009)
! corrected reading in of hourly observables (AAPh): previously, time was out by 1 h
! rtsoilmult and swilt1mult added to UU list
! Updated PEST subroutine to allow parameter optimisation against variables in AAPh
! Tsoilh added to AAPh and ZZPh
!
!CableDyn13 (VH 2-APR-2009)
! add clittequil to VV
! inlcude Ross soil model as switchable option
!
!CableDyn14 (VH 28-APR-2009)
! Changed paths to be unix-compatible (i.e. '/' instead of '\')
! Udated definitions of Wre11:Wrel6; Tsoil1:Tsoil6
! Modified, InitStores to use inifile if it exists, otherwise default value when XXInitFileFlag==1
!
!CableDyn15 (VH 21-MAY-2009)
! Upgrade to CABLE v1.4b (problems with leaf temperature convergence in this version of CABLE)
! Correct so that model runs when UseLocalMetDataFlag=1
! Energy closure for soil and vegetation
!
!CableDyn16 (VH 26-MAY-2009)
! Upgrade to CABLE v1.4d
! add theta1-theta6 and wcol to XX
!
! CableDyn16 (VH 16-JUN-2009)
! CABLE default soil and veg param files passed from Cabledyn to CABLE's default_params
! add FNEE to FF
! include correction to clear-sky l-w radiation
!
! Cabledyn17 (VH 20-JUL-2009:21-JUL-2009)
! Remove evapfbl term (and subsequent modification to latent heat flux) from cable_canopy
! Introduce Lai and Katul (2000) parameterisation for root efficiency function and use to create alternative fwsoil
! Introduce woody (_w) and grassy (_g) veg parameters and reat in from control file
! Restrict soil moisture predictions and obs to 0-8 cm and 0-90 cm
! Add GPP and NPP to ZZP and AAP
! Fix up water balance
! Add fws to DD (factor for limiting stomatal conductance due to soil moisture deficit)
!
!Cabledyn17 (VH 06-AUG-2009)
! New SensTestFlag = 4: permits calculation of relative changes in fluxes when UU varied between UUmin and UUmax
! Reverted from logKsatmult to Ksatmult 
!
! Cabledyn17 (VH 04-SEP-2009)
! Modified WritePestOutput to write soil moisture anomalies, calculate residuals and weights and write weights to ActualObs.dat	
!
!Cabledyn17 (VH 07-SEP-2009)
! Exchanged rtsoilmult for Ksat2mult
! Exchanged Ksatmult for Ksat1mult
!
!Cabledyn20 (VH 10-NOV-2009)
! Incorporated updated soil-litter-iso
! 10 soil layers
! Incorporated updated cable-canopy, incl new root extraction
!
!Cabledyn21 (VH 05-JAN-2010)
! Included sm0030, sm3060 and sm6090 in point scale observations and observables
!
!Cabledyn21 (VH 11-JAN-2010)
! New formulation for Ksat, based on Timlin et al., 1999
!
! Cabledyn 22 (VH)
! set a1 and ds0 as spatially constant params
!
! Cabledyn23 (VH 14-Dec-2010)
! write out hourly Tsoil, Ta, theta
! Set MapOutputFlag to 4 if day lies between StartDate_daily_output and EndDate_daily_output
! introduce integer vector Ihour_output to specify which hours to write out
!
! Cabledyn24 (VH March 2011)
! Introduced cable-working code, last modified 21/3/2011
! Modified InitCable and call to CABLE to allow use of new CABLE code
! Tiling capability introduced
! separate params for Woody and Grassy tiles
! UseLAIFlag=2 => separate LAI for grassy and woody components
! CableParamFlag        (0,1,2) = (single tile estimate, woody/grassy tile, CABLE defaults) 
! Readobs modified to get monthly remote sensing inputs
! New datefunctions to check if months are equal
! New cable-sli-main and cable-sli-solve to allow root extraction form multiple tiles from single soil column
! Root-extraction now in cable-sli-main
! Fixed bug in set_par (cable-_sli_utils)
! New Array RR contains remote sensing drivers
! RR files opened and initialised in InitModel
! RR updated in ReadRR (allows for climatologies and time-series)
! File for input of atmospheric CO2, opened and initialised in InitModel
! CO2A updated in ReadRR
! Create OutputMap_tmp to output monthly state variables to tmp file in InitStores for restart
! Augment DD to accommodate variables required to drive CASA-CNP
! InitMod and ReadRR modified to allow climatolgical data to be used before and after RR time-seires input when 
! time-series do not span full run length
!
! Cabledyn25 (VH June 2011)
! UseLocalMetDataFlag=2: repeat met across all pixels
!
! Cabledyn25 (VH Sept 2011)
! HourlyOutputFlag=2: write cable subdiurnal output to netcdf
! delT : change to variable for CABLE time step (in s), and replace delT elsewhere in code with 1.0 (days)
! FWPT diagnostic enabled
!
! Cabledyn 25a (VH May 2012)
! if no vph available, set vph =esat(Tmin)
! spin-up option
!
! Cabledyn 25b (VH JUN 2012)
! converting screen T and Q to values at ref height
!#******************************************************************************
!*******************************************************************************

MODULE TypeDef
!------------------------------------------------------------------------------- 
! PRB: 18-06-2000
! This module explicitly defines the sizes of variable types
!-------------------------------------------------------------------------------
  implicit none
  save
! Define integer kind parameters to accommodate the range of numbers usually 
! associated with 4, 2, and 1 byte integers. 
  integer,parameter :: i4b = selected_int_kind(9) 
  integer,parameter :: i2b = selected_int_kind(4)
  integer,parameter :: i1b = selected_int_kind(2)
! Define single and double precision real kind parameters: 
! * Kind(1.0)   defines sp as the machine's default size for single precision
! * Kind(1.0d0) defines dp as the machine's default size for double precision
  integer,parameter :: sp  = kind(1.0)    
  integer,parameter :: dp  = kind(1.0d0)
! lgt is set to the default kind required for representing logical values. 
  integer,parameter :: lgt = kind(.true.)

END MODULE TypeDef

!*******************************************************************************
!*******************************************************************************

MODULE Constants
!-------------------------------------------------------------------------------
! * Physical constants and parameters
! * MRR, 2004
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
save
real(sp),parameter:: Pi    = 3.141592653589793238462643383279502884197
real(sp),parameter:: TwoPi = 6.283185307179586476925286766559005768394
real(sp),parameter:: Rgas      = 8.3143   ! universal gas constant    [J/mol/K]
real(sp),parameter:: RMA       = 0.02897  ! molecular wt dry air      [kg/mol]
real(sp),parameter:: RMW       = 0.018016 ! molecular wt of water     [kg/mol]
real(sp),parameter:: RMC       = 0.012000 ! atomic wt of C            [kg/mol]
real(sp),parameter:: SBoltz    = 5.67e-8  ! Stefan-Boltzmann constant [W/m2/K4]
real(sp),parameter:: Rlat      = 44140.0  ! lat heat evap H2O at 20C  [J/molW]
real(sp),parameter:: Capp      = 29.09    ! isobaric spec heat air    [J/molA/K]
real(sp),parameter:: RhoW      = 55506.0  ! liquid water density      [molW/m3]
real(sp),parameter:: ViscW     = 1.002    ! water viscosity, 25 C     [Pa s]
real(sp),parameter:: Grav      = 9.81     ! gravity acceleration      [m/s2]
real(sp),parameter:: VonK      = 0.4      ! von Karman constant       [-]
real(sp),parameter:: SecDay    = 86400.0  ! seconds/day               [-]
real(sp),parameter:: QMJSolar  = 2.0      ! (molPAR)/(MJsolar)        [molQ/MJ]
real(sp),parameter:: SolCMJday = 118.4    ! solar constant            [MJ/m2/d]
real(sp),parameter:: SolCWm2   = 1370.0   ! solar constant            [W/m2]
real(sp),parameter:: REarth    = 6.37e6   ! average earth radius      [m]
REAL(sp), PARAMETER :: rlam =  2.5104e6  ! latent heat of evaporation (J kg-1)  ! edit vh 13/02/08
END MODULE Constants

!*******************************************************************************
!*******************************************************************************

MODULE DateFunctions
!-------------------------------------------------------------------------------
! Date arithmetic (PRB, 12-07-2000; MRR, 11-oct-05)
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
public
save
! Define a type for Australian-style dates built from 3 integer fields.
! Assign as Date = dmydate(iday,imth,iyear)
! Access components as iday=Date%Day, imth=Date%Month, iyear=Date%Year
type dmydate
  integer(i4b):: Day
  integer(i4b):: Month
  integer(i4b):: Year
end type dmydate
! Define interfaces for the +, -, .eq., .ne., .lt., .gt., .le., .ge. operators to:
! * add, subtract a number of days (integer) from a dmydate, with AddDay, SubDay; 
! * logically compare two dates, with EQDates, NEDates etc.
! These explicit interfaces are required.
interface operator (+)
  module procedure AddDay
end interface
interface operator (-)
  module procedure SubDay
end interface
interface operator (.eq.)
  module procedure EQDates
end interface
interface operator (.ne.)
  module procedure NEDates
end interface
interface operator (.lt.)
  module procedure LTDates
end interface
interface operator (.gt.)
  module procedure GTDates
end interface
interface operator (.le.)
  module procedure LEDates
end interface
interface operator (.ge.)
  module procedure GEDates
end interface
!-------------------------------------------------------------------------------

CONTAINS 

!*******************************************************************************
  
  function AddDay(Today,Days2Add)
!-------------------------------------------------------------------------------
! Extends the '+' operator to enable adding days to a date
!-------------------------------------------------------------------------------
  implicit none
  type (dmydate)              :: AddDay    ! Date with days added (function name)
  type (dmydate), intent(in)  :: Today     ! Current date
  integer(i4b),   intent(in)  :: Days2Add  ! Days to add to current date
  real(dp)                    :: JDay      ! Julian Day
!-------------------------------------------------------------------------------
  JDay = JulianDay(Today)
  JDay = JDay + float(Days2Add)
  AddDay = GregDate(JDay)
  end function AddDay

!*******************************************************************************

  function SubDay(Today,Days2Sub)
!-------------------------------------------------------------------------------
! Extends the '-' operator to enable subtracting days from a date
!-------------------------------------------------------------------------------
  implicit none
  type (dmydate)              :: SubDay    ! Date with days added (function name)
  type (dmydate), intent(in)  :: Today     ! Current date
  integer(i4b),   intent(in)  :: Days2Sub  ! Days to add to current date
  real(dp)                    :: JDay      ! Julian Day
!-------------------------------------------------------------------------------
  JDay = JulianDay(Today)
  JDay = JDay - float(Days2Sub)
  SubDay = GregDate(JDay)
  end function SubDay

!*******************************************************************************

  function EQDayofYear(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.eq.' operator to enable comparison of two dates for equality.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: EQDayofYear   ! Result of equality test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------
  EQDayofYear = .true.
  if (Date1%day   .ne. Date2%day   .or. &
      Date1%month .ne. Date2%month ) EQDayofYear = .false.
  end function EQDayofYear
!*******************************************************************************

  function EQDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.eq.' operator to enable comparison of two dates for equality.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: EQDates   ! Result of equality test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------
  EQDates = .true.
  if (Date1%day   .ne. Date2%day   .or. &
      Date1%month .ne. Date2%month .or. &
	  Date1%year  .ne. Date2%year) EQDates = .false.
  end function EQDates

!*******************************************************************************

  function NEDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.ne.' operator to enable comparison of two dates for inequality.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: NEDates   ! Result of inequality test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------  
  NEDates = .not. EQDates(Date1,Date2)
  end function NEDates

!*******************************************************************************

  function LTDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.lt.' operator to enable comparison of two dates for Date1 < Date2.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: LTDates   ! Result of LT test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!------------------------------------------------------------------------------- 
  LTDates = (JulianDay(Date1).lt.JulianDay(Date2))
  end function LTDates

!*******************************************************************************

  function GTDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.gt.' operator to enable comparison of two dates for Date1 > Date2.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: GTDates   ! Result of GT test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------  
  GTDates = (JulianDay(Date1).gt.JulianDay(Date2))
  end function GTDates

!*******************************************************************************

  function LEDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.le.' operator to enable comparison of two dates Date1 <= Date2.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: LEDates   ! Result of LE test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------  
  LEDates = (LTDates(Date1,Date2) .or. EQDates(Date1,Date2))
  end function LEDates

!*******************************************************************************

  function GEDates(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.ge.' operator to enable comparison of two dates Date1 >= Date2.
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: GEDates   ! Result of GE test between dates
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------  
  GEDates = (GTDates(Date1,Date2) .or. EQDates(Date1,Date2))
  end function GEDates
!*******************************************************************************

  function EQMonths(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.eq.' operator to enable comparison of two time-series months for equality.
! VH, 25-03-11
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: EQMonths   ! Result of equality test between months
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------
  EQMonths = .true.
  if (Date1%month .ne. Date2%month .or. &
	  Date1%year  .ne. Date2%year) EQMonths = .false.
  end function EQMonths
!*******************************************************************************

  function EQMonthofYear(Date1,Date2)
!-------------------------------------------------------------------------------
! Extends the '.eq.' operator to enable comparison of two months of year for equality.
! VH, 25-03-11
!-------------------------------------------------------------------------------
  implicit none
  logical                     :: EQMonthofYear   ! Result of equality test between months
  type (dmydate), intent(in)  :: Date1     ! First date for logical comparison
  type (dmydate), intent(in)  :: Date2     ! Second date for logical comparison
!-------------------------------------------------------------------------------
  EQMonthofYear = .true.
  if (Date1%month .ne. Date2%month) EQMonthofYear = .false.
  end function EQMonthofYear
!*******************************************************************************

  function JulianDay(GregDate)
!-------------------------------------------------------------------------------
! Returns a real(sp) Julian Day when given a Gregorian Date.
! Adapted from Date Algorithms of Peter Baum:
! http://vsg.cape.com/~pbaum/date/date0.htm
!-------------------------------------------------------------------------------
  implicit none
  real(dp)                   :: JulianDay
  type (dmydate), intent(in) :: GregDate
  real(dp)                   :: D,M,Y
!-------------------------------------------------------------------------------
  D = dble(GregDate%Day)
  M = dble(GregDate%Month)
  Y = dble(GregDate%Year)
  if (M.lt.3.0_8) then 
    M = M + 12.0_8 
    Y = Y - 1.0_8 
  end if 
  JulianDay = D + int((153.0_8 * M - 457.0_8)/5.0_8) + (365.0_8*Y) +     &
              floor(Y/4.0_8) - floor(Y/100.0_8) + floor(Y/400.0_8) + 1721118.5_8
  end function JulianDay

!*******************************************************************************

  function GregDate(JulianDay)
!-------------------------------------------------------------------------------
! Returns a Gregorian Date when given a Julian Day 
! Modified from Date Algorithms of Peter Baum 
! http://vsg.cape.com/~pbaum/date/date0.htm
!-------------------------------------------------------------------------------
  implicit none
  type (dmydate)             :: GregDate
  real(dp), intent(in)       :: JulianDay
  real(dp)                   :: D,M,Y
  real(dp)                   :: Z,R,A,G,B,C
!-------------------------------------------------------------------------------
  Z = floor(JulianDay - 1721118.5_8) 
  R = JulianDay - 1721118.5_8 - Z
  G = Z - 0.25_8 
  A = floor(G/36524.25_8) 
  B = A - floor(A/4.0_8) 
  Y = floor((B+G)/365.25_8) 
  C = B + Z - floor(365.25_8*Y) 
  M = int(((5.0_8*C) + 456.0_8) / 153.0_8) 
  D = C - int((153.0_8 * M - 457.0_8) / 5.0_8) + R  
  if (M .gt. 12.0_8) then 
    Y = Y + 1.0_8 
    M = M - 12.0_8 
  end if
  GregDate = dmydate(int(D),int(M),int(Y))  ! Keep truncated day number only
  end function GregDate

!*******************************************************************************

  function LeapDay(Year)
!-------------------------------------------------------------------------------
! Returns 1 if leap year, 0 if not.  Add it to the length of February.
!-------------------------------------------------------------------------------
  implicit none
  integer(i4b)               :: LeapDay  
  integer(i4b), intent(in)   :: Year
!-------------------------------------------------------------------------------
  if (mod(Year,4).ne.0) then
    LeapDay = 0
  else if (mod(Year,400).eq.0) then
    LeapDay = 1
  else if (mod(Year,100).eq.0) then
    LeapDay = 0
  else
    LeapDay = 1
  end if
  end function LeapDay

!*******************************************************************************

function DaysInMonth(Date)
!-------------------------------------------------------------------------------
! returns number of days in current month
! MRR, 12-oct-2005: return 0 if month not legal
!-------------------------------------------------------------------------------
type(dmydate),intent(in):: Date
integer(i4b)          :: DaysInMonth  
integer(i4b),parameter:: MonthDays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
!-------------------------------------------------------------------------------
if (Date%Month >= 1 .and. Date%Month <= 12) then
  if (Date%Month.eq.2) then
      DaysInMonth = MonthDays(2) + LeapDay(Date%Year)
  else 
    DaysInMonth = MonthDays(Date%Month)
  end if
else
  DaysInMonth = 0
end if
end function DaysInMonth

!*******************************************************************************

  function YearDay(Date)
!-------------------------------------------------------------------------------
  type(dmydate), intent(in)   :: Date
  integer(i4b)                :: YearDay 
  integer(i4b), dimension(12) :: MonthDays
!-------------------------------------------------------------------------------
  MonthDays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  MonthDays(2) = 28 + LeapDay(Date%Year)
  if (Date%Month.eq.1) then
    YearDay = Date%Day
  else
    YearDay = sum(MonthDays(1:Date%Month-1)) + Date%Day
  end if
  end function YearDay

!*******************************************************************************

FUNCTION DayDifference (Date1,Date0)
!-------------------------------------------------------------------------------
! Returns Date1 - Date0 in integer days.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
type(dmydate),intent(in):: Date1, Date0
integer(i4b):: DayDifference
!-------------------------------------------------------------------------------
DayDifference = nint(JulianDay(Date1) - JulianDay(Date0))
END FUNCTION DayDifference

!*******************************************************************************

FUNCTION LegalDate (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is legal (test D,M only), otherwise .false.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
type(dmydate),intent(in):: Date
logical(lgt):: LegalDate
!-------------------------------------------------------------------------------
LegalDate = .false.
if ( (Date%month >= 1 .and. Date%month <= 12) .and.         &   ! check month
     (Date%day >= 1 .and. Date%day <= DaysInMonth(Date)) )  &   ! check day
   LegalDate = .true.
END FUNCTION LegalDate

!*******************************************************************************

FUNCTION EndMonth (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is last day of month, otherwise .false.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
type(dmydate),intent(in):: Date
logical(lgt):: EndMonth
!-------------------------------------------------------------------------------
EndMonth = .false.
if (Date%day == DaysInMonth(Date)) EndMonth = .true.
END FUNCTION EndMonth
!*******************************************************************************

FUNCTION BeginMonth (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is first day of month, otherwise .false.
! VH, 25-03-11
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
type(dmydate),intent(in):: Date
logical(lgt):: BeginMonth
!-------------------------------------------------------------------------------
BeginMonth = .false.
if (Date%day == 1) BeginMonth = .true.
END FUNCTION BeginMonth

!*******************************************************************************

FUNCTION EndYear (Date)
!-------------------------------------------------------------------------------
! Returns .true. if date is last day of year, otherwise .false.
! MRR, 11-oct-05
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
type(dmydate),intent(in):: Date
logical(lgt):: EndYear
!-------------------------------------------------------------------------------
EndYear = .false.
if (Date%day == 31 .and. Date%month == 12) EndYear = .true.
END FUNCTION EndYear

!*******************************************************************************

END MODULE DateFunctions

!*******************************************************************************
!*******************************************************************************

MODULE Utils
!-------------------------------------------------------------------------------
! * Utility subroutines from library Milib
!-------------------------------------------------------------------------------
CONTAINS

!*******************************************************************************

ELEMENTAL FUNCTION Qsatf(TC,Pmb)
!-------------------------------------------------------------------------------
! At temperature TC [deg C] and pressure Pmb [mb], return saturation specific
! humidity Qsatf [kgH2O / kg moist air] from Teten formula.
! MRR, xx/1987
! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
!                 just like intrinsic functions.
! MRR, 28-mar-05: Use RMA, RMW from MODULE Constants
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
implicit none
real(sp),intent(in) :: TC, Pmb      ! temp [deg C], pressure [mb]
real(sp):: Qsatf                    ! saturation specific humidity
real(sp):: TCtmp                    ! local
real(sp),parameter:: A = 6.106      ! Teten coefficients
real(sp),parameter:: B = 17.27      ! Teten coefficients
real(sp),parameter:: C = 237.3      ! Teten coefficients
!-------------------------------------------------------------------------------
TCtmp = tc                          ! preserve TC
if (TCtmp.gt.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
if (TCtmp.lt.-40.0) TCtmp = -40.0
Qsatf = (RMW/RMA) * (A*EXP(B*TCtmp/(C+TCtmp))) / Pmb
                                    ! sat specific humidity (kg/kg)
END FUNCTION Qsatf

!*******************************************************************************

ELEMENTAL FUNCTION Esatf(TC)
!-------------------------------------------------------------------------------
! At temperature TC [deg C], return saturation water vapour pressure Esatf [mb] 
! from Teten formula.
! MRR, xx/1987
! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
!                 just like intrinsic functions.
! MRR, 12-mar-02: Convert Qsatf (specific humidity routine) to Esatf
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
real(sp), intent(in):: TC           ! temp [deg C]
real(sp):: Esatf                    ! saturation vapour pressure [mb]
real(sp):: TCtmp                    ! local
real(sp),parameter:: A = 6.106      ! Teten coefficients
real(sp),parameter:: B = 17.27      ! Teten coefficients
real(sp),parameter:: C = 237.3      ! Teten coefficients
!-------------------------------------------------------------------------------
TCtmp = TC                          ! preserve TC
if (TCtmp.gt.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
if (TCtmp.lt.-40.0) TCtmp = -40.0
Esatf = A*EXP(B*TCtmp/(C+TCtmp))    ! sat vapour pressure (mb)
 
END FUNCTION Esatf

!*******************************************************************************

ELEMENTAL FUNCTION Epsif(TC,Pmb)
!-------------------------------------------------------------------------------
! At temperature TC [deg C] and pressure Pmb [mb], return 
! epsi = (RLAM/CAPP) * d(sat spec humidity)/dT [(kg/kg)/K], from Teten formula.
! MRR, xx/1987, 27-jan-94
! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
!                 just like intrinsic functions.
! MRR, 28-mar-05: use Rlat [J/molW], Capp [J/molA/K] from MODULE Constants, 
!                 to ensure consistency with other uses of Rlat, Capp
! MRR, 28-mar-05: Remove dependence of Rlat (latent heat vaporisation of water)
!                 on temperature, use value at 20 C
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
implicit none
real(sp),intent(in):: TC, Pmb       ! temp [deg C], pressure [mb]
real(sp):: Epsif                    ! epsi
real(sp):: TCtmp, ES, dESdT         ! local
real(sp),parameter:: A = 6.106      ! Teten coefficients
real(sp),parameter:: B = 17.27      ! Teten coefficients
real(sp),parameter:: C = 237.3      ! Teten coefficients
!-------------------------------------------------------------------------------
TCtmp = TC                          ! preserve TC
if (TCtmp.gt.100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
if (TCtmp.lt.-40.0) TCtmp = -40.0
ES    = A*EXP(B*TCtmp/(C+TCtmp))    ! sat vapour pressure
dESdT = ES*B*C/(C+TCtmp)**2         ! d(sat VP)/dT: (mb/K)
Epsif = (Rlat/Capp) * dESdT / Pmb   ! dimensionless (ES/Pmb = molW/molA)

END FUNCTION Epsif

!*******************************************************************************

SUBROUTINE AstronDailySolar (YearFrac, LatDeg,   DecDeg, DayLtFrac, SolarNorm)
!-------------------------------------------------------------------------------
! Calculate daily solar irradiance (normalised by solar constant), daylight 
! fraction and solar declination angle, for given latitude and year day
! In:  * YearFrac  = year day as year fraction in [0,1] (scalar)
!      * LatDeg    = latitude (degrees) (vector of many points)
! Out: * DecDeg    = Solar declination (deg, +23.5 on 22 June, -23.5 on 22 Dec)
!      * DayLtFrac = Daylight Fraction (0 to 1)
!      * SolarNorm = Daily solar irradiance without atmosphere, normalised 
!                      by solar constant, calculated from solar geometry
! * DecDeg    = func(YearFrac), hence scalar
! * DayLtFrac = func(YearFrac, LatDeg), hence vector
! * SolarNorm = func(YearFrac, LatDeg), hence vector
! HISTORY:
! * 04-sep-03 (MRR): Program SubDiurnalWeather
! * 12-dec-04 (MRR): Subroutine AstronDailySolar created from SubDiurnalWeather,
!                    for program WaterEquil
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
implicit none
! Global variables
real(sp),intent(in) :: YearFrac     ! day of the year as year fraction (0 to 1)
real(sp),intent(in) :: LatDeg(:)    ! latitude [degs, neg for SH]
real(sp),intent(out):: DecDeg       ! declination [deg] (+23.5 deg on 22 June)
real(sp),intent(out):: DayLtFrac(:) ! daylight fraction (dawn:dusk)    (0 to 1)
real(sp),intent(out):: SolarNorm(:) ! (daily solar)/(solar const): geometry
! Local variables
real(sp):: YearRad, DecRad
real(sp):: LatDeg1(size(LatDeg)), LatRad(size(LatDeg)), &
           TanTan(size(LatDeg)), HDLRad(size(LatDeg))
!-------------------------------------------------------------------------------

if (YearFrac < 0 .or. YearFrac > 1) stop "AstronDailySolar: illegal YearFrac"
if (any(abs(LatDeg) > 90.0))        stop "AstronDailySolar: illegal LatDeg"
! DecRad = declination in radians (+23.5 deg on 22 June, -23.5 deg on 22 Dec)
YearRad = 2.0*Pi*YearFrac               ! day of year in radians
DecRad  = 0.006918 - 0.399912*cos(YearRad) + 0.070257*sin(YearRad)  &
          - 0.006758*cos(2.0*YearRad) + 0.000907*sin(2.0*YearRad)   & 
          - 0.002697*cos(3.0*YearRad) + 0.001480*sin(3.0*YearRad)
                                        ! Paltridge and Platt eq [3.7]
DecDeg = (180.0/Pi) * DecRad            ! Declination in degrees

! Daylength: HDLRad = Half Day Length in radians (dawn:noon = noon:dusk)
LatDeg1 = sign(min(abs(LatDeg),89.9), LatDeg)   ! avoid singularity at pole
LatRad  = LatDeg1*Pi/180.0                      ! latitude in radians
TanTan  = -tan(LatRad)*tan(DecRad)
where (TanTan .le. -1.0)
  HDLRad = Pi                           ! polar summer: sun never sets
elsewhere (TanTan .ge. 1.0)
  HDLRad = 0.0                          ! polar winter: sun never rises
elsewhere
  HDLRad = acos(TanTan)                 ! Paltridge and Platt eq [3.21]
end where                               ! (HDLRad = their capital PI)
DayLtFrac = 2.0*HDLRad / (2.0*Pi)       ! Daylight fraction (dawn:dusk)

! Daily solar irradiance without atmosphere, normalised by solar constant,
! with both energy fluxes in MJ/m2/day, calculated from solar geometry:
SolarNorm =  &                          ! Paltridge and Platt eq [3.22]
  (HDLRad*sin(LatRad)*sin(DecRad) + cos(LatRad)*cos(DecRad)*sin(HDLRad)) / Pi

END SUBROUTINE AstronDailySolar

!*******************************************************************************

SUBROUTINE AstronSZAtoTime (SZARad, YearFrac, LatDeg, TimeHr)
!-------------------------------------------------------------------------------
! Find local time TimeHr from solar zenith angle (SZARad) and time of year (YearFrac)
! HISTORY:
! * 13-jul-05 (MRR): written
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
implicit none
! Global variables
real(sp),intent(in) :: SZARad(:)    ! solar zenith angle (rad) (0=zenith, Pi/2=horizon)
real(sp),intent(in) :: YearFrac     ! day of the year as year fraction (0 to 1)
real(sp),intent(in) :: LatDeg(:)    ! latitude [degs, neg for SH]
real(sp),intent(out):: TimeHr(:)    ! local time (0.0 to 24.0, 12.0=noon)
! Local variables
real(sp):: YearRad, DecRad
real(sp):: LatRad(size(LatDeg)), CosThr(size(LatDeg)), SZARad1(size(LatDeg))
!-------------------------------------------------------------------------------

if (YearFrac < 0 .or. YearFrac > 1) stop "AstronSZAtoTime: illegal YearFrac"
SZARad1 = abs(SZARad)                       ! keep SZARad > 0
where (SZARad1 > Pi/2.0) SZARad1 = Pi/2.0   ! and <= Pi/2
! DecRad = declination in radians (+23.5 deg on 22 June, -23.5 deg on 22 Dec)
YearRad = 2.0*Pi*YearFrac               ! day of year in radians
DecRad  = 0.006918 - 0.399912*cos(YearRad) + 0.070257*sin(YearRad)  &
          - 0.006758*cos(2.0*YearRad) + 0.000907*sin(2.0*YearRad)   & 
          - 0.002697*cos(3.0*YearRad) + 0.001480*sin(3.0*YearRad)
                                        ! Paltridge and Platt eq [3.7]
LatRad  = LatDeg*Pi/180.0
CosThr  = (cos(SZARad1) - sin(DecRad)*sin(LatRad)) / (cos(DecRad)*cos(LatRad))
                                        ! Thr = Hour Angle from 1200 (rad)
                                        ! Paltridge and Platt eq [3.4]
CosThr  = min(CosThr, 1.0)              ! don't allow invalid CosThr
CosThr  = max(CosThr, -1.0)
TimeHr  = acos(CosThr) *24.0/(2.0*Pi) + 12.0
                                        ! Time (Hr, 0 to 24)
END SUBROUTINE AstronSZAtoTime

!*******************************************************************************

SUBROUTINE ComSkp (iunit)
!-------------------------------------------------------------------------------
! MRR, 5-AUG-83
! SKIPS COMMENT LINES IN CONTROL DATA FILE, STOPS IF EOF FOUND
! PRB, Sep-99, F90: Use preferred alternative to 'do while' [M&R96:12.3.2] 
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
character (len=1)      :: com
integer(i4b),intent(in):: iunit
!-------------------------------------------------------------------------------

com = '!'
do 
  if (com.ne.'!' .and. com.ne.'C' .and. com.ne.'c') exit ! Exit loop
  read (iunit,'(A1)',end=99) com
end do
backspace iunit
return
99 stop 'CONTROL FILE EOF'

END SUBROUTINE ComSkp

!*******************************************************************************

ELEMENTAL FUNCTION AnStep(x, xs, dx)
!-------------------------------------------------------------------------------
! Function AnStep(x,xs,dx), with small mod(dx) (say < 0.1) is an analytic 
! smooth approximation to Heaviside step function, using hyperbolic tangent.
! For dx > 0, function being approximated is the upward step
! AnStep = 0.0 for x < xs
! AnStep = 1.0 for x > xs
! For dx < 0, the function being approximated is a downward step.
! MRR, 03-jul-04
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
real(sp),intent(in) :: x, xs, dx
real(sp)            :: AnStep
!-------------------------------------------------------------------------------
AnStep = 0.5 * ( 1.0 + tanh((x-xs)/dx) )
END FUNCTION AnStep

!*******************************************************************************

ELEMENTAL FUNCTION AnMin(x1, x2, a)
!-------------------------------------------------------------------------------
! Returns an analytic smooth approximation to the minimum of x1 and x2.
! Tightness of approximation is determined by a: mod(a) should be 10 or larger.
! Note: formulation assumes x1 > 0, x2 > 0.
! MRR, 07-jul-04
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
real(sp),    intent(in) :: x1, x2
integer(i4b),intent(in) :: a
real(sp)                :: AnMin, x1a, x2a
!-------------------------------------------------------------------------------
x1a   = x1**a
x2a   = x2**a
AnMin = ( (x1a*x2a)/(x1a+x2a) )**(1.0/real(a))
END FUNCTION AnMin

END MODULE Utils

!*******************************************************************************
!*******************************************************************************

MODULE PointerModule
!-------------------------------------------------------------------------------
! * This module defines and assigns pointers to specific variables in WaterDynStep
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! * Pointers for time variables (array TTime)
real(sp),pointer:: Time1y       ! [y]           ! 01 run time at end of step
real(sp),pointer:: TTDay        ! [d]           ! 02 current date: day
real(sp),pointer:: TTMonth      ! [mth]         ! 03 current date: month
real(sp),pointer:: TTYear       ! [y]           ! 04 current date: year
real(sp),pointer:: TTEndMth     ! [-]           ! 05 (0,1) = last day of month (N,Y)
! * Pointers for stores (array XX)
real(sp),pointer:: WRelA(:)     ! [0,1]         ! 01 rel soil water (layer A)
real(sp),pointer:: WRelB(:)     ! [0,1]         ! 02 rel soil water (layer 2B)
real(sp),pointer:: CLea(:)      ! [molC/m2]     ! 03 live leaf C
real(sp),pointer:: WRel1(:)     ! [0,1]         ! 04 rel soil water (layer 1)
real(sp),pointer:: WRel2(:)     ! [0,1]         ! 05 rel soil water (layer 2)
real(sp),pointer:: WRel3(:)     ! [0,1]         ! 06 rel soil water (layer 3)
real(sp),pointer:: WRel4(:)     ! [0,1]         ! 07 rel soil water (layer 4)
real(sp),pointer:: WRel5(:)     ! [0,1]         ! 08 rel soil water (layer 5)
real(sp),pointer:: WRel6(:)     ! [0,1]         ! 09 rel soil water (layer 6)
real(sp),pointer:: WRel7(:)     ! [0,1]         ! 10 rel soil water (layer 7)
real(sp),pointer:: WRel8(:)     ! [0,1]         ! 11 rel soil water (layer 8)
real(sp),pointer:: WRel9(:)     ! [0,1]         ! 12 rel soil water (layer 9)
real(sp),pointer:: WRel10(:)     ! [0,1]         ! 13 rel soil water (layer 10)
real(sp),pointer:: Tsoil1(:)    ! [0,1]         ! 14 soil T (layer 1)
real(sp),pointer:: Tsoil2(:)    ! [0,1]         ! 15 soil T (layer 2)
real(sp),pointer:: Tsoil3(:)    ! [0,1]         ! 16 soil T (layer 3)
real(sp),pointer:: Tsoil4(:)    ! [0,1]         ! 17 soil T (layer 4)
real(sp),pointer:: Tsoil5(:)    ! [0,1]         ! 18 soil T (layer 5)
real(sp),pointer:: Tsoil6(:)    ! [0,1]         ! 19 soil T (layer 6)
real(sp),pointer:: Tsoil7(:)    ! [0,1]         ! 20 soil T (layer 7)
real(sp),pointer:: Tsoil8(:)    ! [0,1]         ! 21 soil T (layer 8)
real(sp),pointer:: Tsoil9(:)    ! [0,1]         ! 22 soil T (layer 9)
real(sp),pointer:: Tsoil10(:)    ! [0,1]         ! 23 soil T (layer 10)
real(sp),pointer:: cplant1(:)   ! [gC/m2]       ! 24 leaf C
real(sp),pointer:: cplant2(:)   ! [gC/m2]       ! 25 wood C
real(sp),pointer:: cplant3(:)   ! [gC/m2]       ! 26 root C
real(sp),pointer:: csoil1(:)    ! [gC/m2]       ! 27 soil microbial biomass
real(sp),pointer:: csoil2(:)    ! [gC/m2]       ! 28 soil fast C pool
real(sp),pointer:: cansto(:)    ! [mm]           ! 29 water stored in canopy
real(sp),pointer:: theta1(:)     ! [0,1]         ! 30 vol soil water (layer 1)
real(sp),pointer:: theta2(:)     ! [0,1]         ! 31 vol soil water (layer 2)
real(sp),pointer:: theta3(:)     ! [0,1]         ! 32 vol soil water (layer 3)
real(sp),pointer:: theta4(:)     ! [0,1]         ! 33 vol soil water (layer 4)
real(sp),pointer:: theta5(:)     ! [0,1]         ! 34 vol soil water (layer 5)
real(sp),pointer:: theta6(:)     ! [0,1]         ! 35 vol soil water (layer 6)
real(sp),pointer:: theta7(:)     ! [0,1]         ! 36 vol soil water (layer 7)
real(sp),pointer:: theta8(:)     ! [0,1]         ! 37 vol soil water (layer 8)
real(sp),pointer:: theta9(:)     ! [0,1]         ! 38 vol soil water (layer 9)
real(sp),pointer:: theta10(:)     ! [0,1]         ! 39 vol soil water (layer 10)
real(sp),pointer:: wcol(:)     ! [m]         ! 240 total soil water column (layers 1-10)

! * Pointers for fluxes and dXXdt (array FF)
real(sp),pointer:: dWcolAdt(:)  ! [0,1]/d       ! 01 dXX/dt:  soil water (A)
real(sp),pointer:: dWcolBdt(:)  ! [0,1]/d       ! 02 dXX/dt:  soil water (B)
real(sp),pointer:: FWPrec(:)    ! [m/d]         ! 03 precipitation
real(sp),pointer:: FWTraA(:)    ! [m/d]         ! 04 transpiration (layer 1)
real(sp),pointer:: FWTraB(:)    ! [m/d]         ! 05 transpiration (layer 2)
real(sp),pointer:: FWSoil(:)    ! [m/d]         ! 06 soil evaporation
real(sp),pointer:: FWRun(:)     ! [m/d]         ! 07 runoff
real(sp),pointer:: FWLchA(:)    ! [m/d]         ! 08 leach ex layer 1
real(sp),pointer:: FWLchB(:)    ! [m/d]         ! 09 leach ex layer 2 (deep drainage)
real(sp),pointer:: FWE(:)       ! [m/d]         ! 10 total evaporation
real(sp),pointer:: FWTra(:)     ! [m/d]'        ! 11 total transpiration
real(sp),pointer:: FWDis(:)     ! [m/d]'        ! 12 total point discharge (runoff+lch)
real(sp),pointer:: FWThrough(:) ! [m/d]'        ! 13 throughfall
real(sp),pointer:: FWwc(:)      ! [m/d]'        ! 14 flux from wet canopy
real(sp),pointer:: dWCanopydt(:)! [m/d]'        ! 15 dXX/dt: stored water in canopy
real(sp),pointer:: FCGPP(:)     ! [molC/m2/d]   ! 16 photsynthetic uptake (GPP)
real(sp),pointer:: FCGro(:)     ! [molC/m2/d]   ! 17 Growth C flux (NPP)
real(sp),pointer:: FCNEE(:)     ! [molC/m2/d]   ! 18 Net Ecosystem Exchange (NEE)
real(sp),pointer:: phiE(:)       ! [W/m2]         ! 19 daytime latent heat flux
real(sp),pointer:: phiH(:)       ! [W/m2]         ! 20 daytime sensible heat flux
real(sp),pointer:: FWPt(:)     ! [m/d]         ! 21 Preistley-Taylor (not calc in CABLE)
real(sp),pointer:: dCLeadt(:)     ! [molC/m2/d]         ! 22
real(sp),pointer:: dWcoldt(:)  ! [0,1]/d       ! 23 dXX/dt: change in total soil water 

! * Pointers for derived quantities (array DD)
real(sp),pointer:: DayltFrac(:) ! [-]           ! 01 daylight fraction
real(sp),pointer:: RhoA(:)      ! [molA/m3]     ! 02 air density
real(sp),pointer:: FracVExt(:)  ! [-]           ! 03 veg cover fraction (from external file)
real(sp),pointer:: rLAIExt(:)   ! [-]           ! 04 leaf area index    (from external file)
real(sp),pointer:: FracVCLea(:) ! [-]           ! 05 veg cover fraction (from CLea)
real(sp),pointer:: rLAICLea(:)  ! [-]           ! 06 leaf area index    (from CLea)
real(sp),pointer:: AllocLea(:)  ! [-]           ! 07 leaf allocation coeff
real(sp),pointer:: FCGroL(:)    ! [molC/m2/d]   ! 08 light-lim NPP
real(sp),pointer:: FCGroW(:)    ! [molC/m2/d]   ! 09 water-lim NPP
real(sp),pointer:: ImBalA(:)    ! [m/d]         ! 10 dWA/dt-(FWPrec-FWTraA-FWSoil-FWRun-FWLchA)
real(sp),pointer:: ImBalB(:)    ! [m/d]         ! 11 dWB/dt-(FWLchA-FWTraB-FWLchB)
real(sp),pointer:: ImBalCanopy(:)! [m/d]        ! 12 Canopy Imbalance
real(sp),pointer:: ImBalSoil(:)  ! [m/d]        ! 13   Soil Imbalance
real(sp),pointer:: ImBal(:)      ! [m/d]        ! 14 Total moisture imbalance
real(sp),pointer:: fws(:)        ! [-]          ! 15 limitation to stomatal conductance due to soil moisture deficit
real(sp),pointer:: Tabar(:)      ! [degC]       ! 16 mean air temperature (for CASA-CNP)
real(sp),pointer:: Tsoilbar_g(:) ! [degC]       ! 17 root-weighted mean soil temp, grass (for CASA-CNP)
real(sp),pointer:: Tsoilbar_w(:) ! [degC]       ! 18 root-weighted mean soil temp, wood (for CASA-CNP)
real(sp),pointer:: Sbar_g(:)     ! [0;1]        ! 19 root-weighted soil moisture relative to saturation, grass (for CASA-CNP)
real(sp),pointer:: Sbar_w(:)     ! [0;1]        ! 20 root-weighted soil moisture relative to saturation, wood (for CASA-CNP)
real(sp),pointer:: btran_g(:)    ! [0;1]        ! 21 water limitation factor for carbon allocation, grass (for CASA-CNP)
real(sp),pointer:: btran_w(:)    ! [0;1]        ! 22 water limitation factor for carbon allocation, wood (for CASA-CNP)
real(sp),pointer:: FCGPP_g(:)    ! [molC/m2]/d  ! 23 GPP , grass (per unit tile area) (for CASA-CNP)
real(sp),pointer:: FCGPP_w(:)    ! [molC/m2]/d  ! 24 GPP , wood (per unit tile area) (for CASA-CNP)
real(sp),pointer:: FCLeafR_g(:)  ! [molC/m2]/d  ! 25 day-time leaf respiration , grass (per unit tile area) (for CASA-CNP)
real(sp),pointer:: FCLeafR_w(:)  ! [molC/m2]/d  ! 26 GPP , day-time leaf respiration, grass (per unit tile area) (for CASA-CNP)
real(sp),pointer:: LAI_g(:)  ! [] ! 27
real(sp),pointer:: LAI_w(:)  ! [] ! 28
real(sp),pointer:: fws_g(:)  ! [] ! 29
real(sp),pointer:: fws_w(:)  ! [] ! 30
! * Pointers for predicted point-scale observables (array ZZP)
real(sp),pointer:: ZNDVI(:)      ! [degC]        ! 01 NDVI
real(sp),pointer:: ZphiRnet(:)   ! [MJ m-2 d-1]  ! 02 day-time Rnet
real(sp),pointer:: ZphiH(:)      ! [MJ m-2 d-1]  ! 03 day-time phiH
real(sp),pointer:: ZphiE(:)      ! [MJ m-2 d-1]  ! 04 day-time phiE
real(sp),pointer:: ZphiNEE(:)    ! [moles CO2 m-2 d-1]  ! 05 day-time phiNEE
real(sp),pointer:: ZphiNPP(:)    ! [moles CO2 m-2 d-1]  ! 06 day-time phiNPP
real(sp),pointer:: ZphiGPP(:)    ! [moles CO2 m-2 d-1]  ! 07 day-time phiGPP
real(sp),pointer:: Zsm0008(:)    ! [m3 m-3]       !08 volumetric soil moisture content 0-8 cm
real(sp),pointer:: Zsm0090(:)    ! [m3 m-3]       !09 volumetric soil moisture content 0-90 cm
real(sp),pointer:: Zsm0030(:)    ! [m3 m-3]       !10 volumetric soil moisture content 0-30 cm
real(sp),pointer:: Zsm3060(:)    ! [m3 m-3]       !11 volumetric soil moisture content 30-60 cm
real(sp),pointer:: Zsm6090(:)    ! [m3 m-3]       !12 volumetric soil moisture content 60-90 cm
! * Pointers for actual point-scale observations (array AAP)
real(sp),pointer:: ANDVI(:)      ! [degC]        ! 01 NDVI
real(sp),pointer:: AphiRnet(:)   ! [MJ m-2 d-1]  ! 02 day-time Rnet
real(sp),pointer:: AphiH(:)      ! [MJ m-2 d-1]  ! 03 day-time phiH
real(sp),pointer:: AphiE(:)      ! [MJ m-2 d-1]  ! 04 day-time phiE
real(sp),pointer:: AphiNEE(:)    ! [moles CO2 m-2 d-1]  ! 05 day-time phiNEE
real(sp),pointer:: AphiNPP(:)    ! [moles CO2 m-2 d-1]  ! 06 day-time phiNPP
real(sp),pointer:: AphiGPP(:)    ! [moles CO2 m-2 d-1]  ! 07 day-time phiGPP
real(sp),pointer:: Asm0008(:)    ! [m3 m-3]       !08 volumetric soil moisture content 0-15 cm
real(sp),pointer:: Asm0090(:)    ! [m3 m-3]       !09 volumetric soil moisture content 0-90 cm
real(sp),pointer:: Asm0030(:)    ! [m3 m-3]       !10 volumetric soil moisture content 0-30 cm
real(sp),pointer:: Asm3060(:)    ! [m3 m-3]       !11 volumetric soil moisture content 30-60 cm
real(sp),pointer:: Asm6090(:)    ! [m3 m-3]       !12 volumetric soil moisture content 60-90 cm
real(sp),pointer:: ALAIg(:)       ![] ! 13 grassy leaf area index (per unit tile area)
real(sp),pointer:: ALAIw(:)		![] ! 14 woody leaf area index (per unit tile area)
real(sp),pointer:: AfWoody(:)	![] ! 15 fraction cover woody
real(sp),pointer:: AscattVIS(:)	![] ! 16 monthly scattering coefficient (vis)
real(sp),pointer:: AscattNIR(:)	![] ! 17 monthly scattering coefficient (NIR)

! * Pointers for predicted point-scale observables (array ZZPh)
real(sp),pointer:: ZLST(:,:)      ! [degC]       ! 01 LST at overpass time
real(sp),pointer:: ZLSTTime(:,:)  ! [hrSol]      ! 02 overpass time for LST
real(sp),pointer:: ZLSTAngle(:,:) ! [deg]        ! 03 view angle for LST
real(sp),pointer:: DelTsTaZ(:,:)  ! [degC]       ! 04 DelTsTaA = ZLST - TempAt
real(sp),pointer:: TaZ(:,:)       ! [degC]       ! 05 TempAt
real(sp),pointer:: ZphiRneth(:,:) ! [Wm-2]       ! 06 phiRneth
real(sp),pointer:: ZphiHh(:,:)    ! [Wm-2]       ! 07 phiHh
real(sp),pointer:: ZphiEh(:,:)    ! [Wm-2]       ! 08 phiEh
real(sp),pointer:: ZphiNEEh(:,:)  ! [umol m-2 s-1]! 09 phiNEEh
real(sp),pointer:: ZphiGh(:,:)    ! [Wm-2]       ! 10 phiGh
real(sp),pointer:: ZTsoilh(:,:)    ! [deg C]       ! 11 Tsoilh
real(sp),pointer:: ZphiHhtime(:,:)! [h]          ! 12 phiHhtime
real(sp),pointer:: TrVeg(:,:)     ! [degC]       ! 13 Radiative temperature of veg
real(sp),pointer:: TrSoil(:,:)    ! [degC]       ! 14 Radiative temperature of soil
real(sp),pointer:: TrEff(:,:)     ! [degC]       ! 15 Effective Radiative temperature of surface (depends on view angle)
real(sp),pointer:: TSoil1h(:,:)     ! [degC]       ! 16 Soil T
real(sp),pointer:: TSoil2h(:,:)     ! [degC]       ! 17 Soil T
real(sp),pointer:: TSoil3h(:,:)     ! [degC]       ! 18 Soil T
real(sp),pointer:: TSoil4h(:,:)     ! [degC]       ! 19 Soil T
real(sp),pointer:: TSoil5h(:,:)     ! [degC]       ! 20 Soil T
real(sp),pointer:: TSoil6h(:,:)     ! [degC]       ! 21 Soil T
real(sp),pointer:: TSoil7h(:,:)     ! [degC]       ! 22 Soil T
real(sp),pointer:: TSoil8h(:,:)     ! [degC]       ! 23 Soil T
real(sp),pointer:: TSoil9h(:,:)     ! [degC]       ! 24 Soil T
real(sp),pointer:: TSoil10h(:,:)     ! [degC]       ! 25 Soil T
real(sp),pointer:: theta1h(:,:)     ! []       ! 26 Soil moisture
real(sp),pointer:: theta2h(:,:)     ! []       ! 27 Soil moisture
real(sp),pointer:: theta3h(:,:)     ! []       ! 28 Soil moisture
real(sp),pointer:: theta4h(:,:)     ! []       ! 29 Soil moisture
real(sp),pointer:: theta5h(:,:)     ! []       ! 30 Soil moisture
real(sp),pointer:: theta6h(:,:)     ! []       ! 31 Soil moisture
real(sp),pointer:: theta7h(:,:)     ! []       ! 32 Soil moisture
real(sp),pointer:: theta8h(:,:)     ! []       ! 33 Soil moisture
real(sp),pointer:: theta9h(:,:)     ! []       ! 34 Soil moisture
real(sp),pointer:: theta10h(:,:)     ! []       ! 35 Soil moisture

! * Pointers for actual point-scale observations (array AAPh)
real(sp),pointer:: ALST(:,:)      ! [degC]        ! 01 LST at overpass time
real(sp),pointer:: ALSTTime(:,:)  ! [hrSol]       ! 02 overpass time for LST
real(sp),pointer:: ALSTAngle(:,:) ! [deg]         ! 03 view angle for LST
real(sp),pointer:: DelTsTaA(:,:)  ! [degC]        ! 04 DelTsTaA = ALST - TempAt
real(sp),pointer:: TempAt(:,:)    ! [degC]        ! 05 TempAt
real(sp),pointer:: AphiRneth(:,:) ! [Wm-2]        ! 06 AphiRneth
real(sp),pointer:: AphiHh(:,:)    ! [Wm-2]        ! 07 AphiHh
real(sp),pointer:: AphiEh(:,:)    ! [Wm-2]        ! 08 AphiEh
real(sp),pointer:: AphiNEEh(:,:)  ! [umol m-2 s-1]! 09 AphiNEEh
real(sp),pointer:: AphiGh(:,:)    ! [Wm-2]        ! 10 AphiGh
real(sp),pointer:: ATsoilh(:,:)    ! [degC]        ! 11 AphiTsoilh
real(sp),pointer:: AphiHhtime(:,:)! [h]           ! 12 AphiHhtime
! * Pointers for predicted catchment-scale observables (array ZZC)
real(sp),pointer:: ZDisCD(:)    ! [m/d]         ! 01 catchment daily outflow
real(sp),pointer:: ZDisCM(:)    ! [m/mth]       ! 02 catchment monthly outflow
real(sp),pointer:: ZRunCD(:)    ! [m/d]         ! 03 catchment daily runoff
real(sp),pointer:: ZLchCD(:)    ! [m/d]         ! 04 catchment daily leaching
! * Pointers for actual catchment-scale observations (array AAC)
real(sp),pointer:: ADisCD(:)    ! [m/d]         ! 01 catchment daily outflow
real(sp),pointer:: ADisCM(:)    ! [m/mth]       ! 02 catchment monthly outflow
real(sp),pointer:: ARunCD(:)    ! [m/d]         ! 03 catchment daily runoff
real(sp),pointer:: ALchCD(:)    ! [m/d]         ! 04 catchment daily leaching
! * Pointers for met forcing variables (array MM)
real(sp),pointer:: SolarMJ(:)   ! [MJ/d]        ! 01 incident solar radn
real(sp),pointer:: Precip(:)    ! [m/d]         ! 02 precipitation
real(sp),pointer:: TempMax(:)   ! [degC]        ! 03 24hr-max air temp
real(sp),pointer:: TempMin(:)   ! [degC]        ! 04 24hr-min air temp
real(sp),pointer:: vph09(:)     ! [mb]          ! 05 24hr-min air temp
real(sp),pointer:: vph15(:)     ! [mb]          ! 06 24hr-min air temp
! * Pointers for hourly met  forcing variables (array hMM)
real(sp),pointer:: hFsd(:,:)       ! [W m-2]       ! 01 incident solar radn
real(sp),pointer:: hFld(:,:)       ! [W m-2]       ! 02 down-welling lowg-wave
real(sp),pointer:: hPrecip(:,:)    ! [mm/dt]       ! 03 precip
real(sp),pointer:: hUa(:,:)        ! [ms-1]        ! 04 windspeed at ref height
real(sp),pointer:: hTc(:,:)        ! [degC]        ! 05 T at ref height
real(sp),pointer:: hqv(:,:)        ! [kg/kg]       ! 06 specific humidity at ref height
real(sp),pointer:: hpmb(:,:)       ! [mb]          ! 07 atmospheric pressure
real(sp),pointer:: hcoszen(:,:)       ! [mb]       ! 08 cos (solar zenith angle)
! * Pointers for remote sensing driver variables (array RR)
real(sp),pointer:: LAIg(:)   ! []        ! 01 grassy LAI 
real(sp),pointer:: LAIw(:)    ! []         ! 02 woody LAI
real(sp),pointer:: fWoody(:)   ! []        ! 03 fraction woody cover
real(sp),pointer:: scattVIS(:)   ! []        ! 04 leaf scattering coefficient (vis)
real(sp),pointer:: scattNIR(:)     ! []          ! 05 leaf scattering coeffcient (NIR)
! * Pointers for spatially uniform parameters (array UU)
real(sp),pointer:: CoeffPT      ! [-]           ! 01 Priestley-Taylor coeff
real(sp),pointer:: CoeffBeer    ! [-]           ! 02 Beer Law extinction coeff
real(sp),pointer:: CLea0        ! [-]           ! 03 CLea0 = RhoCLeaf*LeafThick
real(sp),pointer:: RateCLea     ! [1/d]         ! 04 rate constant for CLea decay
real(sp),pointer:: AllocLg      ! []         ! 05 allocation of C to leaves (grassy)
real(sp),pointer:: AllocLw     ! []         ! 06 allocation of C to leaves (woody)
real(sp),pointer:: HyConSat1    ! [m/d]         ! 07 sat hydraulic conductivity (1)
real(sp),pointer:: HyConSat2    ! [m/d]         ! 08 sat hydraulic conductivity (2)
real(sp),pointer:: PwrFWSoil    ! [-]           ! 09 FWSoil ~ WRel**(PwrFWSoil+1)
real(sp),pointer:: PwrFWLch     ! [-]           ! 10 FWLch  ~ WRel**(PwrFWLch+1)
real(sp),pointer:: alfaQ        ! [molC/molQ]   ! 11 Light Use Efficiency
real(sp),pointer:: alfaWpri     ! [molC/W]      ! 12 prior Water Use Eff   (use if > 0)
real(sp),pointer:: alfaWmul     ! [-]           ! 13 mult for WUE from def (use if > 0)
real(sp),pointer:: CO2A         ! [molC/A]      ! 14 air [CO2]
real(sp),pointer:: WRelA0       ! [-]           ! 15 scale for AllocLea(WRel)
real(sp),pointer:: ZSoil1Mult   ! [-]           ! 16 multiplier for ZSoil1
real(sp),pointer:: ZSoil2Mult   ! [-]           ! 17 multiplier for ZSoil2
real(sp),pointer:: rLAImax      ! [-]           ! 18 maximum rLAI
real(sp),pointer:: Gaero        ! [m/s]         ! 19 reference aerodynamic conductance
real(sp),pointer:: TimeTxFrac   ! [-]           ! 20 Tmax time as frac of (TDawn,TDusk)
real(sp),pointer:: cN0          ! [-]           ! 21 ZNDVI = Sum[cNi*FracV^i]
real(sp),pointer:: cN1          ! [-]           ! 22 ZNDVI = Sum[cNi*FracV^i]
real(sp),pointer:: cN2          ! [-]           ! 23 ZNDVI = Sum[cNi*FracV^i]
real(sp),pointer:: TZRunDef     ! [d]           ! 24 default catchment runoff timescale (use if >0)
real(sp),pointer:: TZLchDef     ! [d]           ! 25 default catchment leach  timescale (use if >0)
real(sp),pointer:: CoeffPAR     ! [-]           ! 26 FracV = CoeffPAR*FAPAR
real(sp),pointer:: RatioJV      ! [-]           ! 27 RatioJV=Jmax/Vcmax
real(sp),pointer:: za           ! [m]           ! 28 reference height
real(sp),pointer:: ratecp1      ! [1/year]		! 29 rate constant: plant carbon pool 1
real(sp),pointer:: ratecp2      ! [1/year]		! 30 rate constant: plant carbon pool 2
real(sp),pointer:: ratecp3      ! [1/year]		! 31 rate constant: plant carbon pool 3
real(sp),pointer:: ratecs1      ! [1/year]		! 32 rate constant: soil carbon pool 1
real(sp),pointer:: ratecs2      ! [1/year]		! 33 rate constant: soil carbon pool 2
real(sp),pointer:: zeta      ! [-]		    ! 34 macropore flow parameter
real(sp),pointer:: fsatmax      ! [-]		! 35 Multiplier for litter depth
real(sp),pointer:: B1Mult      ! [-]		    ! 36 Multiplier for b1
real(sp),pointer:: Psie1Mult      ! [-]		! 37 Multiplier for psie1
real(sp),pointer:: Ksat1Mult      ! [-]		! 38 Multiplier for Ksat
real(sp),pointer:: B2Mult      ! [-]		    ! 39 Multiplier for b1
real(sp),pointer:: Psie2Mult      ! [-]		! 40 Multiplier for psie1
real(sp),pointer:: Ksat2Mult      ! [-]		! 41 Multiplier for Ksat
real(sp),pointer:: dleaf_g	    ![m]		! 42 leaf length (grass)
real(sp),pointer:: vcmax_g	    ![mol m-2 s-1] ! 43 maximum RuBP carboxylation rate top leaf (grass)
real(sp),pointer:: hc_g	    ![m] ! 44 canopy ht (grass)
real(sp),pointer:: xfang_g	    ![-] ! 45 leaf angle dist parameter (grass)
real(sp),pointer:: rp20_g	    ![m] ! 46 relative plant respiration coefficient (grass)
real(sp),pointer:: vbeta_g	    ![m] ! 47 stomatal sensitivity to soil water (grass)
real(sp),pointer:: rootbeta_g	    ![m] ! 48 parameter to describe root density distribution (grass)
real(sp),pointer:: F10_g	    ![m] ! 49 fraction of roots in top 10 cm (grass)
real(sp),pointer:: ZR_g	    ![m] ! 50 maximum rooting depth (grass)
real(sp),pointer:: loggamma_g	    ![-] ! 51 parameter in root efficiency function (grass)
real(sp),pointer:: dleaf_w	    ![m]		! 52 leaf length (woody)
real(sp),pointer:: vcmax_w	    ![mol m-2 s-1] ! 53 maximum RuBP carboxylation rate top leaf (woody)
real(sp),pointer:: hc_w	    ![m] ! 54 canopy ht (woody)
real(sp),pointer:: xfang_w	    ![-] ! 55 leaf angle dist parameter (woody)
real(sp),pointer:: rp20_w	    ![m] ! 56 relative plant respiration coefficient (woody)
real(sp),pointer:: vbeta_w	    ![m] ! 57 stomatal sensitivity to soil water (woody)
real(sp),pointer:: rootbeta_w	    ![m] ! 58 parameter to describe root density distribution (woody)
real(sp),pointer:: F10_w	    ![m] ! 59 fraction of roots in top 10 cm (woody)
real(sp),pointer:: ZR_w	    ![m] ! 60 maximum rooting depth (woody)
real(sp),pointer:: loggamma_w	    ![-] ! 61 parameter in root efficiency function (woody)
real(sp),pointer:: a1	    ![-] ! 62 parameter in stomatal conductance function
real(sp),pointer:: ds0	    ![-] ! 63 sensitivity of stomatal conduxtance to VPD


! * Pointers for spatially variable parameters (array VV)
real(sp),pointer:: CatchMap(:)  ! [-]           ! 01 CRCCH catchment ID map
real(sp),pointer:: Albedo(:)    ! [-]           ! 02 albedo
real(sp),pointer:: LatDeg(:)    ! [degN]        ! 03 latitude
real(sp),pointer:: LongDeg(:)    ! [degE]        ! 04 longitude
real(sp),pointer:: Altitude(:)  ! [m]           ! 05 altitude = elevation
real(sp),pointer:: FAPAR(:,:)   ! [-]           ! 06-17 observed FAPAR: Jan-Dec
real(sp),pointer:: canst1(:)    ! [mm/LAI]      ! 18 max intercepted water by canopy (mm/LAI)
real(sp),pointer:: dleaf(:)     ! [m]           ! 19 leaf length (m)
real(sp),pointer:: vcmax(:)     ! [mol m-2 s-1] ! 20 maximum RuBP carboxylation rate top leaf (mol/m2/s)
real(sp),pointer:: hc(:)        ! [m]           ! 21 canopy height (m)
real(sp),pointer:: xfang(:)     ! [-]           ! 22 leaf angle dist
real(sp),pointer:: rp20(:)      ! [-]           ! 23 relative plant respiration coefficient
real(sp),pointer:: rpcoef(:)    ![1/degC]       ! 24 temperature coef nonleaf plant respiration 
real(sp),pointer:: rs20(:)     ! [-]           ! 25 relative soil respiration at 20C (tuning parameter)
real(sp),pointer:: shelrb(:)    ! [-]           ! 26 sheltering factor (dimensionless)
real(sp),pointer:: frac4(:)     ! [-]           ! 27 fraction of c4 plant
real(sp),pointer:: tminvj(:)    ! [degC]        ! 28 min temperature of the start of photosynthesis
real(sp),pointer:: tmaxvj(:)    ! [degC]        ! 29 min temperature of the start of photosynthesis
real(sp),pointer:: vbeta(:)     ! [-]           ! 30 stomatal sensitivity to soil water
real(sp),pointer:: albsoil(:)   ! [-]           ! 31 soil albedo
real(sp),pointer:: swilt1(:)    ! [-]           ! 32 vol H2O @ wilting (A horizon)
real(sp),pointer:: sfc1(:)      ! [-]           ! 33 vol H2O @ field capacity (A horizon)
real(sp),pointer:: WVolSat1(:)  ! [-]           ! 34 vol H2O @ saturation (A horizon)
real(sp),pointer:: bch1(:)      ! [-]           ! 35 parameter b in Campbell equation (A horizon)
real(sp),pointer:: hyds1(:)     ! [m/s]			! 36 hydraulic conductivity @ saturation (A horizon)
real(sp),pointer:: sucs1(:)     ! [m]			! 37 suction at saturation (A horizon)
real(sp),pointer:: rhosoil1(:)  ! [kg/m3]		! 38 soil density (A horizon)
real(sp),pointer:: css1(:)      ! [kJ/kg/K]		! 39 thermal conductivity of dry soil (A horizon)
real(sp),pointer:: clay1(:)    ! [-]			! 40 fraction clay (A horizon)
real(sp),pointer:: silt1(:)    ! [%]			! 41 fraction silt(A horizon)
real(sp),pointer:: ZSoil1(:)    ! [m]			! 42 water-holding soil depth (A horizon)
real(sp),pointer:: swilt2(:)    ! [-]           ! 43 vol H2O @ wilting (B horizon)
real(sp),pointer:: sfc2(:)      ! [-]           ! 44 vol H2O @ field capacity (B horizon)
real(sp),pointer:: WVolSat2(:)  ! [-]           ! 45 vol H2O @ saturation (B horizon)
real(sp),pointer:: bch2(:)      ! [-]           ! 46 parameter b in Campbell equation (B horizon)
real(sp),pointer:: hyds2(:)     ! [m/s]         ! 47 hydraulic conductivity @ saturation (B horizon)
real(sp),pointer:: sucs2(:)     ! [m]			! 48 suction at saturation (B horizon)
real(sp),pointer:: rhosoil2(:)  ! [kg/m3]		! 49 soil density (B horizon)
real(sp),pointer:: css2(:)      ! [kJ/kg/K]		! 50 thermal conductivity of dry soil (B horizon)
real(sp),pointer:: ZSoil2(:)    ! [m]			! 51 water-holding soil depth (layer 1) (B horizon)
real(sp),pointer:: clay2(:)     ! [m]			! 52 fraction clay (B horizon)
real(sp),pointer:: silt2(:)     ! [m]			! 53 fraction silt (B horizon)
real(sp),pointer:: rootbeta(:)	    ![-]		! 54 parameter to describe root density distribution
real(sp),pointer:: extkn(:)		![-]		    ! 55 extinction coefficient for nitrogen with canopy depth
real(sp),pointer:: vegcf(:)		![-]		    ! 56  biome-specific soil respiration rate
real(sp),pointer:: iveg(:)		![-]		    ! 57 IGBP vegetation type
real(sp),pointer:: LAI(:,:)   ! [-]             ! 58-69 LAI: Jan-Dec
real(sp),pointer:: clittequil(:)    ! [tC/ha]        ! 70 litter pool (from BIOS)
real(sp),pointer:: albsoilVIS(:) ! [-]  ! 71 soil albedo (vis)
real(sp),pointer:: albsoilNIR(:) ! [-]  ! 72 soil albedo (NIR)

CONTAINS

!*******************************************************************************

SUBROUTINE PointAll (TargetTTime, TargetXX, TargetFF, TargetDD,     &
                     TargetZZP, TargetAAP,TargetZZPh, TargetAAPh,   &
					 TargetZZC, TargetAAC,    &
                     TargetMM,TargetRR,TargetUU, TargetVV,TargethMM)
!-------------------------------------------------------------------------------
! Assign all pointers from generic to specific array names
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Input variables
real(sp),intent(in),target,optional:: TargetTTime(:)
real(sp),intent(in),target,optional:: TargetXX(:,:)
real(sp),intent(in),target,optional:: TargetFF(:,:)
real(sp),intent(in),target,optional:: TargetDD(:,:)
real(sp),intent(in),target,optional:: TargetZZP(:,:)
real(sp),intent(in),target,optional:: TargetAAP(:,:)
real(sp),intent(in),target,optional:: TargetZZPh(:,:,:)
real(sp),intent(in),target,optional:: TargetAAPh(:,:,:)
real(sp),intent(in),target,optional:: TargetZZC(:,:)
real(sp),intent(in),target,optional:: TargetAAC(:,:)
real(sp),intent(in),target,optional:: TargetMM(:,:)
real(sp),intent(in),target,optional:: TargethMM(:,:,:)
real(sp),intent(in),target,optional:: TargetRR(:,:)
real(sp),intent(in),target,optional:: TargetUU(:)
real(sp),intent(in),target,optional:: TargetVV(:,:)
!-------------------------------------------------------------------------------
! * Pointers for time variables (array TTime)
if (present(TargetTTime)) then
  if (size(TargetTTime) /= 5) stop "PointAll: wrong TTime dimension"
  Time1y     => TargetTTime(01)
  TTDay      => TargetTTime(02)
  TTMonth    => TargetTTime(03)
  TTYear     => TargetTTime(04)
  TTEndMth   => TargetTTime(05)
end if
! * Pointers for stores (array XX)
if (present(TargetXX)) then
  if (size(TargetXX,2) /= 40) stop "PointAll: wrong nxx in XX"
  WRelA      => TargetXX(:,01)
  WRelB      => TargetXX(:,02)
  CLea       => TargetXX(:,03)
  WRel1      => TargetXX(:,04)
  WRel2      => TargetXX(:,05)
  WRel3      => TargetXX(:,06)
  WRel4      => TargetXX(:,07)
  WRel5      => TargetXX(:,08)
  WRel6      => TargetXX(:,09)
  WRel7      => TargetXX(:,10)
  WRel8      => TargetXX(:,11)
  WRel9      => TargetXX(:,12)
  WRel10      => TargetXX(:,13)
  Tsoil1     => TargetXX(:,14)
  Tsoil2     => TargetXX(:,15)
  Tsoil3     => TargetXX(:,16)
  Tsoil4     => TargetXX(:,17)
  Tsoil5     => TargetXX(:,18)
  Tsoil6     => TargetXX(:,19)
  Tsoil7     => TargetXX(:,20)
  Tsoil8     => TargetXX(:,21)
  Tsoil9     => TargetXX(:,22)
  Tsoil10     => TargetXX(:,23)
  cplant1    => TargetXX(:,24)
  cplant2    => TargetXX(:,25)
  cplant3    => TargetXX(:,26)
  csoil1     => TargetXX(:,27)
  csoil2     => TargetXX(:,28)
  cansto     => TargetXX(:,29)
  theta1	 => TargetXX(:,30)
  theta2	 => TargetXX(:,31)
  theta3	 => TargetXX(:,32)
  theta4	 => TargetXX(:,33)
  theta5	 => TargetXX(:,34)
  theta6	 => TargetXX(:,35)
  theta7	 => TargetXX(:,36)
  theta8	 => TargetXX(:,37)
  theta9	 => TargetXX(:,38)
  theta10	 => TargetXX(:,39)
  wcol		 => TargetXX(:,40)
end if
! * Pointers for fluxes and dXXdt (array FF)
if (present(TargetFF)) then
  if (size(TargetFF,2) /= 23) stop "PointAll: wrong nff"
  dWColAdt   => TargetFF(:,01)
  dWColBdt   => TargetFF(:,02)
  FWPrec     => TargetFF(:,03)
  FWTraA     => TargetFF(:,04)
  FWTraB    => TargetFF(:,05)
  FWSoil     => TargetFF(:,06)
  FWRun      => TargetFF(:,07)
  FWLchA     => TargetFF(:,08)
  FWLchB     => TargetFF(:,09)
  FWE        => TargetFF(:,10)
  FWTra      => TargetFF(:,11)
  FWDis      => TargetFF(:,12)
  FWThrough  => TargetFF(:,13)
  FWwc       => TargetFF(:,14)
  dWcanopydt => TargetFF(:,15)
  FCGPP      => TargetFF(:,16)
  FCGro      => TargetFF(:,17)
  FCNEE      => TargetFF(:,18)
  phiE       => TargetFF(:,19)
  phiH       => TargetFF(:,20)
  FWPT       => TargetFF(:,21)
  dCLeadt       => TargetFF(:,22)
  dWColdt   => TargetFF(:,23)
end if
if (present(TargetDD)) then
  if (size(TargetDD,2) /= 30) stop "PointAll: wrong ndd"
! * Pointers for derived quantities (array DD)
  DayltFrac   => TargetDD(:,01)
  RhoA        => TargetDD(:,02)
  FracVExt    => TargetDD(:,03)
  rLAIExt     => TargetDD(:,04)
  FracVCLea   => TargetDD(:,05)
  rLAICLea    => TargetDD(:,06)
  AllocLea    => TargetDD(:,07)
  FCGroL      => TargetDD(:,08)
  FCGroW      => TargetDD(:,09)
  ImBalA      => TargetDD(:,10)
  ImBalB      => TargetDD(:,11)
  ImBalCanopy => TargetDD(:,12)
  ImBalSoil   => TargetDD(:,13)
  ImBal       => TargetDD(:,14)
  fws      => TargetDD(:,15)
  Tabar     => TargetDD(:,16)   
  Tsoilbar_g    => TargetDD(:,17)
  Tsoilbar_w    => TargetDD(:,18)
  Sbar_g      => TargetDD(:,19)
  Sbar_w      => TargetDD(:,20) 
  btran_g     => TargetDD(:,21) 
  btran_w     => TargetDD(:,22)
  FCGPP_g     => TargetDD(:,23) 
  FCGPP_w     => TargetDD(:,24)  
  FCLeafR_g    => TargetDD(:,25)
  FCLeafR_w    => TargetDD(:,26)
  LAI_g => TargetDD(:,27)
  LAI_w => TargetDD(:,28)
  fws_g => TargetDD(:,29)
  fws_w => TargetDD(:,30)
end if
! * Pointers for predicted point observables (array ZZP)
if (present(TargetZZP)) then
  if (size(TargetZZP,2) /= 12) stop "PointAll: wrong nzzP"
  ZNDVI       => TargetZZP(:,01)
  ZphiRnet    => TargetZZP(:,02)
  ZphiH       => TargetZZP(:,03)
  ZphiE       => TargetZZP(:,04)
  ZphiNEE     => TargetZZP(:,05)
  ZphiNPP     => TargetZZP(:,06)
  ZphiGPP     => TargetZZP(:,07)
  Zsm0008    => TargetZZP(:,08)
  Zsm0090     => TargetZZP(:,09)
  Zsm0030     => TargetZZP(:,10)
  Zsm3060     => TargetZZP(:,11)
  Zsm6090     => TargetZZP(:,12)
end if
  ! * Pointers for predicted point observables (array AAP)
if (present(TargetAAP)) then
  if (size(TargetAAP,2) /= 12) stop "PointAll: wrong naaP"
  ANDVI       => TargetAAP(:,01)
  AphiRnet    => TargetAAP(:,02)
  AphiH       => TargetAAP(:,03)
  AphiE       => TargetAAP(:,04)
  AphiNEE     => TargetAAP(:,05)
  AphiNPP     => TargetAAP(:,06)
  AphiGPP     => TargetAAP(:,07)
  Asm0008     => TargetAAP(:,08)
  Asm0090     => TargetAAP(:,09)
  Asm0030     => TargetAAP(:,10)
  Asm3060     => TargetAAP(:,11)
  Asm6090     => TargetAAP(:,12)
  ALAIg     => TargetAAP(:,13)
  ALAIw     => TargetAAP(:,14)
  AfWoody     => TargetAAP(:,15)
  AscattVIS     => TargetAAP(:,16)
  AscattNIR     => TargetAAP(:,17)

end if
! * Pointers for predicted point observables (array ZZPh)
if (present(TargetZZPh)) then
  if (size(TargetZZPh,2) /= 35) stop "PointAll: wrong nzzPh"
  ZLST        => TargetZZPh(:,01,:)
  ZLSTTime    => TargetZZPh(:,02,:)
  ZLSTAngle   => TargetZZPh(:,03,:)
  DelTsTaZ    => TargetZZPh(:,04,:)
  TaZ         => TargetZZPh(:,05,:)
  ZphiRneth   => TargetZZPh(:,06,:)
  ZphiHh      => TargetZZPh(:,07,:)
  ZphiEh      => TargetZZPh(:,08,:)
  ZphiNEEh    => TargetZZPh(:,09,:)
  ZphiGh      => TargetZZPh(:,10,:)
  ZTsoilh      => TargetZZPh(:,11,:)
  ZphiHhTime  => TargetZZPh(:,12,:)
  TrVeg       => TargetZZPh(:,13,:)
  TrSoil      => TargetZZPh(:,14,:)
  TrEff       => TargetZZPh(:,15,:)
  TSoil1h       => TargetZZPh(:,16,:)
  TSoil2h       => TargetZZPh(:,17,:)
  TSoil3h       => TargetZZPh(:,18,:)
  TSoil4h       => TargetZZPh(:,19,:)
  TSoil5h       => TargetZZPh(:,20,:)
  TSoil6h       => TargetZZPh(:,21,:)
  TSoil7h       => TargetZZPh(:,22,:)
  TSoil8h       => TargetZZPh(:,23,:)
  TSoil9h       => TargetZZPh(:,24,:)
  TSoil10h       => TargetZZPh(:,25,:)
  theta1h       => TargetZZPh(:,26,:)
  theta2h       => TargetZZPh(:,27,:)
  theta3h       => TargetZZPh(:,28,:)
  theta4h       => TargetZZPh(:,29,:)
  theta5h       => TargetZZPh(:,30,:)
  theta6h       => TargetZZPh(:,31,:)
  theta7h       => TargetZZPh(:,32,:)
  theta8h       => TargetZZPh(:,33,:)
  theta9h       => TargetZZPh(:,34,:)
  theta10h       => TargetZZPh(:,35,:)
end if
! * Pointers for actual point observations (array AAPh)
if (present(TargetAAPh)) then
  if (size(TargetAAPh,2) /= 12) stop "PointAll: wrong naaPh"
  ALST       => TargetAAPh(:,01,:)
  ALSTTime   => TargetAAPh(:,02,:)
  ALSTAngle  => TargetAAPh(:,03,:)
  DelTsTaA   => TargetAAPh(:,04,:)
  TempAt     => TargetAAPh(:,05,:)
  AphiRneth  => TargetAAPh(:,06,:)
  AphiHh     => TargetAAPh(:,07,:)
  AphiEh     => TargetAAPh(:,08,:)
  AphiNEEh   => TargetAAPh(:,09,:)
  AphiGh     => TargetAAPh(:,10,:)
  ATsoilh     => TargetAAPh(:,11,:)
  AphiHhTime => TargetAAPh(:,12,:)
end if
! * Pointers for predicted catchment observables (array ZZC)
if (present(TargetZZC)) then
  if (size(TargetZZC,2) /= 4) stop "PointAll: wrong nzzC"
  ZDisCD     => TargetZZC(:,01)
  ZDisCM     => TargetZZC(:,02)  
  ZRunCD     => TargetZZC(:,03)
  ZLchCD     => TargetZZC(:,04)
end if
! * Pointers for actual catchment observables (array AAC)
if (present(TargetAAC)) then
  if (size(TargetAAC,2) /= 4) stop "PointAll: wrong naaC"
  ADisCD     => TargetAAC(:,01)
  ADisCM     => TargetAAC(:,02)  
  ARunCD     => TargetAAC(:,03)
  ALchCD     => TargetAAC(:,04)
end if
! * Pointers for met forcing variables (array MM)
if (present(TargetMM)) then
  if (size(TargetMM,2) /= 6) stop "PointAll: wrong nmm"
  SolarMJ    => TargetMM(:,01)
  Precip     => TargetMM(:,02)
  TempMax    => TargetMM(:,03)
  TempMin    => TargetMM(:,04)
  vph09      => TargetMM(:,05)
  vph15      => TargetMM(:,06)
end if
! * Pointers for met forcing variables (array hMM)
if (present(TargethMM)) then
  if (size(TargethMM,2) /= 8) stop "PointAll: wrong nhMM"
  hFsd	    => TargethMM(:,01,:)
  hFld      => TargethMM(:,02,:)
  hPrecip	=> TargethMM(:,03,:)
  hUa       => TargethMM(:,04,:)
  hTc	    => TargethMM(:,05,:)
  hqv	    => TargethMM(:,06,:)
  hpmb	    => TargethMM(:,07,:)
  hcoszen	=> TargethMM(:,08,:)
end if
! * Pointers for remote seinsing forcing variables (array RR)
if (present(TargetRR)) then
  if (size(TargetRR,2) /= 5) stop "PointAll: wrong nRR"
  LAIg		=> TargetRR(:,01)
  LAIw		=> TargetRR(:,02)
  fWoody    => TargetRR(:,03)
  scattVIS  => TargetRR(:,04)
  scattNIR  => TargetRR(:,05)
end if
! * Pointers for spatially uniform parameters (array UU)
if (present(TargetUU)) then
  if (size(TargetUU) /= 63) stop "PointAll: wrong nuu"
  CoeffPT    => TargetUU(01)
  CoeffBeer  => TargetUU(02)
  CLea0      => TargetUU(03)
  RateCLea   => TargetUU(04)
  AllocLg    => TargetUU(05)
  AllocLw    => TargetUU(06)
  HyConSat1  => TargetUU(07)
  HyConSat2  => TargetUU(08)
  PwrFWSoil  => TargetUU(09)
  PwrFWLch   => TargetUU(10)
  alfaQ      => TargetUU(11)
  alfaWpri   => TargetUU(12)
  alfaWmul   => TargetUU(13)
  CO2A       => TargetUU(14)
  WRelA0     => TargetUU(15)
  ZSoil1Mult => TargetUU(16)
  ZSoil2Mult => TargetUU(17)
  rLAImax    => TargetUU(18)
  Gaero      => TargetUU(19)
  TimeTxFrac => TargetUU(20)
  cN0        => TargetUU(21)
  cN1        => TargetUU(22)
  cN2        => TargetUU(23)
  TZRunDef   => TargetUU(24)
  TZLchDef   => TargetUU(25)
  CoeffPAR   => TargetUU(26)
  RatioJV    => TargetUU(27)
  za         => TargetUU(28)
  ratecp1    => TargetUU(29)
  ratecp2    => TargetUU(30)
  ratecp3    => TargetUU(31)
  ratecs1    => TargetUU(32)
  ratecs2    => TargetUU(33)
  zeta    => TargetUU(34)
  fsatmax    => TargetUU(35)
  B1Mult    => TargetUU(36)
  Psie1Mult    => TargetUU(37)
  Ksat1Mult    => TargetUU(38)
  B2Mult    => TargetUU(39)
  Psie2Mult    => TargetUU(40)
  Ksat2Mult => TargetUU(41)
  dleaf_g    => TargetUU(42)
  vcmax_g    => TargetUU(43)
  hc_g       => TargetUU(44)
  xfang_g    => TargetUU(45)
  rp20_g     => TargetUU(46)
  vbeta_g    => TargetUU(47)
  rootbeta_g => TargetUU(48)
  F10_g      => TargetUU(49)
  ZR_g       => TargetUU(50)
  loggamma_g    => TargetUU(51)
  dleaf_w    => TargetUU(52)
  vcmax_w    => TargetUU(53)
  hc_w       => TargetUU(54)
  xfang_w    => TargetUU(55)
  rp20_w     => TargetUU(56)
  vbeta_w    => TargetUU(56)
  rootbeta_w => TargetUU(58)
  F10_w      => TargetUU(59)
  ZR_w       => TargetUU(60)
  loggamma_w    => TargetUU(61)
  a1		    => TargetUU(62)
  ds0		    => TargetUU(63)
end if
! * Pointers for spatially variable parameters (array VV)
if (present(TargetVV)) then
  if (size(TargetVV,2) /= 72) stop "PointAll: wrong nvv"
  CatchMap   => TargetVV(:,01)
  Albedo     => TargetVV(:,02)
  LatDeg     => TargetVV(:,03)
  LongDeg     => TargetVV(:,04)
  Altitude   => TargetVV(:,05)
  FAPAR      => TargetVV(:,06:17)
  canst1     => TargetVV(:,18)
  dleaf		 => TargetVV(:,19)
  vcmax		 => TargetVV(:,20)
  hc		 => TargetVV(:,21)
  xfang		 => TargetVV(:,22)
  rp20		 => TargetVV(:,23)
  rpcoef	 => TargetVV(:,24)
  rs20		 => TargetVV(:,25)
  shelrb	 => TargetVV(:,26)
  frac4		 => TargetVV(:,27)
  tminvj	 => TargetVV(:,28)
  tmaxvj	 => TargetVV(:,29)
  vbeta		 => TargetVV(:,30)
  albsoil	 => TargetVV(:,31)
  swilt1	 => TargetVV(:,32)
  sfc1		 => TargetVV(:,33)
  WVolSat1   => TargetVV(:,34)
  bch1		 => TargetVV(:,35)
  hyds1		 => TargetVV(:,36)
  sucs1		 => TargetVV(:,37)
  rhosoil1	 => TargetVV(:,38)
  css1		 => TargetVV(:,39)
  ZSoil1     => TargetVV(:,40)
  clay1      => TargetVV(:,41)
  silt1      => TargetVV(:,42)
  swilt2	 => TargetVV(:,43)
  sfc2		 => TargetVV(:,44)
  WVolSat2   => TargetVV(:,45)
  bch2		 => TargetVV(:,46)
  hyds2		 => TargetVV(:,47)
  sucs2		 => TargetVV(:,48)
  rhosoil2	 => TargetVV(:,49)
  css2		 => TargetVV(:,50)
  ZSoil2     => TargetVV(:,51)
  clay2      => TargetVV(:,52)
  silt2      => TargetVV(:,53)
  rootbeta	 => TargetVV(:,54)
  extkn      => TargetVV(:,55)
  vegcf      => TargetVV(:,56)
  iveg       => TargetVV(:,57)
  LAI         => TargetVV(:,58:69)
  clittequil     => TargetVV(:,70)
  albsoilVIS     => TargetVV(:,71)
  albsoilNIR     => TargetVV(:,72)
end if

END SUBROUTINE PointAll

!*******************************************************************************

END MODULE PointerModule


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

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

!*******************************************************************************
!*******************************************************************************

MODULE SubDiurnalMetModuleOld
CONTAINS
SUBROUTINE SubDiurnalMet (ntime, SolarMJDay, YearDay, LatDeg,               &
           PrecipDay, WindDay, TempMinDay, TempMaxDay, VapPmbDay, PmbDay,   &
           DecDeg, DayLength, SolarNorm, SolarFracObs, SolarTotChk,         &
           WindAvgChk, RelHumTempMin,                                       &
           PhiSd, PhiLd, Precip, Wind, Temp, VapPmb, Pmb,coszen)
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
!   * Fit a cosine wave passing through TempMaxDay and TempMinDay, with peak at
!     TimeMaxTemp (to be set, eg 16.00 hours) and trough 12 hours earlier.
! * Water Vapour Pressure, Air Pressure:
!   * Hold constant. Return rel humidity at TempMinDay as a diagnostic on obs. 
!
! HISTORY:
! * 04-sep-2003 (MRR): Written and tested in Program SubDiurnalWeather
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
USE constants
implicit none
! Parameters
!real(sp),parameter:: Pi         = 3.14159265        ! Pi
!real(sp),parameter:: SecDay     = 86400.0           ! Seconds/day
real(sp),parameter:: SolarConst = 1370*SecDay/1e6   ! Solar constant [MJ/m2/day]
real(sp),parameter:: epsilon  = 0.736
! Global variables
integer(i4b),intent(in) :: ntime    ! number of subdiurnal steps in 24 hr
real(sp),intent(in) :: SolarMJday   ! 24hr-total incident solar radiation [MJ/day]
integer(i4b),intent(in) :: YearDay  ! Day of the year (1-366)
real(sp),intent(in) :: LatDeg       ! Latitude [degs, neg for SH]
real(sp),intent(in) :: PrecipDay    ! 24hr-total precipitation [m/day] 
real(sp),intent(in) :: WindDay      ! 24hr-av wind [m/s] 
real(sp),intent(in) :: TempMinDay   ! Daily minimum air temp [degC]
real(sp),intent(in) :: TempMaxDay   ! Daily maximum air temp [degC]
real(sp),intent(in) :: VapPmbDay    ! 24hr-av water vapour pressure [mb]
real(sp),intent(in) :: PmbDay       ! 24hr-av pressure [mb]
real(sp),intent(out):: DecDeg           ! declination [deg] (+23.5 deg on 22 June)
real(sp),intent(out):: DayLength        ! day length (dawn:dusk)    [hr]
real(sp),intent(out):: SolarNorm        ! (daily solar)/(solar const): geometry
real(sp),intent(out):: SolarFracObs     ! (obs daily solar)/(solar geometry value)
real(sp),intent(out):: SolarTotChk      ! day-integrated solar [MJ/m2/day] (check)
real(sp),intent(out):: WindAvgChk       ! day-average wind (check)
real(sp),intent(out):: RelHumTempMin    ! rel humidity at minimum air temp (check)
real(sp),intent(out):: PhiSd(ntime)     ! downward solar irradiance [W/m2]
real(sp),intent(out):: coszen(ntime)     ! cos(zentih angle) (vh 12/2/08)
real(sp),intent(out):: PhiLd(ntime)     ! down longwave irradiance  [W/m2]
real(sp),intent(out):: Precip(ntime)    ! precip [m/s]
real(sp),intent(out):: Wind(ntime)      ! wind   [m/s]
real(sp),intent(out):: Temp(ntime)      ! temp   [degC]
real(sp),intent(out):: VapPmb(ntime)    ! vapour pressure [mb]
real(sp),intent(out):: Pmb(ntime)       ! pressure [mb]
! Local variables
real(sp):: LatDeg1, LatRad, YearRad, DecRad, TanTan, HDLRad, TimeNoon, TimeRad
real(sp):: TempAmp, TempAvg, TimeMaxTemp, TimeRadMaxTemp
real(sp):: RatioWindLiteDark, WindDark, WindLite, adjust_fac(ntime)
integer(i4b):: it
!-------------------------------------------------------------------------------

! -------------------------
! Downward solar irradiance
! -------------------------
LatDeg1 = sign(min(abs(LatDeg),89.9), LatDeg)   ! avoid singularity at pole
LatRad  = LatDeg1*Pi/180.0                      ! latitude in radians
YearRad = 2.0*Pi*(YearDay-1)/365.0              ! day of year in radians

! DecRad = Declination in radians (+23.5 deg on 22 June, -23.5 deg on 22 Dec):
DecRad = 0.006918 - 0.399912*cos(YearRad) + 0.070257*sin(YearRad)     &
         - 0.006758*cos(2.0*YearRad) + 0.000907*sin(2.0*YearRad)      &
         - 0.002697*cos(3.0*YearRad) + 0.001480*sin(3.0*YearRad)
                                        ! Paltridge and Platt eq [3.7]
DecDeg = (180.0/Pi) * DecRad            ! Declination in degrees

! Daylength: HDLRad = Half Day Length in radians (dawn:noon = noon:dusk):
TanTan = -tan(LatRad)*tan(DecRad)
if (TanTan .le. -1.0) then
  HDLRad = Pi                           ! polar summer: sun never sets
else if (TanTan .ge. 1.0) then
  HDLRad = 0.0                          ! polar winter: sun never rises
else
  HDLRad = acos(TanTan)                 ! Paltridge and Platt eq [3.21]
end if                                  ! (HDLRad = their capital PI)
DayLength = 24.0*2.0*HDLRad / (2.0*Pi)  ! Daylength (dawn:dusk) in hours

! Daily solar irradiance without atmosphere, normalised by solar constant,
! with both energy fluxes in MJ/m2/day, calculated from solar geometry:
SolarNorm =  &                          ! Paltridge and Platt eq [3.22]
  (HDLRad*sin(LatRad)*sin(DecRad) + cos(LatRad)*cos(DecRad)*sin(HDLRad)) / Pi
! Observed daily solar irradiance as frac of value from solar geometry (should be < 1)
SolarFracObs = SolarMJDay / (SolarNorm*SolarConst + 1.0e-6)

! Instantaneous downward hemispheric solar irradiance PhiSd
DO it=1,ntime                           ! 0=0000, ntime=24000
  TimeNoon = real(it)/real(ntime) - 0.5 ! Time in day frac (-0.5 to 0.5, zero at noon)
  TimeRad  = 2.0*Pi*TimeNoon            ! Time in day frac (-Pi to Pi, zero at noon)
  if (abs(TimeRad) .lt. HDLRad) then    ! day: sun is up
    PhiSd(it) = (SolarMJDay/SolarNorm)  &   ! PhiSd [MJ/m2/day]
      * ( sin(DecRad)*sin(LatRad) + cos(DecRad)*cos(LatRad)*cos(TimeRad) )
                                        ! Paltridge and Platt eq [3.4]
	coszen(it) = ( sin(DecRad)*sin(LatRad) + cos(DecRad)*cos(LatRad)*cos(TimeRad) )
  else                                  ! night: sun is down
    PhiSd(it) = 0.0
	coszen(it) = 0.0
  end if
END DO
SolarTotChk = sum(PhiSd)/real(ntime)    ! check: day-integrated solar [MJ/m2/day] 
PhiSd(:)    = PhiSd(:)*1e6/SecDay       ! Convert PhiSd: [MJ/m2/day] to [W/m2]

! -------------
! Precipitation
! -------------
Precip(:) = PrecipDay/24.0/3600.       ! correct units (vh 12/02/08)

! ----
! Wind
! ----
Wind(:)   = WindDay
RatioWindLiteDark = 3.0                 ! (daytime wind) / (nighttime wind)
WindDark  = WindDay / ( (24.0-DayLength)/24.0 + RatioWindLiteDark*DayLength/24.0 )
WindLite  = WindDay / ( (1.0/RatioWindLiteDark)*(24.0-DayLength)/24.0 + DayLength/24.0 )
DO it=1,ntime                           ! 0=0000, ntime=24000
  TimeNoon = (it-0.5)/real(ntime) - 0.5 ! Time in day frac (-0.5 to 0.5, zero at noon)
  TimeRad  = 2.0*Pi*TimeNoon            ! Time in day frac (-Pi to Pi, zero at noon)
  if (abs(TimeRad) .le. HDLRad) then    ! day: sun is up
    Wind(it) = WindLite
  else                                  ! night: sun is down
    Wind(it) = WindDark
  end if
END DO
WindAvgChk = sum(Wind)/real(ntime)

! -----------
! Temperature
! -----------
TempAvg        = 0.5 * (TempMaxDay + TempMinDay)
TempAmp        = 0.5 * (TempMaxDay - TempMinDay)
TimeMaxTemp    = 16.00                  ! tunable parameter
TimeRadMaxTemp = 2.0*Pi * ((TimeMaxTemp-12.0)/24.0) ! Time of max temp (day radians)
DO it=1,ntime                           ! 0=0000, ntime=24000
  TimeNoon = real(it)/real(ntime) - 0.5 ! Time in day frac (-0.5 to 0.5, zero at noon)
  TimeRad  = 2.0*Pi*TimeNoon            ! Time in day frac (-Pi to Pi, zero at noon)
  Temp(it) = TempAvg + TempAmp * cos(TimeRad-TimeRadMaxTemp)
END DO

! -----------------------------------
! Water Vapour Pressure, Air Pressure
! -----------------------------------
VapPmb(:) = VapPmbDay
Pmb(:)    = PmbDay
RelHumTempMin = VapPmbDay / ESatf(TempMinDay)

! ----------------------------
! Downward longwave irradiance
! ----------------------------
PhiLd(:) = 335.97 * (((Temp(:) + 273.16) / 293.0)**6)   ! [W/m2] (Swinbank 1963)
!



END SUBROUTINE SubDiurnalMet

END MODULE SubDiurnalMetModuleOld

!*******************************************************************************
!*******************************************************************************
MODULE WaterDynModule
!-------------------------------------------------------------------------------
! * This module contains all subroutines for WaterDyn dynamic and observation
!   models at a single time step.
! * Everything in WaterDynModule is private unless otherwise declared
!-------------------------------------------------------------------------------
USE TypeDef
USE DateFunctions
USE Utils
USE define_types                     ! for CABLE 
USE input_module
!USE parameter_module			! for CABLE
USE output_module, ONLY: create_restart,open_output_file, &
       write_output,close_output_file
implicit none
PRIVATE
! * Public subroutines
PUBLIC:: WaterDynStep, ReadControlFile, InitStores, InitVV
! * Public sizes for XX, FF, DD, ZZ, AA, MM,RR, UU, VV arrays
integer(i4b),public:: np                    ! number of spatial points (computed)
integer(i4b),public:: ntile                    ! number of tiles per spatial point (specified)
integer(i4b),public:: nc                    ! number of catchments     (computed)
integer(i4b),parameter,public:: ntt  = 5    ! number of time variables
integer(i4b),parameter,public:: nxx  = 40    ! number of stores
integer(i4b),parameter,public:: nff  = 23   ! number of fluxes and dXXdt
integer(i4b),parameter,public:: ndd  = 30   ! number of derived quantities
integer(i4b),parameter,public:: nzzP = 12   ! number of predicted point obs types
integer(i4b),parameter,public:: naaP = 12    ! number of actual point obs  types
integer(i4b),parameter,public:: nzzPh = 35   ! number of predicted point obs types (hourly)
integer(i4b),parameter,public:: naaPh = 12    ! number of actual point obs  types (hourly)
integer(i4b),parameter,public:: nzzC = 4    ! number of predicted catchment obs types
integer(i4b),parameter,public:: naaC = 4    ! number of actual catchment obs types
integer(i4b),parameter,public:: nmm  = 6    ! number of met variables
integer(i4b),parameter,public:: nrr  = 5    ! number of remote sensing driver variables
integer(i4b),parameter,public:: nhMM  = 8   ! number of hourly met variables
integer(i4b),parameter,public:: nuu  = 63   ! number of spatially uniform parameters
integer(i4b),parameter,public:: nvv  = 72   ! number of spatially variable parameters
integer(i4b),parameter,public:: nAvType = 3 ! number of time av types (month, ann, run)
! * Public parameters read from control file, available to host routine
real(sp),      public:: DelT            ! time step (s)
integer(i4b),  public:: nEnsemble       ! number of ensemble members
integer(i4b),  public:: nspin           ! number of spin-up cycles
integer(i4b),  public:: ispin           !  spin-up cycles counter
real(sp),      public:: TTstart(ntt)    ! values of TT at start of run
real(sp),      public:: UUinitCtl(nuu)  ! spatially uniform (UU) params     (CtlFile values)
real(sp),      public:: UUpinitCtl(nuu) ! error (sqrt(diag(PP)) for UU      (CtlFile values)
real(sp),      public:: UUqCtl(nuu)     ! model error for UU                (CtlFile values)
real(sp),      public:: UUminCtl(nuu)   ! lower bound for PEST for UU       (CtlFile values)
real(sp),      public:: UUmaxCtl(nuu)   ! upper bound for PEST for UU       (CtlFile values)
integer(i4b),  public:: UUSensFlag(nuu) ! flag: parameters selected for sensitivity runs
integer(i4b),  public:: VVSensFlag(nvv) ! flag: parameters selected for sensitivity runs (spatially variable params)
integer(i4b),  public:: UUPESTFlag(nuu) ! flag: parameters selected for PEST runs
integer(i4b),  public:: UUKFFlag(nuu)   ! flag: parameters selected for KF runs
character,     public:: UUabsrel(nuu)   ! 'a' or 'r' for absolute or relative error
real(sp),      public:: UUscale(nuu)    ! scale factor for UU (divide by this)
character(3),  public:: UUtransf(nuu)   ! transformation for UU ('nil' or 'log')
real(sp),      public:: MMmult(nmm)     ! multiplier for MM
real(sp),      public:: MMoffset(nmm)   ! offset for MM (MM = MMmult*MM + MMoffset)
real(sp),      public:: XXinitCtl(nxx)  ! XXinit(np,nxx)                    (CtlFile values)
real(sp),      public:: XXpinitCtl(nxx) ! error (sqrt(diag(PP)) for XXinit  (CtlFile values)
real(sp),      public:: XXqCtl(nxx)     ! model error for XX                (CtlFile values)
character(1),  public:: XXabsrel(nxx)   ! 'a' or 'r' for absolute or relative error
real(sp),      public:: XXscale(nxx)    ! scale factor for XX (divide by this)
character(3),  public:: XXtransf(nxx)   ! transformation for XX ('nil' or 'log')
real(sp),      public:: ZZPq(nzzP)      ! observation model error (as sd) for ZZP
real(sp),      public:: ZZPhq(nzzPh)      ! observation model error (as sd) for ZZPh
real(sp),      public:: ZZCq(nzzC)      ! observation model error (as sd) for ZZC
real(sp),      public:: AAPr(naaP)      ! actual observation error (as sd) for AAP
real(sp),      public:: AAPhr(naaPh)      ! actual observation error (as sd) for AAPh
real(sp),      public:: AACr(naaC)      ! actual observation error (as sd) for AAC
character(1),  public:: AAPabsrel(naaP) ! 'a' or 'r' for absolute or relative error
character(1),  public:: AAPhabsrel(naaPh) ! 'a' or 'r' for absolute or relative error
character(1),  public:: AACabsrel(naaC) ! 'a' or 'r' for absolute or relative error
real(sp),      public:: AAPrlow(naaP)   ! lower threshold for relative obs error (point)
real(sp),      public:: AAPhrlow(naaPh)   ! lower threshold for relative obs error (point)
real(sp),      public:: AACrlow(naaC)   ! lower threshold for relative obs error (catchment)
real(sp),      public:: AAPscale(naaP)  ! scale factor for AAP (divide by this)
real(sp),      public:: AAPhscale(naaPh)  ! scale factor for AAPh (divide by this)
real(sp),      public:: AACscale(naaC)  ! scale factor for AAC (divide by this)
character(3),  public:: AAPtransf(naaP) ! transformation for AAP ('nil' or 'log')
character(3),  public:: AAPhtransf(naaPh) ! transformation for AAP ('nil' or 'log')
character(3),  public:: AACtransf(naaC) ! transformation for AAC ('nil' or 'log')
integer(i4b),  public:: AAPflag(naaP)   ! observation flag (0=NotUsed; 1=Used; 2=TimeAv)
integer(i4b),  public:: AAPhflag(naaPh)   ! observation flag (0=NotUsed; 1=Used; 2=TimeAv)
integer(i4b),  public:: AACflag(naaC)   ! observation flag (0=NotUsed; 1=Used; 2=TimeAv)
integer(i4b),  public:: RRflag(nrr)     ! remote senising driver flag (0=Not Used;, 1 = time series; 2= climatology)
character(1),  public:: AAPdaymth(naaP) ! flag for day or month or climatology ('d' or 'm' or 'c')
character(1),  public:: AAPhdaymth(naaPh) ! flag for day or month ('d' or 'm')
character(1),  public:: AACdaymth(naaC) ! flag for day or month ('d' or 'm')
character(16), public:: TTname(ntt,2)   ! names, units: time variables
character(16), public:: UUName(nuu,2)   ! names, units: UU params
character(16), public:: VVName(nvv,2)   ! names, units: VV params
character(16), public:: MMName(nmm,2)   ! names, units: Met variables
character(16), public:: RRName(nrr,2)   ! names, units: Remote sensing driver variables
character(16), public:: XXName(nxx,2)   ! names, units: State variables
character(16), public:: FFName(nff,2)   ! names, units: Flux variables
integer(i4b),  public:: FFSensFlag(nff) ! flag: fluxes selected for sensitivity runs
character(16), public:: DDName(ndd,2)   ! names, units: Diag variables
character(16), public:: ZZPName(nzzP,2) ! names, units: Predicted point observables
character(16), public:: ZZPhName(nzzPh,2) ! names, units: Predicted point observables
character(16), public:: ZZCName(nzzC,2) ! names, units: Predicted catchment observables
character(16), public:: AAPName(naaP,2) ! names, units: Actual point obs
character(16), public:: AAPhName(naaPh,2) ! names, units: Actual point obs
character(16), public:: AACName(naaC,2) ! names, units: Actual catchment obs
! * Public model control variables, read from control file
real(sp),      public:: MMConst(nmm)        ! default uniform met
real(sp),      public:: RRConst(nmm)        ! default uniform remote sensing driver
real(sp),      public:: VVConst(nvv)        ! default uniform VV params
character(200),public:: DomName             ! domain name (eg Adelong)
character(3),  public:: RunID               ! 3-char run ID (eg 05a)
character(200),public:: Comment             ! comment (second title line)
character(200),public:: DataInRootDir       ! input  data root directory
character(200),public:: DataOutRootDir      ! output data root directory
character(200),public:: MetFileName(nmm)    ! met filenames (imm=1,2,3,4,5,6 for Solar,Rain,Tmax,Tmin,vph09,vph15)
character(200),public:: RRFileName(nrr)    ! remote sensing driver filenames
character(200),public:: NDVIFileName        ! NDVI filename
character(200),public:: CO2FileName         ! CO2 filename
character(200),public:: LSTFileName         ! LST filename
character(200),public:: LSTTimeFileName     ! LSTTime filename
character(200),public:: LSTAngleFileName    ! LSTAngle filename
character(200),public:: phiRnethFileName       ! phiRneth filename
character(200),public:: phiHhFileName       ! phiHh filename
character(200),public:: phiEhFileName       ! phiEh filename
character(200),public:: phiGhFileName       ! phiGh filename
character(200),public:: TsoilhFileName       ! Tsoilh filename
character(200),public:: phiNEEhFileName      ! phiNEEh filename
character(200),public:: phiHhTimeFileName   ! phiHhtime filename
character(200),public:: phiRnetFileName     ! phiRnet filename
character(200),public:: phiHFileName        ! phiH filename
character(200),public:: phiEFileName        ! phiE filename
character(200),public:: phiNEEFileName      ! phiNEE filename
character(200),public:: sm0008FileName      ! sm0008 filename
character(200),public:: sm0090FileName      ! sm0090 filename
character(200),public:: sm0030FileName      ! sm0030 filename
character(200),public:: sm3060FileName      ! sm3060 filename
character(200),public:: sm6090FileName      ! sm6090 filename
character(200),public:: LocalMetFileName    ! Local Met data filename

integer(i4b),  public:: ForwardModelFlag    ! (0,1,2) = (Waterdyn, CABLE, CABLE + soil-litter)
integer(i4b),  public:: CableParamFlag      ! (0,1,2) = (single tile estimate, woody/grassy tile, CABLE defaults) 
integer(i4b),  public:: SensTestFlag        ! (0,1,2,3) = (forward run, param sensitivity test (UU),param sensitivity test (VV),monthly param sensitivity test (VV) )
integer(i4b),  public:: DiagsFlag           ! (0,1) = screen diags (off, on)
integer(i4b),  public:: ipDiag              ! (0, >=1) for TS output = (domain-av, point ipDiag)
integer(i4b),  public:: DomTSOutputFlag     ! (0,1,2,3,4) = (no,run,ann,mth,day) domain-av TS
integer(i4b),  public:: MapTSOutputFlag     ! (0,1,2,3,4) = (no,run,ann,mth,day) DD TS (dated binary file, spatial)
integer(i4b),  public:: CatchTSOutputFlag   ! (0,1,2,3,4) = (no,run,ann,mth,day) catchment TS
integer(i4b),  public:: MapOutputFlag       ! (0,1,2,3,4) = (no,run,ann,mth,day) maps
integer(i4b),  public:: MapOutputFlag0       !  original flag: can be overwritten if daily output required for specified period
integer(i4b),  public:: Ihour_output(24)    ! specifies which hours to output for AAPh and ZZPh
integer(i4b),  public:: PESTOutputFlag      ! (0=no, 2,3,4=PESToutput)
integer(i4b),  public:: WriteInitStoresFlag ! (0,1) = (N,Y) to write last XX to file for later init
integer(i4b),  public:: XXInitFileFlag      ! (0,1) = (N,Y) to initialise XX from file
integer(i4b),  public:: FracVegFlag         ! (0,1) = FracV from (external VegFPC, internal CLea)
integer(i4b),  public:: HourlyOutputFlag    ! (0,1) = (N,Y) hourly output from CABLE to file
integer(i4b),  public:: UseLocalMetDataFlag !(0,1,2) = (N,Y) use local met forcing data if available (2== repeat met over all pixels)
integer(i4b),  public:: UseLAIFlag          !(0,1,2) = LAI derived from Seawifs FAPAR, MODIS monthly LAI from VV param file, woody/grassy LAI from obs file
character(20), public:: MMFileID(nmm)       ! spatial MM data filenames: string ID
character(50), public:: RRFileID(nrr)       ! spatial RR data filenames: string ID
character(200),public:: VVFileName(nvv)     ! spatial VV data filenames
character(200),public:: XXInitFileName(nxx) ! spatial XX data filenames for initialisation
character(200),public:: CtlFileName         ! control file name (to be seen by host program)
integer(i4b),  public:: VVSpatialFlag(nvv)  ! (0,1) = (VV from VVConst, VV from vector map)
integer(i4b),  public:: MMMapOPFlag(nmm)    ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: RRMapOPFlag(nrr)    ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: XXMapOPFlag(nxx)    ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: FFMapOPFlag(nff)    ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: DDMapOPFlag(ndd)    ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: MMTSOPFlag(nmm)    ! output flags for specific time-series (0=N,1=run,2=ann,3= monthly 4=daily 5 = hourly)
integer(i4b),  public:: DDTSOPFlag(ndd)    ! output flags for specific time-series (0=N,1=run,2=ann,3= monthly 4=daily)
integer(i4b),  public:: XXTSOPFlag(nxx)    ! output flags for specific time-series (0=N,1=run,2=ann,3= monthly 4=daily)
integer(i4b),  public:: FFTSOPFlag(nff)    ! output flags for specific time-series (0=N,1=run,2=ann,3= monthly 4=daily)
integer(i4b),  public:: ZZPMapOPFlag(nzzP)  ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: AAPMapOPFlag(naaP)  ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: ZZPhMapOPFlag(nzzPh)  ! output flags for specific maps (0=N,1=Y)
integer(i4b),  public:: AAPhMapOPFlag(naaPh)  ! output flags for specific maps (0=N,1=Y)
! * Housekeeping variables
integer(i4b),  public:: nSteps              ! number of (daily) time steps
integer(i4b),  public:: ntime               ! number of subdiurnal timesteps 
integer(i4b),  public:: nMonths             ! number of months
integer(i4b),  public:: InitFlag            ! flag to call InitModel (0,1 = Y,N)
integer(i4b),  public:: InitFlagCO2            ! flag to call Init CO2 file (0,1 = Y,N)
integer(i4b),  public:: EndFlag             ! flag to indicate last timestep in model run (0,1 = Y,N)
type(dmydate), public:: InitDate, StartDate, EndDate
integer(i4b),  public:: nSteps_spin             ! number of (daily) time steps in spin cycle
integer(i4b),  public:: nMonths_spin             ! number of months in spin-cycle
integer(i4b),  public:: nSteps_run             ! number of (daily) time steps in spin cycle
integer(i4b),  public:: nMonths_run             ! number of months in spin-cycle
                                            ! init, start, end dates (InitDate = StartDate-1)
type(dmydate), public:: StartDate_spin, EndDate_spin  ! start and end spin dates
type(dmydate), public:: StartDate_run, EndDate_run  ! start and end run dates
type(dmydate), public:: StartDate_DailyOutput, EndDate_DailyOutput                                            ! start, end dates for hourly output
character(200),public:: GenFileName         ! generic filename for building char strings
character(200),public:: GenFileName2        ! generic filename for building char strings
character(200),public:: GenFileName3        ! generic filename for building char strings
character(200),public:: GenFileName4        ! generic filename for building char strings
character(200),public:: GenFileName5        ! generic filename for building char strings
character(200),public:: GenFileName6       ! generic filename for building char strings
character(200),public:: GenFileName7       ! generic filename for building char strings
logical(lgt),  public:: GenFileExistFlag    ! generic flag for existence of data file
character(8),  public:: TTDateChar          ! date as character string (YYYYMMDD)
character(10),  public:: TTDateTimeChar          ! date as character string (YYYYMMDDHH)
real(sp),      public:: ErrVal = -999.0     ! error code
real(sp),      public:: cpuTimeBegin, cpuTimeNow
                                            ! CPU time markers
type(dmydate),allocatable:: NextRRDate(:)      ! dates of next RR records
type(dmydate),allocatable:: NextRRDatetmp(:)      ! dates of first RR records when after StartDate
type(dmydate),public :: NextCO2Date              ! date of next CO2 record
type(dmydate)  :: NextNDVIDate              ! date of next NDVI record
type(dmydate)  :: NextLSTDate               ! date of next LST record
type(dmydate)  :: NextLSTTimeDate           ! date of next LSTTime record
type(dmydate)  :: NextLSTAngleDate          ! date of next LSTAngle record
type(dmydate)  :: NextAphiRnethDate         ! date of next phiRneth record
type(dmydate)  :: NextAphiHhDate            ! date of next phiHh record
type(dmydate)  :: NextAphiEhDate            ! date of next phiEh record
type(dmydate)  :: NextAphiNEEhDate          ! date of next phiNEEh record
type(dmydate)  :: NextAphiGhDate            ! date of next phiGh record
type(dmydate)  :: NextATsoilhDate            ! date of next Tsoilh record
type(dmydate)  :: NextAphiHhTimeDate        ! date of nextphiHhTime record
type(dmydate)  :: NextAphiRnetDate          ! date of next Rnet record
type(dmydate)  :: NextAphiHDate             ! date of next phiH record
type(dmydate)  :: NextAphiEDate             ! date of next phiE record
type(dmydate)  :: NextAphiNEEDate           ! date of next NEE record
type(dmydate)  :: NextAsm0008Date           ! date of next sm0008 record
type(dmydate)  :: NextAsm0090Date           ! date of next sm0090 record
type(dmydate)  :: NextAsm0030Date           ! date of next sm0030 record
type(dmydate)  :: NextAsm3060Date           ! date of next sm3060 record
type(dmydate)  :: NextAsm6090Date           ! date of next sm6090 record
real(sp)       :: SolarUTCOffset = 10.0     ! Solar time offset from UTC (hours)
real(sp)           :: NextCO2              ! next CO2 record
real(sp),allocatable:: NextRR(:,:)			! np*nrr array: data in next RR records
real(sp),allocatable:: RRtmp(:,:)			! np*nrr array: data in first RR record when after StartDate
real(sp),allocatable:: NextANDVI(:)         ! np-vector: data in next NDVI record
real(sp),allocatable:: NextALST(:)          ! np-vector: data in next LST record
real(sp),allocatable:: NextALSTTime(:)      ! np-vector: data in next LSTTime record
real(sp),allocatable:: NextALSTAngle(:)     ! np-vector: data in next LSTAngle record
real(sp),allocatable:: NextAphiRneth(:)     ! np-vector: data in next phiRneth record
real(sp),allocatable:: NextAphiHh(:)        ! np-vector: data in next phiHh record
real(sp),allocatable:: NextAphiEh(:)        ! np-vector: data in next phiEh record
real(sp),allocatable:: NextAphiNEEh(:)        ! np-vector: data in next phiNEEh record
real(sp),allocatable:: NextAphiGh(:)        ! np-vector: data in next phiGh record
real(sp),allocatable:: NextATsoilh(:)        ! np-vector: data in next phiGh record
real(sp),allocatable:: NextAphiHhtime(:)    ! np-vector: data in next phiHhTime record
real(sp),allocatable:: NextAphiRnet(:)      ! np-vector: data in next phiRnet record
real(sp),allocatable:: NextAphiH(:)         ! np-vector: data in next phiH record
real(sp),allocatable:: NextAphiE(:)         ! np-vector: data in next phiE record
real(sp),allocatable:: NextAphiNEE(:)       ! np-vector: data in next phiNEE record
real(sp),allocatable:: NextAsm0008(:)       ! np-vector: data in next sm0008 record
real(sp),allocatable:: NextAsm0090(:)       ! np-vector: data in next sm0090 record
real(sp),allocatable:: NextAsm0030(:)       ! np-vector: data in next sm0030 record
real(sp),allocatable:: NextAsm3060(:)       ! np-vector: data in next sm3060 record
real(sp),allocatable:: NextAsm6090(:)       ! np-vector: data in next sm6090 record
! * Accumulators for space-time averages of np-arrays for TS output
!   * First dimension is 1 as for MMSpStat
!   * Moved from WaterDynStep to WaterDynModule and made PUBLIC, for access by host
!     program to permit sensitivity tests (14-mar-07, at V14a)
real(sp),public:: MMSpTAv(1,nmm,nAvType),   MMSpTAvCount(1,nmm,nAvType)
real(sp),public:: RRSpTAv(1,nrr,nAvType),   RRSpTAvCount(1,nrr,nAvType)
real(sp),public:: XXSpTAv(1,nxx,nAvType),   XXSpTAvCount(1,nxx,nAvType)
real(sp),public:: FFSpTAv(1,nff,nAvType),   FFSpTAvCount(1,nff,nAvType)
real(sp),public:: FFSpTAvTemp(1,nff,nAvType)
real(sp),public:: DDSpTAv(1,ndd,nAvType),   DDSpTAvCount(1,ndd,nAvType)
real(sp),public:: ZZPSpTAv(1,nzzP,nAvType), ZZPSpTAvCount(1,nzzP,nAvType)
real(sp),public:: AAPSpTAv(1,naaP,nAvType), AAPSpTAvCount(1,naaP,nAvType)
!real(sp),public:: ZZPhSpTAv(1,nzzP,ntime,nAvType), ZZPhSpTAvCount(1,nzzP,ntime,nAvType)
!real(sp),public:: AAPhSpTAv(1,naaP,ntime,nAvType), AAPhSpTAvCount(1,naaP,ntime,nAvType)
! * Accumulators for catchment time averages for TS output over nc catchments
real(sp),public,allocatable:: FFCTAv(:,:,:), FFCTAvCount(:,:,:)   ! (np,nFF,nAvType)
real(sp),public,allocatable:: ZZCTAv(:,:,:), ZZCTAvCount(:,:,:)   ! (np,nzzC,nAvType)
real(sp),public,allocatable:: AACTAv(:,:,:), AACTAvCount(:,:,:)   ! (np,naaC,nAvType)
real(sp),public,allocatable:: ZZPhSpTAv(:,:,:,:), ZZPhSpTAvCount(:,:,:,:)
real(sp),public,allocatable:: AAPhSpTAv(:,:,:,:), AAPhSpTAvCount(:,:,:,:)
! * Catchment properties for nc catchments
integer(i4b):: iCatchID(999)        ! integer unique catchment ID numbers
real(sp)    :: TZRunCat(999)       ! smoothing time for catchment daily runoff
real(sp)    :: TZLchCat(999)       ! smoothing time for catchment daily leaching
logical(lgt):: DisCMDataExist(999)  ! existence flag for monthly discharge data
logical(lgt):: DisCDDataExist(999)  ! existence flag for daily discharge data

real(sp),public,allocatable :: tile_area(:,:)  ! tile area for each veg type in each grid cell (np,ntile)
integer(i4b),public,allocatable :: tile_index(:,:)  ! index for each tile within a vector dimension np*ntile (ntile*np)

TYPE (met_type),public 	:: met_old

CONTAINS

!*******************************************************************************

SUBROUTINE WaterDynStep (Time0,Time0_met, UU,VV,XX0,   XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
!-------------------------------------------------------------------------------
! * Advance model one predition step from T0 to T1, and return:
!   (1) priors XX1 and (ZZP1, ZZPh1, ZZC1)
!   (2) actual observations (AAP1, AAPh1, AAC1)
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
USE DateFunctions
USE PointerModule
USE SubDiurnalMetModule
USE Constants

USE define_dimensions, ONLY:r_1,r_2,i_d,ms,mp_patch,mp, max_vegpatches           ! for CABLE
USE main_cable_module					     ! for CABLE
USE checks_module, ONLY: ranges					 ! for CABLE
USE physical_constants, ONLY: tfrz
USE parameter_module				 ! for CABLE
USE canopy_module					 ! for CABLE
USE define_types					 ! for CABLE
USE io_variables, ONLY: logn,filename,nvegt,nsoilt,leaps, &
       verbose,fixedCO2,output,check,patchout,model_structure_flag,&
	   longitude, latitude, smoy,defaultLAI,landpt,patch,nsoilt,nvegt, &
	   xdimsize, ydimsize, timeunits, time_coord, lat_all, lon_all

implicit none

  TYPE (air_type)       :: air  ! air property variables
  TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
  TYPE (canopy_type)    :: canopy ! vegetation variables
  TYPE (met_type)       :: met  ! met input variables
  TYPE (balances_type)  :: bal  ! energy and water balance variables
  TYPE (radiation_type) :: rad  ! radiation variables
  TYPE (roughness_type) :: rough ! roughness varibles
  TYPE (soil_parameter_type) :: soil ! soil parameters	
  TYPE (soil_snow_type) :: ssoil ! soil and snow variables
  TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
  TYPE (veg_parameter_type) :: veg  ! vegetation parameters	 
  TYPE (model_structure_type) :: model_structure  ! vegetation parameters
 ! TYPE (units_type) :: units 
! Global variables, seen by the host program, with arrays passed as assumed-shape
real(sp),intent(in) :: Time0        ! time [days] at start of step, from host program
real(sp),intent(in) :: Time0_met       ! time [days] at start of step, from host program
real(sp),intent(inout) :: UU(:)        ! UU(nuu) = spatially uniform params (adjustable)
real(sp),intent(in) :: VV(:,:)        ! VV(np,nvv) = spatially varaiable params (adjustable)
real(sp),intent(in) :: XX0(:,:)     ! XX0(np,nxx) = state variables at T0
real(sp),intent(out):: XX1(:,:)     ! XX1(np,nxx) = state variables at T1
real(sp),intent(out):: ZZP1(:,:)    ! ZZP1(np,nzzP,ntime) = predicted point obs at T1
real(sp),intent(out):: ZZPh1(:,:,:)    ! ZZP1(np,nzzP,ntime) = predicted point obs at T1
real(sp),intent(out):: ZZC1(:,:)    ! ZZC1(nc,nzzC) = predicted catchment obs at T1
real(sp),intent(out):: AAP1(:,:)    ! AAP1(np,naaP,ntime) = actual point obs at T1
real(sp),intent(out):: AAPh1(:,:,:)    ! AAP1(np,naaP,ntime) = actual point obs at T1
real(sp),intent(out):: AAC1(:,:)    ! AAC1(nc,naaC) = actual catchment obs at T1
real(sp),intent(out):: TTime(:)     ! TTime(ntt)    = model time variables
! Local variables (private to model)
! * Generic arrays passed explicitly to subroutines of WaterDynStep, used as 
!   targets within these subroutines for specific variable names (via PointAll)
! * These arrays have the save attribute so their values are held between calls
!real(sp),allocatable,save:: VV(:,:) ! VV(np,nvv) = spatially variable params
real(sp),allocatable,save:: FF(:,:) ! FF(np,nff) = fluxes and dXXdt
real(sp),allocatable,save:: DD(:,:) ! DD(np,ndd) = derived quantities
real(sp),allocatable,save:: MM(:,:) ! MM(np,nmm) = met variables
real(sp),allocatable,save:: RR(:,:) ! RR(np,nrr) = met variables
real(sp),allocatable,save:: MMprev(:,:) ! MM(np,nmm) = met variables
real(sp),allocatable,save:: MMnext(:,:) ! MM(np,nmm) = met variables
real(sp),allocatable,save:: hMM(:,:,:) !hMM(np,nhmm,ntime) = met variables
real(sp),allocatable,save:: hMMnext(:,:,:) !hMM(np,nhmm,ntime) = met variables
real(sp),allocatable,save:: hMMlocal(:,:,:) !hMM(np,nhmm,ntime) = met variables
! Arrays for spatial statistics of np-arrays (mean, sd, max, min, count)
! * First dimension (1) is a quasi-spatial dimension so that MM(1,nmm,1) has shape
!   compatible with MM(np,nmm), for call to TimeAverage
real(sp):: VVSpStat(1,nvv,5), XXinitSpStat(1,nxx,5),    &
           MMSpStat(1,nmm,5),RRSpStat(1,nrr,5), XXSpStat(1,nxx,5),  &
           FFSpStat(1,nff,5), DDSpStat(1,ndd,5),        &
           ZZPSpStat(1,nzzP,5), AAPSpStat(1,naaP,5) ,    &
		   ZZPhSpStat(1,nzzP,ntime,5), AAPhSpStat(1,naaP,ntime,5)
! Accumulators for space-time averages of np-arrays for TS output
! * declared in WaterDynModule to permit sensitivity tests
! Accumulators for point time averages for np-map output
real(sp),allocatable,save:: MMTAv(:,:,:),  MMTAvCount(:,:,:)    ! (np,nmm,nAvType)
real(sp),allocatable,save:: RRTAv(:,:,:),  RRTAvCount(:,:,:)    ! (np,nrr,nAvType)
real(sp),allocatable,save:: XXTAv(:,:,:),  XXTAvCount(:,:,:)    ! (np,nxx,nAvType)
real(sp),allocatable,save:: FFTAv(:,:,:),  FFTAvCount(:,:,:)    ! (np,nff,nAvType)
real(sp),allocatable,save:: DDTAv(:,:,:),  DDTAvCount(:,:,:)    ! (np,ndd,nAvType)
real(sp),allocatable,save:: ZZPTAv(:,:,:), ZZPTAvCount(:,:,:)   ! (np,nzzP,nAvType)
real(sp),allocatable,save:: AAPTAv(:,:,:), AAPTAvCount(:,:,:)   ! (np,naaP,nAvType)
real(sp),allocatable,save:: ZZPhTAv(:,:,:,:), ZZPhTAvCount(:,:,:,:)   ! (np,nzzP,ntime,nAvType)
real(sp),allocatable,save:: AAPhTAv(:,:,:,:), AAPhTAvCount(:,:,:,:)   ! (np,naaP,ntime,nAvType)
! Other variables
type(dmydate):: TTDate			    ! current date
type(dmydate):: TTDateNext			! Next date
real(sp)            :: TTimeNext(ntt)   ! model time variables for next step
integer(i4b) :: AvType              ! type of time av: (0,1,2,3) = (day,month,year,run)
integer(i4b) :: ixx, ip             ! indices
integer(i4b),save:: Time0Last       ! host time (Time0) when WaterDynStep last called
! Other variables for CABLE
integer(i_d),save ::ktau      ! cumulative counter for daily timestep
integer(i_d),save ::kstart, kend !star and end time-steps
integer(i_d) ::k     ! integer for hourly iteration loop
integer(i_d) ::j     ! integer for tile loop
integer(i_d) ::i     ! integer for tile loop
real(r_1), save:: dels ! CABLE time step (in s)
real(sp):: doy, xphi1(np), xphi2(np), fsoil(np), canopy_wbal(np)
CHARACTER(LEN=99),save :: filename_out ! name of file for CABLE output
CHARACTER(LEN=99),save :: filename_log ! name of file for CABLE log
CHARACTER(LEN=99),save :: filename_hmet ! name of file for CABLE hourly met (testing purpose only)
CHARACTER(LEN=99),save :: filename_LST ! name of file for CABLE hourly LST (testing purpose only)
CHARACTER(LEN=99),save :: filename_roughness ! name of file for CABLE hourly LST (testing purpose only)
INTEGER(i_d),save  :: hmetn                  ! met file unit number
!INTEGER(i_d),save  :: logn                  ! log file unit number
INTEGER(i_d),save  :: LSTn                  ! log file unit number
INTEGER(i_d),save  :: Roughn                  ! log file unit number
INTEGER(i_d) ::       count_daylight(np)        ! counter for daylight hours (required for PhiE and PhiH)
! Local variables for calculation of LAI from FAPAR
real(sp):: FracV(size(XX0,1)), FracVExt1(size(XX0,1)) ! vegetation fraction cover
real(sp):: FracVmax                                 ! maximum FracV
! Local variables for calculating precicted LST
real(sp):: ALSTTime2(np,ntime), ALSTAngle2(np,ntime)
real(sp):: Rtr(np*ntile), Rtg(np*ntile), Rbs(np*ntile), Rbl(np*ntile)  ! special for Tumba resp
INTEGER(i_d) :: irr,iunit,dum
LOGICAL:: istat,isopen

!-------------------------------------------------------------------------------
!   * Set variables in TTime using Time0 = time [days] at start of step, from host program 
TTime(1) = (Time0_met + 1.0)/365.25	! time (y) at end of current step
TTDate   = StartDate + nint(Time0_met)	! current date as dmytype (use date arithmetic)
TTime(2) = TTDate%day
TTime(3) = TTDate%month
TTime(4) = TTDate%year
TTime(5) = 0.0                      ! End of month flag
if (EndMonth(TTDate)) TTime(5) = 1.0
doy = YearDay(TTDate)

TTimeNext(1) = (Time0_met + 1.0)/365.25	! time (y) at end of current step
TTDateNext   = StartDate + nint(Time0_met) +nint(1.0)	! next date as dmytype (use date arithmetic)
TTimeNext(2) = TTDateNext%day
TTimeNext(3) = TTDateNext%month
TTimeNext(4) = TTDateNext%year
TTimeNext(5) = 0.0 
if (EndMonth(TTDateNext)) TTimeNext(5) = 1.0
! * PREPARATION AND DATA INPUT
!   * Model initialisation (first step only)
if (InitFlag == 0) then
  CALL InitModel
 
  if (ForwardModelFlag.ge.1) then
	  ! * Initialise data structures for CABLE
	  if (nspin==0.or.ispin==1) CALL InitCABLE(UU)	  
	  hmetn = 99
	  LSTn=66
	  Roughn = 77
	  filename_hmet    = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // 'hmet_cable.txt'
	  OPEN(hmetn,FILE=filename_hmet)       ! open hmetfile
	  filename_LST    = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // 'LST_cable.csv'
	  filename_roughness    = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // 'Roughness_cable.csv'
	  OPEN(LSTn,FILE=filename_LST)       ! open hmetfile
	  OPEN(Roughn,FILE=filename_roughness)       ! open Roughfile
	  dels = delT                      ! CABLE time step (in s)
	  ktau = 0                          ! reset counter for houly timesteps (CABLE)  
	  kstart = 1
	  kend = ntime   ! number of timesteps in one day
	  mp = np                            ! number of grid cells (CABLE)
	  filename_out = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // 'out_cable.txt'
  endif
  
end if

!   * If host time has changed, read met data 
if (ForwardModelFlag.eq.0) then
	if (Time0 .gt. Time0Last) then 
			CALL ReadObs(TTime, AAP1,AAPh1, AAC1)
			CALL ReadMet (TTime, MM, hMM)
		if (initFlag==0) then
			InitFlag = 1                      ! reset InitFlag
		endif
	endif
elseif (ForwardModelFlag.ge.1) then
	if (Time0 .gt. Time0Last) then 
		if (initFlag==0) then
			CALL ReadObs(TTime, AAP1,AAPh1, AAC1)
			CALL ReadRR(TTime,RR)
			MMprev = MM
			write(*,*) TTime
			CALL ReadMet (TTime, MM, hMMlocal)
			CALL ReadMet (TTimeNext, MMnext, hMMnext)
			
			InitFlag = 1
		elseif (initFlag==1.and.EndFlag==1) then
			MMprev = MM
			MM = MMnext
			hMMlocal = hMMnext
			CALL ReadObs (TTime,  AAP1,AAPh1, AAC1)
			CALL ReadRR(TTime,RR)
			CALL ReadMet (TTimeNext, MMnext, hMMnext)
	
		elseif (initFlag==1.and.EndFlag==0) then
			MMprev = MM
			MM = MMnext
			hMMlocal = hMMnext
			CALL ReadObs (TTime,  AAP1,AAPh1, AAC1)
			CALL ReadRR(TTime,RR)
			!CALL ReadMet (TTimeNext, MMnext, hMMnext)
		endif
	endif
endif

!   * Reset variables in TTime using Time0 = time [days] from start of model run
TTime(1) = (Time0 + 1.0)/365.25	! time (y) at end of current step
TTDate   = StartDate_run + nint(Time0)	! current date as dmytype (use date arithmetic)
TTime(2) = TTDate%day
TTime(3) = TTDate%month
TTime(4) = TTDate%year
TTime(5) = 0.0                      ! End of month flag
if (EndMonth(TTDate)) TTime(5) = 1.0
doy = YearDay(TTDate)

if (ForwardModelFlag.ge.1) then
	if (Time0 .gt. Time0Last) then 
			CALL ReadCO2(TTime,UU)
	endif
endif


! set MapOutputFlag to 4 if required
if (TTDate.ge.StartDate_DailyOutput.and.TTDate.le.EndDate_DailyOutput) then
	MapOutputFlag = 4
else
	MapOutputFlag = MapOutputFlag0
endif		

if (ForwardModelFlag.ge.1) then
	CALL PointAll (TargetFF=FF, TargetXX=XX1, TargetDD=DD, TargetAAPh = AAPh1, TargetZZPh = ZZPh1, TargetZZP = ZZP1)
	! initialise fluxes
	FF = 0.0
	DD = 0.0
	XX1 = XX0
	count_daylight=0
	fsoil = 0.0
	! initialise model observables
	ZZPh1 = errVal
	ZphiRnet = 0.0
	ZphiH = 0.0
	ZphiE=0.0
	ZphiNEE=0.0
	ZphiNPP=0.0
	ZphiRneth = 0.0
	ZphiHh = 0.0
	ZphiEh=0.0
	ZphiNEEh=0.0
	ZphiGh=0.0
	ZTsoilh = 0.0
	
	if (Time0 .gt. Time0Last) then

		if (UseLocalMetDataFlag==0) then   

			CALL GetSubdiurnalMet(TTime,MM,MMprev,MMnext,VV,   hMM)
		elseif (UseLocalMetDataFlag>=1) then                  ! replace hourly met data for weather generator with local met data if available
           if (ALL(hMMLocal.gt.(ErrVal +1)).eq.1) then
				hMM = hMMlocal
		   else
				CALL GetSubdiurnalMet(TTime,MM,MMprev,MMnext,VV,   hMM)
				where (hMMlocal.gt.(ErrVal +1))
					hMM = hMMlocal
				endwhere
		   endif
	    endif


	  ! * Update tile area
	  	if (CableParamFlag.eq.1) then
			tile_area(:,1) = 1.-fWoody
			tile_area(:,2) = fWoody
			patch(1:mp_patch:ntile)%frac = tile_area(:,1)
			patch(2:mp_patch:ntile)%frac = tile_area(:,2)
		else
			tile_area(:,1)=1.0
			patch%frac=1.0
		endif
	
	  ! * Find LAI and veg cover fraction, related by CoeffBeer. Pre-apply rLAIMult to LAI.
		FracVmax  = 1.0 - exp(-CoeffBeer*rLAImax)
		!   * FracV from external VegFPC
		FracVExt1 = coeffPAR * FAPAR(:,nint(TTMonth))       ! externally specified FracV
		where (FracVExt1 < FracVmax)
			rLAIExt = -log(1.0-FracVExt1) / CoeffBeer
		elsewhere
			rLAIExt = rLAImax
		end where
		FracVExt = 1.0 - exp(-CoeffBeer*rLAIExt)            ! and calculate consistent FracV
		!   * FracV from internal CLea (with CLea0 = RhoCLeaf*LeafThick)
		rLAICLea  = min(max(CLea/CLea0, 0.0), rLAImax)      ! constrain to 0 <= rLAI <= rLAImax
		FracVCLea = 1.0 - exp(-CoeffBeer*rLAICLea)          ! and calculate consistent FracV
		!   * select FracV for computation
		if (FracVegFlag <= 0) then
			FracV = FracVExt
			do k=1,ntile
				veg%vlai(k:mp_patch:ntile) = rLAIext 
			
			enddo
		else
			FracV = FracVCLea
			do k=1,ntile
				veg%vlai(k:mp_patch:ntile) = rLAICLea
			enddo
		end if 

		if (UseLAIFlag==1) then
		    do k=1,ntile
				veg%vlai(k:mp_patch:ntile)= max(LAI(:,nint(TTMonth)),0.01)
			enddo
		elseif (UseLAIFlag==2) then  ! separate values for grassy and woody tiles
			where (tile_area(:,2)>1.e-5) ! time series LAI available for woody & grassy components
				veg%vlai(1:mp_patch:ntile) = LAIg/tile_area(:,1)
				veg%vlai(2:mp_patch:ntile) = LAIw/tile_area(:,2)
			elsewhere 
				veg%vlai(1:mp_patch:ntile) = LAIg
				veg%vlai(2:mp_patch:ntile) = LAIw

			endwhere

			where (veg%vlai>6.)
				veg%vlai=6.
			endwhere
		end if


	
	end if ! (Time0 .gt. Time0Last)
	
	where (ALSTTime >= 0.0 .and. ALSTTime <= 24.0)  ! compute ZLST for all points
		ALSTTime2 = ALSTTime
	elsewhere
		ALSTTime2 = 12.0                              ! by setting ALSTTime to default if NA
	end where

	where (ALSTAngle >= 0.0 .and. ALSTAngle <= 90.0)  ! compute ZLST for all points
		ALSTAngle2 = ALSTAngle
	elsewhere
		ALSTAngle2 = 0.0                              ! by setting ALSTAngle to default if NA
	end where


	! Call CABLE 24 times
	do k=kstart,kend
		ktau = ktau + 1

		! assign met variables, same for each tile
		do j=1,ntile
			met%year(j:mp_patch:ntile) = TTDate%year
			met%moy(j:mp_patch:ntile) = TTDate%month
			met%doy(j:mp_patch:ntile) = doy
			met%hod(j:mp_patch:ntile) = k*dels/3600.
			met%fsd(j:mp_patch:ntile) = hFsd(:,k)
			met%fld(j:mp_patch:ntile) = hFld(:,k)
			met%precip(j:mp_patch:ntile) = hPrecip(:,k)*dels/3600. ! convert hourly precip to mm/timestep
			met%precip_s(j:mp_patch:ntile) = 0.0
			met%tk(j:mp_patch:ntile) = hTc(:,k) + 273.15
			met%ua(j:mp_patch:ntile) = hUa(:,k)
			met%qv(j:mp_patch:ntile) = hQv(:,k)
			met%pmb(j:mp_patch:ntile) =  hPmb(:,k)
			met%coszen(j:mp_patch:ntile) = hcoszen(:,k)
		enddo

        where(met%precip.lt.0.0)
			met%precip=0.0
		endwhere
		
		met%tvair = met%tk
		met%tvrad = met%tk
	
	
		
		met%ca = co2a   
		
		met%da = 0.0
		met%dva = 0.0
		veg%extkn = 0.01
		if (1.eq.0) then
		write(hmetn,"(4(i8,','),14(e15.6,','))") TTDATE%year,TTDATE%month, TTDATE%day,k-1, met%fsd(1), met%fld(1), &
		      met%precip(1), met%tk(1), met%ua(1), met%qv(1)*100., met%pmb(1)
		endif
		if (1.eq.0) then
		write(roughn,"(4(i8,','),14(e15.6,','))") TTDATE%year,TTDATE%month, TTDATE%day,k-1, rough%disp(1), rough%z0m(1), &
		      rough%z0soil(1), rough%rt1usa(1), ssoil%rtsoil(1), ssoil%rlitt(1), canopy%tv(1)-273.15, canopy%tscrn(1), met%tk(1)-273.15
		endif



		CALL cable_main(ktau, kstart, kend, dels, air, bgc, canopy, met, &
             & bal, rad, rough, soil, ssoil, sum_flux, veg, nvegt, nsoilt &
			 , model_structure)

! special for Tumba: known soil and plant respiration rates

	 WHERE(hFsd(:,k)<5.)                 ! at night-time
          !all respiration rates in umolm-2s-1
	     Rtr = 0.34*(met%tvair-273.15)+3.93
	     Rtg = 4.1 * Rtr
	     Rbs = 6.9 * Rtr
	     Rbl = 4.1 * Rtr
	
      ELSEWHERE

	     Rtr = 0.34*(met%tvair-273.15)+3.93
	     Rtg = 2.0 * Rtr
	     Rbs = 5.2 * Rtr
	     Rbl = 4.1 * Rtr     
	  ENDWHERE
      
      
	  ! special for Tumba
		!where ((rad%latitude(:).gt.(-35.65-0.01).and.rad%latitude(:).lt.(-35.65+0.01)).and.(rad%longitude(:).gt.(148.15-0.01).and.rad%longitude(:).lt.(148.15+0.01))) 
		!	canopy%frp = (Rtr * 0.0086 + Rtg * 0.0028 + Rbs * 0.010 + Rbl * 0.011) *12./1.e6
		!	canopy%frs=exp(5.822+0.0323*(ssoil%Tgg(:,1)-Tfrz)-2.587*ssoil%wb(:,1)+0.1445*(ssoil%Tgg(:,1)-Tfrz)*ssoil%wb(:,1)) &
         !  *(1000.0/44.0)/3600.0*12./1.e6
		!endwhere
		
		
          ! calculate Surface Temperature as seen by radiometer

			TrVeg(:,k) = (/(sum(canopy%tv(tile_index(i,:))**4*tile_area(i,:)),i=1,np)/)**0.25-Tfrz
			TrSoil(:,k) = (/(sum(ssoil%tss(tile_index(i,:))**4*tile_area(i,:)),i=1,np)/)**0.25-Tfrz

			!Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
			xphi1 = 0.5 - 0.633*veg%xfang(1:mp_patch:ntile) -0.33*veg%xfang(1:mp_patch:ntile)*veg%xfang(1:mp_patch:ntile)
			xphi2 = 0.877 * (1.0 - 2.0*xphi1)
			! needs corecting for multiple tiles
			fsoil = exp(-(xphi1 / cos(ALSTAngle2(:,k)*pi/180.) + xphi2)*veg%vlai(1:mp_patch:ntile))
			fsoil = exp(-(xphi1 / cos(0*pi/180.) + xphi2)*veg%vlai)                ! special for Tumba 
			TrEff(:,k) = ( (1.-fsoil)*canopy%tv**4 + fsoil * ssoil%tss**4 )**0.25 - 273.16
			
			TaZ(:,k) = met%tk(1:mp_patch:ntile)-273.15
			TempAt(:,k) = met%tk(1:mp_patch:ntile)-273.15
			ZphiHh(:,k) = (/(sum(canopy%fh(tile_index(i,:))*tile_area(i,:)),i=1,np)/)
			ZphiRneth(:,k) = (/(sum(canopy%fns(tile_index(i,:))*tile_area(i,:)),i=1,np)/) &
						+ (/(sum(canopy%fnv(tile_index(i,:))*tile_area(i,:)),i=1,np)/) 
            ZphiEh(:,k) = (/(sum(canopy%fe(tile_index(i,:))*tile_area(i,:)),i=1,np)/)
		    !ZphiNEEh(:,k) = (/(sum((canopy%fpn(tile_index(i,:))+canopy%frp(tile_index(i,:)) &
			!                +canopy%frs(tile_index(i,:)))*tile_area(i,:)),i=1,np)/)*1.e6/12.
			ZphiNEEh(:,k) =	(/(sum((canopy%fpn(tile_index(i,:))-canopy%frday(tile_index(i,:))) &
		                   *tile_area(i,:)),i=1,np)/)*1.e6/12.
			ZphiGh(:,k) = (/(sum(canopy%ga(tile_index(i,:))*tile_area(i,:)),i=1,np)/)
			ZTsoilh(:,k) = ((ssoil%Tgg(1:mp_patch:ntile,2)-Tfrz)+(ssoil%Tgg(1:mp_patch:ntile,4)-Tfrz))/2. ! special for Otway

			! assume here equal soil temperature and moisture across tiles
			Tsoil1h(:,k) = ssoil%Tgg(1:mp_patch:ntile,1)-Tfrz
			Tsoil2h(:,k) = ssoil%Tgg(1:mp_patch:ntile,2)-Tfrz
			Tsoil3h(:,k) = ssoil%Tgg(1:mp_patch:ntile,3)-Tfrz
			Tsoil4h(:,k) = ssoil%Tgg(1:mp_patch:ntile,4)-Tfrz
			Tsoil5h(:,k) = ssoil%Tgg(1:mp_patch:ntile,5)-Tfrz
			Tsoil6h(:,k) = ssoil%Tgg(1:mp_patch:ntile,6)-Tfrz
			Tsoil7h(:,k) = ssoil%Tgg(1:mp_patch:ntile,7)-Tfrz
			Tsoil8h(:,k) = ssoil%Tgg(1:mp_patch:ntile,8)-Tfrz
			Tsoil9h(:,k) = ssoil%Tgg(1:mp_patch:ntile,9)-Tfrz
			Tsoil10h(:,k) = ssoil%Tgg(1:mp_patch:ntile,10)-Tfrz

			theta1h(:,k) = ssoil%wb(1:mp_patch:ntile,1)
			theta2h(:,k) = ssoil%wb(1:mp_patch:ntile,2)
			theta3h(:,k) = ssoil%wb(1:mp_patch:ntile,3)
			theta4h(:,k) = ssoil%wb(1:mp_patch:ntile,4)
			theta5h(:,k) = ssoil%wb(1:mp_patch:ntile,5)
			theta6h(:,k) = ssoil%wb(1:mp_patch:ntile,6)
			theta7h(:,k) = ssoil%wb(1:mp_patch:ntile,7)
			theta8h(:,k) = ssoil%wb(1:mp_patch:ntile,8)
			theta9h(:,k) = ssoil%wb(1:mp_patch:ntile,9)
			theta10h(:,k) = ssoil%wb(1:mp_patch:ntile,10)

	
		if (HourlyOutputFlag.eq.1) then
			CALL text_output(ktau, kstart, kend, dels,filename_out, &
                       met, canopy, ssoil,rad,bal,air,soil,veg, rough )

		elseif (HourlyOutputFlag.eq.2) then
			CALL write_output(ktau,dels,met,canopy,ssoil,rad,bal,air,soil,veg)
		endif
			
		! accumulate daily fluxes
		FWprec = FWprec + met%precip(1:mp_patch:ntile)/1000.                          ! precip (m d-1)
		FWthrough = FWthrough + (/(sum(canopy%through(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.                            ! throughfall (m d-1)
    	FWwc = FWwc + (/(sum(canopy%fevw(tile_index(i,:))/air%rlam(tile_index(i,:))*tile_area(i,:)),i=1,np)/)* &
		       dels /1000. 
	
		FWTra = FWTra + (/(sum(canopy%fevc(tile_index(i,:))/air%rlam(tile_index(i,:))*tile_area(i,:)),i=1,np)/)* &
                 dels/1000.


		
		FWTraA = FWTraA + (/(sum(ssoil%rexA(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000. 
		FWTraB = FWTraB + (/(sum(ssoil%rexB(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000. 
		if (ForwardModelFlag.eq.1) then
			FWsoil = FWsoil + (/(sum(canopy%fes(tile_index(i,:))/ssoil%cls(tile_index(i,:)) &
			                 *tile_area(i,:)),i=1,np)/)*dels/1000./air%Rlam                ! soil evaporation (md-1)
		elseif (ForwardModelFlag.eq.2) then
			FWsoil = FWsoil + (/(sum(ssoil%evap(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.                    ! soil evaporation (md-1)
		endif
		FWrun =  FWrun  + (/(sum(ssoil%rnof1(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.                                       ! surface runoff (m d-1)
		FWlchA = FWlchA + (/(sum(ssoil%leachAB(tile_index(i,:))*tile_area(i,:)),i=1,np)/)*dels/1000.                             ! water flux through bottom of 3rd layer (m d-1)
		FWlchB = FWlchB + (/(sum(ssoil%rnof2(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.										  ! deep drainaige
		FCGPP = FCGPP -(/(sum((canopy%fpn(tile_index(i,:))-canopy%frday(tile_index(i,:))) &
		                   *tile_area(i,:)),i=1,np)/)/12.0*dels	
		FCGro =  FCGro	-(/(sum((canopy%fpn(tile_index(i,:))+canopy%frp(tile_index(i,:)))* &
		                  tile_area(i,:)),i=1,np)/)/12.0*dels					  ! plant photosynthesis - plant respiration 	(mol C m-2 d-1)											!NPP
		FCNEE = FCNEE + ZphiNEEh(:,k)/1.e6*dels          ! NEE (plant respiration + soil respiration - photosynthesis) (mol C m-2 d-1)

			PhiE =   PhiE   + ZphiEh(:,k)										  ! latent heat flux (Wm-2)
			PhiH =   PhiH   + ZphiHh(:,k)										  ! sensible heat flux (Wm-2)
			ZPhiRnet = ZPhiRnet + ZphiRneth(:,k)*dels/1.e6                           ! net radiation flux (MJ m-2 d-1)
			ZPhiH = ZPhiH + ZphiHh(:,k)*dels/1.e6                                   ! sensible heat flux (MJ m-2 d-1)
			ZPhiE = ZPhiE + ZphiEh(:,k)*dels/1.e6                                   ! latent heat flux (MJ m-2 d-1)            ! daytime GPP (mol CO2 m-2 d-1)
			ZPhiNPP = ZPhiNPP - (/(sum((canopy%fpn(tile_index(i,:))+canopy%frp(tile_index(i,:)))* &
		                  tile_area(i,:)),i=1,np)/)/12.0*dels	            ! daytime NPP (mol CO2 m-2 d-1)
			ZPhiNEE = ZPhiNEE + + ZphiNEEh(:,k)/1.e6*dels              ! daytime NEE (mol CO2 m-2 d-1)

			ZPhiGPP=FCGPP


		FWPT = FWPT + ((/(sum(canopy%fns(tile_index(i,:))*tile_area(i,:)),i=1,np)/)+ &
		(/(sum(canopy%fnv(tile_index(i,:))*tile_area(i,:)),i=1,np)/)) &
									*dels/1000./air%Rlam &
		                            *(air%dsatdk/(air%dsatdk+air%psyc*0.01))*CoeffPT
		FWE = FWTra + FWsoil +FWwc	  ! total evaporation (m d-1)

		FWDis = FWrun + FWlchB                                                    ! total runoff (m d-1)
		dWCanopydt = dWCanopydt + (/(sum(canopy%delwc(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.                              ! change in canopy storage (m d-1)
		dWcoldt = dWcoldt + (/(sum(ssoil%delwcol(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.                                             ! change in stored water (m d-1)	
	    dWcolAdt = dWcolAdt + (/(sum(ssoil%delwcolA(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.
		dWcolBdt = dWcolBdt + (/(sum(ssoil%delwcolB(tile_index(i,:))*tile_area(i,:)),i=1,np)/)/1000.

		! DD varaibles for use in CASA-CNP
		Tabar         = Tabar + met%tk(tile_index(:,1)) 

		Tsoilbar_g    = Tsoilbar_g + sum(ssoil%tgg(tile_index(:,1),:)*veg%froot(tile_index(:,1),:),2)
		btran_g     = btran_g + sum(max((min(soil%sfc_vec(tile_index(:,1),:),ssoil%wb(tile_index(:,1),:)) &
		                          -soil%swilt_vec(tile_index(:,1),:)),0.0) &
								  /(soil%sfc_vec(tile_index(:,1),:)-soil%swilt_vec(tile_index(:,1),:)) &
								  *veg%froot(tile_index(:,1),:),2)
		FCGPP_g     = FCGPP_g  -(canopy%fpn(tile_index(:,1))-canopy%frday(tile_index(:,1)))/12.0*dels &
		                         *tile_area(:,1)	
		FCLeafR_g    = FCLeafR_g + canopy%frday(tile_index(:,1))/12.0*dels &
		                             *tile_area(:,1) 
		LAI_g = veg%vlai(tile_index(:,1))*tile_area(:,1) 
		fws_g = canopy%fwsoil(tile_index(:,1))


		if (CableParamFlag.eq.1) then  ! separate veg params for woody and grassy veg tiles
			Tsoilbar_w    = Tsoilbar_w + sum(ssoil%tgg(tile_index(:,2),:)*veg%froot(tile_index(:,2),:),2)
			
			btran_w    = btran_w + sum(max((min(soil%sfc_vec(tile_index(:,2),:),ssoil%wb(tile_index(:,2),:)) &
									  -soil%swilt_vec(tile_index(:,2),:)),0.0) &
									  /(soil%sfc_vec(tile_index(:,2),:)-soil%swilt_vec(tile_index(:,2),:)) &
									  *veg%froot(tile_index(:,2),:),2)
			FCGPP_w     = FCGPP_w -(canopy%fpn(tile_index(:,2))-canopy%frday(tile_index(:,2)))/12.0*dels &
			                       *tile_area(:,2)
			FCLeafR_w    = FCLeafR_w + canopy%frday(tile_index(:,2))/12.0*dels &
			                            *tile_area(:,2) 
			LAI_w = veg%vlai(tile_index(:,2))*tile_area(:,2)
			fws_w = canopy%fwsoil(tile_index(:,2))
		endif


	enddo   ! end call CABLE kend times

	! leaf NPP for comparison with VAST
	ZphiNPP = (FCGPP_g*0.63*AllocLg + FCGPP_w*0.46*AllocLw)*12.0*365.25*0.01 

	write(hmetn,"(4(i8,','),14(e15.6,','))") TTDATE%year,TTDATE%month, TTDATE%day,k-1,met%ca(1)*1e6
	
	! DD varaibles for use in CASA-CNP: take daily average
	Tabar = Tabar/real(ntime) -Tfrz
	Tsoilbar_g = Tsoilbar_g/real(ntime) - Tfrz
	btran_g     = btran_g/real(ntime)
	
	if (CableParamFlag.eq.1) then  ! separate veg params for woody and grassy veg tiles
		Tsoilbar_w = Tsoilbar_w/real(ntime) - Tfrz
		btran_w     = btran_w/real(ntime)
	endif

	Sbar_g      =  sum(ssoil%S(tile_index(:,1),1:3)*SPREAD(soil%zse(1:3),1,np),2) &
		                        /sum(SPREAD(soil%zse(1:3),1,np),2)

	
	Sbar_w      =  sum(ssoil%S(tile_index(:,2),1:6)*SPREAD(soil%zse(1:6),1,np),2) &
		                        /sum(SPREAD(soil%zse(1:6),1,np),2)

	
	PhiE = PhiE/real(kend)							! average latent and sensible heat fluxes over daylight hours
	PhiH = PhiH/real(kend)


	! soil-moisture observables
	Zsm0008 = (ssoil%wb(tile_index(:,1),1)*soil%zse(1) + ssoil%wb(tile_index(:,1),2)* &
	           soil%zse(2))/(soil%zse(1)+soil%zse(2))
	Zsm0090 = (ssoil%wb( tile_index(:,1),1)*soil%zse(1) + ssoil%wb(tile_index(:,1),2)*soil%zse(2)+ssoil%wb(tile_index(:,1),3)*soil%zse(3) &
			+ ssoil%wb(tile_index(:,1),4)*soil%zse(4)+ ssoil%wb(tile_index(:,1),5)*soil%zse(5)+ ssoil%wb(tile_index(:,1),6)*soil%zse(6))/ &
			(soil%zse(1)+soil%zse(2)+soil%zse(3)+ &
			soil%zse(4)+soil%zse(5)+soil%zse(6))

	Zsm0030 = (ssoil%wb(tile_index(:,1),1)*soil%zse(1) + ssoil%wb(tile_index(:,1),2)*soil%zse(2) &
	           +ssoil%wb(tile_index(:,1),3)*soil%zse(3) &
	          +ssoil%wb(tile_index(:,1),4)*soil%zse(4))/ &
			(soil%zse(1)+soil%zse(2)+soil%zse(3)+soil%zse(4))

	Zsm3060 = ssoil%wb(tile_index(:,1),5)

	Zsm6090 = ssoil%wb(tile_index(:,1),6)

	! special for Tumba
	where ((rad%latitude(tile_index(:,1)).gt.(-35.65-0.01).and.rad%latitude(tile_index(:,1)).lt.(-35.65+0.01)).and.(rad%longitude(tile_index(:,1)).gt.(148.15-0.01).and.rad%longitude(tile_index(:,1)).lt.(148.15+0.01))) 

		Zsm0008 =	(ssoil%wb(tile_index(:,1),1)*soil%zse(1) + ssoil%wb(tile_index(:,1),2)*soil%zse(2)+ssoil%wb(tile_index(:,1),3)*soil%zse(3))/&
			(soil%zse(1)+soil%zse(2)+soil%zse(3))
	endwhere
	
	! update state variables
	WRel1 = ssoil%S(tile_index(:,1),1)
	WRel2 = ssoil%S(tile_index(:,1),2)
	WRel3 = ssoil%S(tile_index(:,1),3)
	WRel4 = ssoil%S(tile_index(:,1),4)
	WRel5 = ssoil%S(tile_index(:,1),5)
	WRel6 = ssoil%S(tile_index(:,1),6)
	WRel7 = ssoil%S(tile_index(:,1),7)
	WRel8 = ssoil%S(tile_index(:,1),8)
	WRel9 = ssoil%S(tile_index(:,1),9)
	WRel10 = ssoil%S(tile_index(:,1),10)
	theta1 = ssoil%wb(tile_index(:,1),1)
	theta2 = ssoil%wb(tile_index(:,1),2)
	theta3 = ssoil%wb(tile_index(:,1),3)
	theta4 = ssoil%wb(tile_index(:,1),4)
	theta5 = ssoil%wb(tile_index(:,1),5)
	theta6 = ssoil%wb(tile_index(:,1),6)
	theta7 = ssoil%wb(tile_index(:,1),7)
	theta8 = ssoil%wb(tile_index(:,1),8)
	theta9 = ssoil%wb(tile_index(:,1),9)
	theta10 = ssoil%wb(tile_index(:,1),10)
	wcol = ssoil%wbtot(tile_index(:,1))/1000.
	WRelA = ssoil%SA(tile_index(:,1))
	WRelB = ssoil%SB(tile_index(:,1))
	Tsoil1 = ssoil%Tgg(tile_index(:,1),1) -Tfrz
	Tsoil2 = ssoil%Tgg(tile_index(:,1),2)-Tfrz 
	Tsoil3 = ssoil%Tgg(tile_index(:,1),3)-Tfrz 
	Tsoil4 = ssoil%Tgg(tile_index(:,1),4)-Tfrz 
	Tsoil5 = ssoil%Tgg(tile_index(:,1),5)-Tfrz 
	Tsoil6 = ssoil%Tgg(tile_index(:,1),6)-Tfrz
    Tsoil7 = ssoil%Tgg(tile_index(:,1),7) -Tfrz
	Tsoil8 = ssoil%Tgg(tile_index(:,1),8) -Tfrz
	Tsoil9 = ssoil%Tgg(tile_index(:,1),9) -Tfrz
	Tsoil10 = ssoil%Tgg(tile_index(:,1),10)  -Tfrz
	Clea   = bgc%cplant(tile_index(:,1),1)/12.0
	cansto = canopy%cansto(tile_index(:,1))
	cplant1 = (/(sum(bgc%cplant(tile_index(i,:),1)*tile_area(i,:)),i=1,np)/)
	cplant2 = (/(sum(bgc%cplant(tile_index(i,:),2)*tile_area(i,:)),i=1,np)/)
	cplant3 = (/(sum(bgc%cplant(tile_index(i,:),3)*tile_area(i,:)),i=1,np)/)
	fws = (/(sum(canopy%fwsoil(tile_index(i,:))*tile_area(i,:)),i=1,np)/)

	! daily water balance on soil

	Imbalsoil = FWthrough - FWtra-FWsoil-FWDis -dwColdt


	! daily water balance on canopy
	Imbalcanopy = FWprec - FWthrough -dWCanopydt - FWwc

	Imbal = Imbalsoil+Imbalcanopy

	CALL WaterDynObs (TTime, XX0, MM, UU, VV, AAP1,AAPh1, DD, FF,   ZZP1,ZZPh1,ZZC1)
endif

if (ForwardModelFlag.eq.0) then
! * MODEL COMPUTATION FOR ONE TIME STEP
!   * Calculate fluxes and increment XX
	CALL WaterDynFlux (TTime, XX0, MM,RR, UU, VV,   FF, DD)
	!XX1 = XX0 + DelT*FF(:,1:nxx)        ! FF(:,1:nxx) = first nxx entries of FF(np,nff)
	XX1(:,1:3) = XX0(:,1:3) + 1.0*FF(:,1:3)        ! FF(:,1:3) = first 3 entries of FF(np,nff)

	!   * Diagnosis for negative XX1 values
	if (DiagsFlag == 1) then
	  do ip=1,np
		do ixx=1,3
		  if (XX1(ip,ixx) <= 0.0) then
			write(*,"('Negative XX1(ip,ixx), ip, ixx:',e12.4,5x,2i8)") XX1(ip,ixx), ip, ixx
		  end if
		end do
	  end do
	end if

	!   * Calculate predicted observations (ZZP1, ZZC1)
	!     Use XX at either start of step (XX0) or end of step (XX1) %%%%
	!     Use forcing (MM) averaged through step
	! CALL WaterDynObs (TTime, XX1, MM, UU, VV, AAP1,   ZZP1, ZZC1)    ! XX at END of step
	CALL WaterDynObs (TTime, XX0, MM, UU, VV, AAP1,AAPh1, DD, FF,   ZZP1,ZZPh1, ZZC1)    ! XX at START of step
endif
! * TIME SERIES OUTPUT of spatial averages across whole domain, and ZZC, AAC for each catchment
!   * Find spatial statistics of np-arrays at current time step over whole domain
!     array dimensions: XX(np,nxx), XXSpStat(1,nxx,5) (5 for mean, sd, max, min, count)
CALL SpatialStats (MM,   MMSpStat)      ! spatial average over np points
CALL SpatialStats (RR,   RRSpStat)      ! spatial average over np points
CALL SpatialStats (XX1,  XXSpStat)      ! spatial average over np points
CALL SpatialStats (FF,   FFSpStat)      ! spatial average over np points
CALL SpatialStats (DD,   DDSpStat)      ! spatial average over np points
CALL SpatialStats (ZZP1, ZZPSpStat)     ! spatial average over np points
CALL SpatialStats (AAP1, AAPSpStat)     ! spatial average over np points
!   * Find cumulative time averages of np-spatial averages and catchment arrays
!     array dimensions for np-arrays: XXSpStat(1,nxx,5), XXSpTAv(1,nxx,nAvType)
!     array dimensions for nc-arrays: ZZC1(nc,nzzC), ZZCTAv(nc,nzzC,nAvType)
do AvType = 1,3                         ! AvType = (1,2,3) gives (month,ann,run averages)
  CALL TimeAverage (1, MMSpStat(:,:,1),  MMSpTAv(:,:,AvType),  MMSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, RRSpStat(:,:,1),  RRSpTAv(:,:,AvType),  RRSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, XXSpStat(:,:,1),  XXSpTAv(:,:,AvType),  XXSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, FFSpStat(:,:,1),  FFSpTAv(:,:,AvType),  FFSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, DDSpStat(:,:,1),  DDSpTAv(:,:,AvType),  DDSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, ZZPSpStat(:,:,1), ZZPSpTAv(:,:,AvType), ZZPSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, AAPSpStat(:,:,1), AAPSpTAv(:,:,AvType), AAPSpTAvCount(:,:,AvType))
  CALL TimeAverage (1, ZZC1, ZZCTAv(:,:,AvType), ZZCTAvCount(:,:,AvType))
  CALL TimeAverage (1, AAC1, AACTAv(:,:,AvType), AACTAvCount(:,:,AvType))
end do

  FFSpTAvTemp = FFSpTAv

!   * Write daily, monthly, annual, whole-of-run domain TS output (AvType = 0,1,2,3)
!       * DomTSOutputFlag = 0,1,2,3,4 -> no, run, run+ann, run+ann+mth, run+ann+mth+day

!     * daily domain time series output
if (DomTSOutputFlag == 4 .and. ipDiag == 0) then
  AvType = 0
  CALL TimeSeriesOutput                                                         &
         (1, AvType, ipDiag, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), TTime,   &
          MMSpStat(1,:,1),RRSpStat(1,:,1),XXSpStat(1,:,1), FFSpStat(1,:,1),     &
		   DDSpStat(1,:,1)   &
          ,ZZPSpStat(1,:,1), AAPSpStat(1,:,1))
end if

!     * monthly domain time series output and reinitialisation of accumulators
if (DomTSOutputFlag >= 3 .and. ipDiag == 0 .and. EndMonth(TTDate)) then
  AvType = 1
  CALL TimeSeriesOutput                                                         &
       (1, AvType, ipDiag, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), TTime,     &
        MMSpTAv(1,:,AvType),RRSpTAv(1,:,AvType), XXSpTAv(1,:,AvType),                               &
        FFSpTAv(1,:,AvType), DDSpTAv(1,:,AvType)                                &
        ,ZZPSpTAv(1,:,AvType), AAPSpTAv(1,:,AvType))
end if
if (EndMonth(TTDate)) then
  AvType = 1
  CALL TimeAverage (0, MMSpStat(:,:,1),  MMSpTAv(:,:,AvType),  MMSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, RRSpStat(:,:,1),  RRSpTAv(:,:,AvType),  RRSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, XXSpStat(:,:,1),  XXSpTAv(:,:,AvType),  XXSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, FFSpStat(:,:,1),  FFSpTAv(:,:,AvType),  FFSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, DDSpStat(:,:,1),  DDSpTAv(:,:,AvType),  DDSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, ZZPSpStat(:,:,1), ZZPSpTAv(:,:,AvType), ZZPSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, AAPSpStat(:,:,1), AAPSpTAv(:,:,AvType), AAPSpTAvCount(:,:,AvType))
end if


!     * annual domain time series output and reinitialisation of accumulators
if (DomTSOutputFlag >= 2 .and. ipDiag == 0 .and. EndYear(TTDate)) then
  AvType = 2
  CALL TimeSeriesOutput                                                         &
       (1, AvType, ipDiag, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), TTime,     &
        MMSpTAv(1,:,AvType),RRSpTAv(1,:,AvType), XXSpTAv(1,:,AvType),                               &
        FFSpTAv(1,:,AvType), DDSpTAv(1,:,AvType)                                &
        ,ZZPSpTAv(1,:,AvType), AAPSpTAv(1,:,AvType))
end if
if (EndYear(TTDate)) then
  AvType = 2
  CALL TimeAverage (0, MMSpStat(:,:,1),  MMSpTAv(:,:,AvType),  MMSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, RRSpStat(:,:,1),  RRSpTAv(:,:,AvType),  RRSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, XXSpStat(:,:,1),  XXSpTAv(:,:,AvType),  XXSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, FFSpStat(:,:,1),  FFSpTAv(:,:,AvType),  FFSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, DDSpStat(:,:,1),  DDSpTAv(:,:,AvType),  DDSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, ZZPSpStat(:,:,1), ZZPSpTAv(:,:,AvType), ZZPSpTAvCount(:,:,AvType))
  CALL TimeAverage (0, AAPSpStat(:,:,1), AAPSpTAv(:,:,AvType), AAPSpTAvCount(:,:,AvType))
end if

!     * whole-of-run domain time series output; reinitialise accumulators in InitModel
if (DomTSOutputFlag >= 1 .and. ipDiag == 0 .and. (TTDate .ge. EndDate)) then
  AvType = 3
  CALL TimeSeriesOutput                                                         &
       (1, AvType, ipDiag, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), TTime,     &
        MMSpTAv(1,:,AvType),RRSpTAv(1,:,AvType), XXSpTAv(1,:,AvType),                               &
        FFSpTAv(1,:,AvType), DDSpTAv(1,:,AvType)                                &
        ,ZZPSpTAv(1,:,AvType), AAPSpTAv(1,:,AvType))
end if

!   * Write daily, monthly, annual, whole-run catchment time series
!       * CatchTSOutputFlag = 0,1,2,3,4 -> no, run, run+ann, run+ann+mth, run+ann+mth+day

!     * daily catchment TS output and reinitialisation of accumulators
if (CatchTSOutputFlag == 4) then
  AvType = 0
  CALL CatchmentTSOutput                                                        &
         (1, AvType, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), iCatchID,        &
          1, 2, TTime, ZZC1, AAC1)
end if

!     * monthly catchment TS output and reinitialisation of accumulators
if (CatchTSOutputFlag >= 3 .and. EndMonth(TTDate)) then
  AvType = 1
  CALL CatchmentTSOutput                                                        &
       (1, AvType, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), iCatchID,          &
        1, 2, TTime, ZZCTAv(:,:,AvType), AACTAv(:,:,AvType))
end if
if (EndMonth(TTDate)) then
  AvType = 1
  CALL TimeAverage (0, ZZC1, ZZCTAv(:,:,AvType), ZZCTAvCount(:,:,AvType))
  CALL TimeAverage (0, AAC1, AACTAv(:,:,AvType), AACTAvCount(:,:,AvType))
  end if

!     * annual catchment TS output and reinitialisation of accumulators
if (CatchTSOutputFlag >= 2 .and. EndYear(TTDate)) then
  AvType = 2
  CALL CatchmentTSOutput                                                        &
       (1, AvType, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), iCatchID,          &
        1, 2, TTime, ZZCTAv(:,:,AvType), AACTAv(:,:,AvType))
end if
if (EndYear(TTDate)) then
  AvType = 2
  CALL TimeAverage (0, ZZC1, ZZCTAv(:,:,AvType), ZZCTAvCount(:,:,AvType))
  CALL TimeAverage (0, AAC1, AACTAv(:,:,AvType), AACTAvCount(:,:,AvType))
  end if

!     * whole-of-run catchment averages; reinitialise accumulators in InitModel
if (CatchTSOutputFlag >= 1 .and. (TTDate .ge. EndDate)) then
  AvType = 3
  CALL CatchmentRunAvgOutput (UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1),         &
       iCatchID, 1, 2, ZZCTAv(:,:,AvType), AACTAv(:,:,AvType))
end if

! * SPATIAL (MAP) OUTPUT and POINT TIME SERIES (DIAGNOSTIC) OUTPUT
!     * These are done together becasue they need the same accumulators
!     * MapOutputFlag = 0,1,2,3,4 -> no, run, run+ann, run+ann+mth, run+ann+mth+day

!   * Find cumulative time averages of np-maps at each point
!     * Array dimensions: XX(np,nxx), XXTAv(np,nxx,nAvType)
do AvType = 1,3                     ! AvType = (1,2,3) gives (monthly, annual, whole run)
                                    ! arrays: XX(np,nxx), XXTAv(np,nxx,nAvType)
  CALL TimeAverage (1, MM,   MMTAv(:,:,AvType),  MMTAvCount(:,:,AvType))
  CALL TimeAverage (1, RR,   RRTAv(:,:,AvType),  RRTAvCount(:,:,AvType))
  CALL TimeAverage (1, XX1,  XXTAv(:,:,AvType),  XXTAvCount(:,:,AvType))
  CALL TimeAverage (1, FF,   FFTAv(:,:,AvType),  FFTAvCount(:,:,AvType))
  CALL TimeAverage (1, DD,   DDTAv(:,:,AvType),  DDTAvCount(:,:,AvType))
  CALL TimeAverage (1, ZZP1, ZZPTAv(:,:,AvType), ZZPTAvCount(:,:,AvType))
  CALL TimeAverage (1, AAP1, AAPTAv(:,:,AvType), AAPTAvCount(:,:,AvType))
  CALL TimeAverageh (1, ZZPh1, ZZPhTAv(:,:,:,AvType), ZZPhTAvCount(:,:,:,AvType))
  CALL TimeAverageh (1, AAPh1, AAPhTAv(:,:,:,AvType), AAPhTAvCount(:,:,:,AvType))
end do

! overwrite AphiGPP, ZPhiGPP with sum over monthly mean diurnal cycle of NEE, with base (night-time) value subtracted
if (EndMonth(TTDate)) then 
	where((maxval(AAPhTAvCount(:,9,:,1),2).eq.0).or.(AAPhTAvCount(:,9,1,1).eq.0)) 
    AAPTAv(:,7,1) = ErrVal
	ZZPTAv(:,7,1) = ErrVal
	elsewhere	
		AAPTAv(:,7,1) = 0
		ZZPTAv(:,7,1) = 0
		AAPTAv(:,7,1) =- SUM((AAPhTAv(:,9,:,1)-SPREAD(AAPhTAv(:,9,1,1),2,ntime)),2,AAPhTAv(:,9,:,1)>ErrVal)/1.e6*dels
		ZZPTAv(:,7,1) = -SUM((ZZPhTAv(:,9,:,1)-SPREAD(ZZPhTAv(:,9,1,1),2,ntime)),2,AAPhTAv(:,9,:,1)>ErrVal)/1.e6*dels
	endwhere

! overwrite AphiE,ZphiE,AphiH,ZphiH,AphiRnet,ZphiRnet with sum over monthly mean diurnal cycle 

	where(maxval(AAPhTAvCount(:,8,:,1),2).eq.0) 
		AAPTAv(:,4,1) = ErrVal
		ZZPTAv(:,4,1) = ErrVal
	elsewhere	
		AAPTAv(:,4,1) = 0
		ZZPTAv(:,4,1) = 0
		AAPTAv(:,4,1) = SUM(AAPhTAv(:,8,:,1),2,AAPhTAv(:,8,:,1)>ErrVal)/1.e6*dels
		ZZPTAv(:,4,1) = SUM(ZZPhTAv(:,8,:,1),2,AAPhTAv(:,8,:,1)>ErrVal)/1.e6*dels
	endwhere

	where(maxval(AAPhTAvCount(:,7,:,1),2).eq.0) 
		AAPTAv(:,3,1) = ErrVal
		ZZPTAv(:,3,1) = ErrVal
	elsewhere	
		AAPTAv(:,3,1) = 0
		ZZPTAv(:,3,1) = 0
		AAPTAv(:,3,1) = SUM(AAPhTAv(:,7,:,1),2,AAPhTAv(:,7,:,1)>ErrVal)/1.e6*dels
		ZZPTAv(:,3,1) = SUM(ZZPhTAv(:,7,:,1),2,AAPhTAv(:,7,:,1)>ErrVal)/1.e6*dels
	endwhere

	where(maxval(AAPhTAvCount(:,6,:,1),2).eq.0)
		AAPTAv(:,2,1) = ErrVal
		ZZPTAv(:,2,1) = ErrVal
	elsewhere	
		AAPTAv(:,2,1) = 0
		ZZPTAv(:,2,1) = 0
		AAPTAv(:,2,1) = SUM(AAPhTAv(:,6,:,1),2,AAPhTAv(:,6,:,1)>ErrVal)/1.e6*dels
		ZZPTAv(:,2,1) = SUM(ZZPhTAv(:,6,:,1),2,AAPhTAv(:,6,:,1)>ErrVal)/1.e6*dels
	endwhere
endif
if (EndYear(TTDate)) then  ! annual
	where((maxval(AAPhTAvCount(:,9,:,2),2).eq.0).or.(AAPhTAvCount(:,9,1,2).eq.0)) 
    AAPTAv(:,7,2) = ErrVal
	ZZPTAv(:,7,2) = ErrVal
	elsewhere	
		AAPTAv(:,7,2) = 0
		ZZPTAv(:,7,2) = 0
		AAPTAv(:,7,2) =- SUM((AAPhTAv(:,9,:,2)-SPREAD(AAPhTAv(:,9,1,2),2,ntime)),2,AAPhTAv(:,9,:,2)>ErrVal)/1.e6*dels
		ZZPTAv(:,7,2) = -SUM((ZZPhTAv(:,9,:,2)-SPREAD(ZZPhTAv(:,9,1,2),2,ntime)),2,AAPhTAv(:,9,:,2)>ErrVal)/1.e6*dels
	endwhere

! overwrite AphiE,ZphiE,AphiH,ZphiH,AphiRnet,ZphiRnet with sum over monthly mean diurnal cycle 

	where(maxval(AAPhTAvCount(:,8,:,2),2).eq.0) 
		AAPTAv(:,4,2) = ErrVal
		ZZPTAv(:,4,2) = ErrVal
	elsewhere	
		AAPTAv(:,4,2) = 0
		ZZPTAv(:,4,2) = 0
		AAPTAv(:,4,2) = SUM(AAPhTAv(:,8,:,2),2,AAPhTAv(:,8,:,2)>ErrVal)/1.e6*dels
		ZZPTAv(:,4,2) = SUM(ZZPhTAv(:,8,:,2),2,AAPhTAv(:,8,:,2)>ErrVal)/1.e6*dels
	endwhere

	where(maxval(AAPhTAvCount(:,7,:,2),2).eq.0) 
		AAPTAv(:,3,2) = ErrVal
		ZZPTAv(:,3,2) = ErrVal
	elsewhere	
		AAPTAv(:,3,2) = 0
		ZZPTAv(:,3,2) = 0
		AAPTAv(:,3,2) = SUM(AAPhTAv(:,7,:,2),2,AAPhTAv(:,7,:,2)>ErrVal)/1.e6*dels
		ZZPTAv(:,3,2) = SUM(ZZPhTAv(:,7,:,2),2,AAPhTAv(:,7,:,2)>ErrVal)/1.e6*dels
	endwhere

	where(maxval(AAPhTAvCount(:,6,:,2),2).eq.0)
		AAPTAv(:,2,2) = ErrVal
		ZZPTAv(:,2,2) = ErrVal
	elsewhere	
		AAPTAv(:,2,2) = 0
		ZZPTAv(:,2,2) = 0
		AAPTAv(:,2,2) = SUM(AAPhTAv(:,6,:,2),2,AAPhTAv(:,6,:,2)>ErrVal)/1.e6*dels
		ZZPTAv(:,2,2) = SUM(ZZPhTAv(:,6,:,2),2,AAPhTAv(:,6,:,2)>ErrVal)/1.e6*dels
	endwhere


endif ! end overwrite AphiGPP, ZPhiGPP, AphiE,ZphiE,AphiH,ZphiH,AphiRnet,ZphiRnet with sum over mean diurnal cycle 

! * PEST OUTPUT 
if (PESTOutputFlag >= 2) then
!   * Write observations or observables for PEST run
  CALL WritePESTOutput (Time0,TTDate,AAP1,AAC1,ZZP1,ZZC1,FFTAv,VV,1,AAPH1,ZZPH1,AAPTAv,ZZPTAv)
!   * Close PEST files if end of run
!  if (TTime .eq. EndDate) CALL WritePESTOutput (Time0,TTDate,AAP1,AAC1,ZZP1,ZZC1,2)
  if (Time0 .eq. nSteps*1.0) CALL WritePESTOutput (Time0,TTDate,AAP1,AAC1,ZZP1,ZZC1,FFTAv,VV,2,AAPH1,ZZPH1,AAPTAv,ZZPTAv)
end if


!     * daily map output, point TS output
if (MapOutputFlag == 4) then        ! daily point map OP disabled as a safety measure
  CALL MapOutput (MM(:,:),  MMname(:,1),  MMMapOPFlag,  'day', TTime)
  CALL MapOutput (RR(:,:),  RRname(:,1),  RRMapOPFlag,  'day', TTime)
  CALL MapOutput (XX1(:,:), XXname(:,1),  XXMapOPFlag,  'day', TTime)
  CALL MapOutput (FF(:,:),  FFname(:,1),  FFMapOPFlag,  'day', TTime)
  CALL MapOutput (DD(:,:),  DDname(:,1),  DDMapOPFlag,  'day', TTime)
  CALL MapOutput (ZZP1(:,:), ZZPname(:,1), ZZPMapOPFlag, 'day', TTime)
  CALL MapOutput (AAP1(:,:), AAPname(:,1), AAPMapOPFlag, 'day', TTime)
  CALL MapOutputh (ZZPh1(:,:,:), ZZPhname(:,1), ZZPhMapOPFlag,Ihour_output,'day', TTime)
  CALL MapOutputh (AAPh1(:,:,:), AAPhname(:,1), AAPhMapOPFlag,Ihour_output,'day', TTime)
end if
!     * daily point TS output
if (ipDiag > 0) then
  AvType = 0
 CALL TimeSeriesOutput (1, AvType, ipDiag, UU, VV(ipDiag,:), XX0(ipDiag,:),    &
          TTime, MM(ipDiag,:),RR(ipDiag,:), XX1(ipDiag,:), FF(ipDiag,:), DD(ipDiag,:)       &
          ,ZZP1(ipDiag,:), AAP1(ipDiag,:))
end if

! * DD output (spatial TS)
if (MapTSOutputFlag == 4) then ! DD TS output
	AvType = 0
	CALL TimeSeriesOutput_bin (1, DDTSOPFlag,AvType,'day', TTime,DD,DDName(:,1),1)
	CALL TimeSeriesOutput_bin (1, XXTSOPFlag,AvType,'day', TTime,XX1,XXName(:,1),2)
	CALL TimeSeriesOutput_bin (1, FFTSOPFlag,AvType,'day', TTime,FF,FFName(:,1),3)
endif
if (MapTSOutputFlag >= 3 .and. EndMonth(TTDate)) then 
	AvType = 1
	CALL TimeSeriesOutput_bin (1, DDTSOPFlag,AvType,'mth', TTime,DDTAv(:,:,AvType),DDName(:,1),1)
	CALL TimeSeriesOutput_bin (1, XXTSOPFlag,AvType,'mth', TTime,XXTAv(:,:,AvType),XXName(:,1),2)
	CALL TimeSeriesOutput_bin (1, FFTSOPFlag,AvType,'mth', TTime,FFTAv(:,:,AvType),FFName(:,1),3)
endif
if (MapTSOutputFlag >= 2 .and. EndYear(TTDate)) then          ! TS OP
   AvType = 2
   CALL TimeSeriesOutput_bin (1, DDTSOPFlag,AvType,'ann', TTime,DDTAv(:,:,AvType),DDName(:,1),1)
   CALL TimeSeriesOutput_bin (1, XXTSOPFlag,AvType,'ann', TTime,XXTAv(:,:,AvType),XXName(:,1),2)
   CALL TimeSeriesOutput_bin (1, FFTSOPFlag,AvType,'ann', TTime,FFTAv(:,:,AvType),FFName(:,1),3)
endif

if (MapTSOutputFlag >= 1 .and. (TTDate .ge. EndDate)) then          ! TS OP
   AvType = 3
   CALL TimeSeriesOutput_bin (1, DDTSOPFlag,AvType,'run', TTime,DDTAv(:,:,AvType),DDName(:,1),1)
   CALL TimeSeriesOutput_bin (1, XXTSOPFlag,AvType,'run', TTime,XXTAv(:,:,AvType),XXName(:,1),2)
   CALL TimeSeriesOutput_bin (1, FFTSOPFlag,AvType,'run', TTime,FFTAv(:,:,AvType),FFName(:,1),3)
endif


!     * monthly map output, point TS output, and reinitialisation of accumulators
if (MapOutputFlag >= 3 .and. EndMonth(TTDate)) then         ! Map OP
  AvType = 1
  CALL MapOutput (MMTAv(:,:,AvType),  MMname(:,1),  MMMapOPFlag,  'mth', TTime)
  CALL MapOutput (RRTAv(:,:,AvType),  RRname(:,1),  RRMapOPFlag,  'mth', TTime)
  CALL MapOutput (XXTAv(:,:,AvType),  XXname(:,1),  XXMapOPFlag,  'mth', TTime)
  CALL MapOutput (FFTAv(:,:,AvType),  FFname(:,1),  FFMapOPFlag,  'mth', TTime)
  CALL MapOutput (DDTAv(:,:,AvType),  DDname(:,1),  DDMapOPFlag,  'mth', TTime)
  CALL MapOutput (ZZPTAv(:,:,AvType), ZZPname(:,1), ZZPMapOPFlag, 'mth', TTime)
  CALL MapOutput (AAPTAv(:,:,AvType), AAPname(:,1), AAPMapOPFlag, 'mth', TTime)
  CALL MapOutputh (ZZPhTAv(:,:,:,AvType), ZZPhname(:,1), ZZPhMapOPFlag,Ihour_output,'mth', TTime)
  CALL MapOutputh (AAPhTAv(:,:,:,AvType), AAPhname(:,1), AAPhMapOPFlag,Ihour_output,'mth', TTime)
end if
if (ipDiag > 0 .and. EndMonth(TTDate)) then                 ! Point TS OP
  AvType = 1
  CALL TimeSeriesOutput (1, AvType, ipDiag, UU, VV(ipDiag,:), XX0(ipDiag,:),    &
          TTime, MMTAv(ipDiag,:,AvType),RRTAv(ipDiag,:,AvType), XXTAv(ipDiag,:,AvType),                & 
          FFTAv(ipDiag,:,AvType), DDTAv(ipDiag,:,AvType)                        &
          ,ZZPTAv(ipDiag,:,AvType), AAPTAv(ipDiag,:,AvType))
end if
if (EndMonth(TTDate)) then                                  ! reinitialise accumulators
  AvType = 1
  CALL TimeAverage (0, MM,   MMTAv(:,:,AvType),  MMTAvCount(:,:,AvType))
  CALL TimeAverage (0, RR,   RRTAv(:,:,AvType),  RRTAvCount(:,:,AvType))
  CALL TimeAverage (0, XX1,  XXTAv(:,:,AvType),  XXTAvCount(:,:,AvType))
  CALL TimeAverage (0, FF,   FFTAv(:,:,AvType),  FFTAvCount(:,:,AvType))
  CALL TimeAverage (0, DD,   DDTAv(:,:,AvType),  DDTAvCount(:,:,AvType))
  CALL TimeAverage (0, ZZP1, ZZPTAv(:,:,AvType), ZZPTAvCount(:,:,AvType))
  CALL TimeAverage (0, AAP1, AAPTAv(:,:,AvType), AAPTAvCount(:,:,AvType))
  CALL TimeAverageh (0, ZZPh1, ZZPhTAv(:,:,:,AvType), ZZPhTAvCount(:,:,:,AvType))
  CALL TimeAverageh (0, AAPh1, AAPhTAv(:,:,:,AvType), AAPhTAvCount(:,:,:,AvType))
end if

!     * annual map output, point TS output, and reinitialisation of accumulators
if (MapOutputFlag >= 2 .and. EndYear(TTDate)) then          ! Map OP
  AvType = 2
  CALL MapOutput (MMTAv(:,:,AvType),  MMname(:,1),  MMMapOPFlag,  'ann', TTime)
  CALL MapOutput (RRTAv(:,:,AvType),  RRname(:,1),  RRMapOPFlag,  'ann', TTime)
  CALL MapOutput (XXTAv(:,:,AvType),  XXname(:,1),  XXMapOPFlag,  'ann', TTime)
  CALL MapOutput (FFTAv(:,:,AvType),  FFname(:,1),  FFMapOPFlag,  'ann', TTime)
  CALL MapOutput (DDTAv(:,:,AvType),  DDname(:,1),  DDMapOPFlag,  'ann', TTime)
  CALL MapOutput (ZZPTAv(:,:,AvType), ZZPname(:,1), ZZPMapOPFlag, 'ann', TTime)
  CALL MapOutput (AAPTAv(:,:,AvType), AAPname(:,1), AAPMapOPFlag, 'ann', TTime)
  CALL MapOutputh (ZZPhTAv(:,:,:,AvType), ZZPhname(:,1), ZZPhMapOPFlag,Ihour_output,'ann', TTime)
  CALL MapOutputh (AAPhTAv(:,:,:,AvType), AAPhname(:,1), AAPhMapOPFlag,Ihour_output,'ann', TTime)
end if
if (ipDiag > 0 .and. EndYear(TTDate)) then                  ! Point TS OP
  AvType = 2
  CALL TimeSeriesOutput (1, AvType, ipDiag, UU, VV(ipDiag,:), XX0(ipDiag,:),    &
          TTime, MMTAv(ipDiag,:,AvType), RRTAv(ipDiag,:,AvType),XXTAv(ipDiag,:,AvType),                & 
          FFTAv(ipDiag,:,AvType), DDTAv(ipDiag,:,AvType)                        &
          ,ZZPTAv(ipDiag,:,AvType), AAPTAv(ipDiag,:,AvType))
end if
if (EndYear(TTDate)) then                                   ! reinitialise accumulators
  AvType = 2
  CALL TimeAverage (0, MM,   MMTAv(:,:,AvType),  MMTAvCount(:,:,AvType))
  CALL TimeAverage (0, RR,   RRTAv(:,:,AvType),  RRTAvCount(:,:,AvType))
  CALL TimeAverage (0, XX1,  XXTAv(:,:,AvType),  XXTAvCount(:,:,AvType))
  CALL TimeAverage (0, FF,   FFTAv(:,:,AvType),  FFTAvCount(:,:,AvType))
  CALL TimeAverage (0, DD,   DDTAv(:,:,AvType),  DDTAvCount(:,:,AvType))
  CALL TimeAverage (0, ZZP1, ZZPTAv(:,:,AvType), ZZPTAvCount(:,:,AvType))
  CALL TimeAverage (0, AAP1, AAPTAv(:,:,AvType), AAPTAvCount(:,:,AvType))
  CALL TimeAverageh (0, ZZPh1, ZZPhTAv(:,:,:,AvType), ZZPhTAvCount(:,:,:,AvType))
  CALL TimeAverageh (0, AAPh1, AAPhTAv(:,:,:,AvType), AAPhTAvCount(:,:,:,AvType))
end if

!     * whole-of-run map output, point TS output; reinitialise accumulators in InitModel
if (MapOutputFlag >= 1 .and. (TTDate .ge. EndDate)) then   ! Map OP
  AvType = 3
  CALL MapOutput (MMTAv(:,:,AvType),  MMname(:,1),  MMMapOPFlag,  'run', TTime)
  CALL MapOutput (RRTAv(:,:,AvType),  RRname(:,1),  RRMapOPFlag,  'run', TTime)
  CALL MapOutput (XXTAv(:,:,AvType),  XXname(:,1),  XXMapOPFlag,  'run', TTime)
  CALL MapOutput (FFTAv(:,:,AvType),  FFname(:,1),  FFMapOPFlag,  'run', TTime)
  CALL MapOutput (DDTAv(:,:,AvType),  DDname(:,1),  DDMapOPFlag,  'run', TTime)
  CALL MapOutput (ZZPTAv(:,:,AvType), ZZPname(:,1), ZZPMapOPFlag, 'run', TTime)
  CALL MapOutput (AAPTAv(:,:,AvType), AAPname(:,1), AAPMapOPFlag, 'run', TTime)
  CALL MapOutputh (ZZPhTAv(:,:,:,AvType), ZZPhname(:,1), ZZPhMapOPFlag,Ihour_output,'run', TTime)
  CALL MapOutputh (AAPhTAv(:,:,:,AvType), AAPhname(:,1), AAPhMapOPFlag,Ihour_output,'run', TTime)
end if
if (ipDiag > 0 .and. (TTDate .ge. EndDate)) then           ! Point TS OP
  AvType = 3
  CALL TimeSeriesOutput (1, AvType, ipDiag, UU, VV(ipDiag,:), XX0(ipDiag,:),    &
          TTime, MMTAv(ipDiag,:,AvType),RRTAv(ipDiag,:,AvType), &
		   XXTAv(ipDiag,:,AvType),                & 
          FFTAv(ipDiag,:,AvType), DDTAv(ipDiag,:,AvType)                        &
          ,ZZPTAv(ipDiag,:,AvType), AAPTAv(ipDiag,:,AvType))
end if

! * Write monthly XX1(np,nxx) to supply temporary initial stores for next run
if (WriteInitStoresFlag >= 1 .and. EndMonth(TTDate)) then
  CALL MapOutput_tmp (XX1(:,:), XXname(:,1),  XXMapOPFlag,  'ini', TTime)
end if


! * Write last XX1(np,nxx) to supply initial stores for next run
if (WriteInitStoresFlag >= 1 .and. (TTDate .ge. EndDate)) then
  CALL MapOutput (XX1(:,:), XXname(:,1),  XXMapOPFlag,  'ini', TTime)
end if


Time0Last = Time0

if (DiagsFlag == 1 .and. EndMonth(TTDate)) then
  CALL CPU_TIME (cpuTimeNow)
  write(*,"('BIOS2Step: completed',3i6,'  at cpuTime (s) =',f12.2)")         & 
    nint(TTime(2:4)), cpuTimeNow - cpuTimeBegin
end if

if (Time0 .eq. nSteps*1.0.and. HourlyOutputFlag.eq.2) then
CALL close_output_file(bal, air, bgc, canopy, met, rad, &
       rough, soil, ssoil, sum_flux, veg, model_structure)
endif

CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE InitModel
!-------------------------------------------------------------------------------
! Initialise quantities with local scope to WaterDynStep:
! * Allocate arrays with dimension np or nc 
! * Initialise spatial arrays
! * Called from WaterDynStep on first call only
!-------------------------------------------------------------------------------
USE define_dimensions, ONLY:r_1,r_2,i_d,ms,mp_patch,mp           ! for CABLE
implicit none
! Local variables (only needed here)
integer(i4b)  :: ivv, ixx, imm,irr, ic, i, iunit,istat
real(sp)      :: VVtest                  ! test element for reading VV
integer(i4b)  :: nVVtest(nvv)            ! test number of elements in VV data
logical(lgt)  :: nVVtestFlag(nvv)        ! flag for nVVtest = np
character(6)  :: CatchName
integer(i4b)  :: iday, iSkipDays
type(dmydate) :: FirstMetDate, MetDate,RRDate, DischargeDate
integer(i4b)  :: ObsDay, ObsMonth, ObsYear
character(1)  :: YN
character(100):: dumchar
real(sp)      :: dum(20), CatchIDdum
integer(i4b)  :: EndMetFlag
logical(lgt)   :: GenFileExistFlag
!-------------------------------------------------------------------------------

! * ARRAY ALLOCATION for arrays private to WaterDynStep
!if (.not.allocated(VV)) allocate (VV(np,nvv))           ! spatially variable params
if (.not.allocated(FF)) allocate (FF(np,nff))           ! fluxes and dXXdt
if (.not.allocated(DD)) allocate (DD(np,ndd))           ! derived quantities
if (.not.allocated(MM)) allocate (MM(np,nmm))           ! met variables
if (.not.allocated(RR)) allocate (RR(np,nrr))           ! remote senisng driver variables
if (.not.allocated(MMprev)) allocate (MMprev(np,nmm))           ! met variables
if (.not.allocated(MMnext)) allocate (MMnext(np,nmm))           ! met variables
if (.not.allocated(hMM)) allocate (hMM(np,nhmm,ntime))           ! hourly met variables
if (.not.allocated(hMMnext)) allocate (hMMnext(np,nhmm,ntime))           ! hourly met variables
if (.not.allocated(hMMlocal)) allocate (hMMlocal(np,nhmm,ntime))           ! hourly met variables
if (.not.allocated(MMTAv))       allocate (MMTAv(np,nmm,nAvType))           ! accumulators
if (.not.allocated(MMTAvCount))  allocate (MMTAvCount(np,nmm,nAvType))      ! counters
if (.not.allocated(RRTAv))       allocate (RRTAv(np,nrr,nAvType))           ! accumulators
if (.not.allocated(RRTAvCount))  allocate (RRTAvCount(np,nrr,nAvType))      ! counters
if (.not.allocated(XXTAv))       allocate (XXTAv(np,nxx,nAvType))           ! accumulators
if (.not.allocated(XXTAvCount))  allocate (XXTAvCount(np,nxx,nAvType))      ! counters
if (.not.allocated(FFTAv))       allocate (FFTAv(np,nff,nAvType))           ! accumulators
if (.not.allocated(FFTAvCount))  allocate (FFTAvCount(np,nff,nAvType))      ! counters
if (.not.allocated(DDTAv))       allocate (DDTAv(np,ndd,nAvType))           ! accumulators
if (.not.allocated(DDTAvCount))  allocate (DDTAvCount(np,ndd,nAvType))      ! counters
if (.not.allocated(ZZPTAv))      allocate (ZZPTAv(np,nzzP,nAvType))         ! accumulators
if (.not.allocated(ZZPTAvCount)) allocate (ZZPTAvCount(np,nzzP,nAvType))    ! counters
if (.not.allocated(ZZPhTAv))      allocate (ZZPhTAv(np,nzzPh,ntime,nAvType))         ! accumulators
if (.not.allocated(ZZPhTAvCount)) allocate (ZZPhTAvCount(np,nzzPh,ntime,nAvType))    ! counters
if (.not.allocated(AAPTAv))      allocate (AAPTAv(np,naaP,nAvType))         ! accumulators
if (.not.allocated(AAPTAvCount)) allocate (AAPTAvCount(np,naaP,nAvType))    ! counters
if (.not.allocated(AAPhTAv))      allocate (AAPhTAv(np,naaPh,ntime,nAvType))         ! accumulators
if (.not.allocated(AAPhTAvCount)) allocate (AAPhTAvCount(np,naaPh,ntime,nAvType))    ! counters
if (.not.allocated(ZZCTAv))      allocate (ZZCTAv(nc,nzzC,nAvType))         ! accumulators
if (.not.allocated(ZZCTAvCount)) allocate (ZZCTAvCount(nc,nzzC,nAvType))    ! counters
if (.not.allocated(AACTAv))      allocate (AACTAv(nc,naaC,nAvType))         ! accumulators
if (.not.allocated(AACTAvCount)) allocate (AACTAvCount(nc,naaC,nAvType))    ! counters
if (.not.allocated(ZZPhSpTAv)) allocate (ZZPhSpTAv(1,nzzP,ntime,nAvType))
if (.not.allocated(ZZPhSpTAvCount)) allocate (ZZPhSpTAvCount(1,nzzP,ntime,nAvType))
if (.not.allocated(AAPhSpTAv)) allocate (AAPhSpTAv(1,nzzP,ntime,nAvType))
if (.not.allocated(AAPhSpTAvCount)) allocate (AAPhSpTAvCount(1,nzzP,ntime,nAvType))
if (.not.allocated(tile_area)) allocate (tile_area(np,ntile))
if (.not.allocated(tile_index)) allocate (tile_index(np,ntile))
! * Assign pointers to generic arrays
CALL PointAll (TargetUU=UU, TargetVV=VV, TargetAAP=AAP1,TargetAAPh=AAPh1, TargetZZC=ZZC1)

! * CONSISTENCY CHECKS on UU parameters
if ( .not. ( (alfaWpri > 0.0 .and. alfaWmul < 0.0) .or. &
             (alfaWpri < 0.0 .and. alfaWmul > 0.0) ) )  &
      STOP "InitModel: illegal alfaWpri, alfaWmul combination"

!   * Apply multiplier to ZSoil1, ZSoil2, and hardwired defaults to ZSoil and WVolSat
where (ZSoil1(:) > ErrVal+1.0)
  ZSoil1(:) = ZSoil1(:) * ZSoil1Mult
elsewhere
  ZSoil1(:) = 0.2
end where
where (ZSoil2(:) > ErrVal+1.0)
  ZSoil2(:) = ZSoil2(:) * ZSoil2Mult
elsewhere
  ZSoil2(:) = 0.2
end where
where (WVolSat1(:) < ErrVal+1.0) WVolSat1(:) = 0.4
where (WVolSat2(:) < ErrVal+1.0) WVolSat2(:) = 0.4


! * OPEN AND INITIALISE MET FILES
Time0Last = ErrVal                              ! initialise Time0Last
do imm = 1,nmm-2                                  ! go through met files, except vp (one per variable)
  GenFileName = trim(DataInRootDir) // trim(DomName) // '/met/' // trim(MetFileName(imm))
  iunit = 30 + imm
  CLOSE (iunit)  
  INQUIRE (file=GenFileName, exist=GenFileExistFlag)  ! ensure file closed before opening
  if (GenFileExistFlag) then
	  OPEN (iunit, file=GenFileName, form='binary', status='old', action='read')
	  read (iunit)  FirstMetDate, MM(:,imm)         ! read first date
	  REWIND (iunit)                                ! reset file at start
	  iSkipDays = DayDifference(StartDate, FirstMetDate)    ! StartDate-FirstMetDate (days)
	  if (iSkipDays < 0) then
		write(*,"('imm,StartDate,FirstMetDate:',7i6)") imm, StartDate, FirstMetDate
		write(*,"('MetFileName: ',a)") trim(MetFileName(imm))
		STOP 'InitModel: StartDate before FirstMetDate' 
	  end if
	  if (DiagsFlag == 1) then
		write(*,"('Initialise: ',a)") trim(GenFileName)
	  end if
	  do iday=1,iSkipDays                           ! skip until StartDate reached
		read (iunit, end=20) MetDate, MM(:,imm)
	  end do
  endif
end do

do imm = nmm-1,nmm                                  ! go vp (one per variable)
  GenFileName = trim(DataInRootDir) // trim(DomName) // '/met/' // trim(MetFileName(imm))
  iunit = 30 + imm
  CLOSE (iunit)  
  INQUIRE (file=GenFileName, exist=GenFileExistFlag)  ! ensure file closed before opening
  if (GenFileExistFlag) then
	  OPEN (iunit, file=GenFileName, form='binary', status='old', action='read')
	  read (iunit)  FirstMetDate, MM(:,imm)         ! read first date
	  REWIND (iunit)                                ! reset file at start
	  iSkipDays = DayDifference(StartDate, FirstMetDate)    ! StartDate-FirstMetDate (days)
	  if (iSkipDays < 0) then
		write(*,"('imm,StartDate,FirstMetDate:',7i6)") imm, StartDate, FirstMetDate
		write(*,"('MetFileName: ',a)") trim(MetFileName(imm))
		STOP 'InitModel: StartDate before FirstMetDate' 
	  end if
	  if (DiagsFlag == 1) then
		write(*,"('Initialise: ',a)") trim(GenFileName)
	  end if
	  do iday=1,iSkipDays                           ! skip until StartDate reached
		read (iunit, end=20) MetDate, MM(:,imm)
	  end do
  else
	  MM(:,imm)= esatf(MM(:,4))  ! set vp to saturated vp at Tmin
  end if
end do

goto 21

20 if (UseLocalMetDataFlag == 0) then
	stop 'InitMod: end of daily met file'
  else
    EndMetFlag=1
	continue
  end if

! * OPEN AND INITIALISE LOCAL MET FILES (optional)
21 if (UseLocalMetDataFlag>=1) then
  Time0Last = ErrVal                              ! initialise Time0Last
  GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(LocalMetFileName)
  iunit = 30 + nmm +1
  CLOSE (iunit)                                 ! ensure file closed before opening
  OPEN (iunit, file=GenFileName, form='binary', status='old', action='read')
  !read (iunit)  FirstMetDate, hMM(:,:,:)         ! read first date
  if (UseLocalMetDataFlag == 1) then
		read (iunit)  FirstMetDate, hMM(:,:,:)
  elseif (UseLocalMetDataFlag == 2) then
		read (iunit)  FirstMetDate, hMM(1,:,:)
		hMM = SPREAD(hMM(1,:,:),1,np)
  endif
  REWIND (iunit)                                ! reset file at start
  iSkipDays = DayDifference(StartDate, FirstMetDate)    ! StartDate-FirstMetDate (days)
  if (iSkipDays < 0) then
    write(*,"('StartDate,FirstLocalMetDate:',7i6)") imm, StartDate, FirstMetDate
    write(*,"('LocalMetFileName: ',a)") trim(LocalMetFileName)
	STOP 'InitModel: StartDate before FirstLocalMetDate' 
  end if
  if (DiagsFlag == 1) then
    write(*,"('Initialise: ',a)") trim(GenFileName)
  end if
  do iday=1,iSkipDays                           ! skip until StartDate reached
  if (UseLocalMetDataFlag == 1) then
    read (iunit) MetDate, hMM(:,:,:)
   else
	read (iunit) MetDate, hMM(1,:,:)
	hMM = SPREAD(hMM(1,:,:),1,np)
  endif

  end do
endif

! * OPEN and INITIALISE atmospheric CO2 file
if (InitFlagCO2==0) then
NextCO2Date%day = 99; NextCO2Date%month = 99; NextCO2Date%year = 9999
GenFileName = trim(DataInRootDir)   &
                 // trim(CO2FileName)
	INQUIRE (file=GenFileName, exist=GenFileExistFlag)
	if (DiagsFlag == 1) then
		if (.not.GenFileExistFlag) write(*,"('File not found:',/,a)") GenFileName
	endif
	if (GenFileExistFlag) then
		iunit =90
		 OPEN (iunit, file=GenFileName, status='old', form='binary', action='read')
		 REWIND (iunit)
		 read(iunit) NextCO2Date, NextCO2
		 write(*,*) NextCO2Date, NextCO2
		 do while(.not.EqMonths(NextCO2Date,StartDate_run)) 
			read(iunit,IOSTAT=istat) NextCO2Date, NextCO2
			if (istat.eq.-1) then
				write(*,*) 'InitModel: StartDate before FirstCO2Date' 
				NextCO2 = CO2A*1.e6
				exit ! end of file
			endif
		 enddo
		 CO2A = NextCO2/1.e6
	endif

	if (DiagsFlag == 1) then
		write(*,"('Initial CO2 ppm: ', f8.2)") CO2A*1e6
		write(*,"('CO2 start date: ', 3i5)") NextCO2Date
	endif
endif



! * OPEN and INITIALISE REMOTE SENSING DRIVER FILES
if (.not.allocated(NextRR)) allocate (NextRR(np,nrr))
if (.not.allocated(RRtmp)) allocate (RRtmp(np,nrr))
if (.not.allocated(NextRRDate)) allocate (NextRRDate(nrr))
if (.not.allocated(NextRRDatetmp)) allocate (NextRRDatetmp(nrr))
 do irr = 1,nrr
    RR(:,irr) = RRConst(irr) ! initialise all RR to constant values
	RRtmp(:,irr) = RRConst(irr) ! initialise all RR to constant values
	NextRR(:,irr) = RR(:,irr)
	NextRRDate(irr)%day = 99; NextRRDate(irr)%month = 99; NextRRDate(irr)%year = 9999
	NextRRDatetmp(irr)%day = 99; NextRRDatetmp(irr)%month = 99; NextRRDatetmp(irr)%year = 9999
   if (RRFlag(irr) .gt. 0) then
    GenFileName = trim(DataInRootDir) //  trim(DomName) // '/'  &
                  // 'met/' // trim(RRFileName(irr))
    INQUIRE (file=GenFileName, exist=GenFileExistFlag)
		if (GenFileExistFlag) then

		  if (RRFlag(irr).eq.2 ) then ! climatology	
			  iunit = 900 + irr	  
			  CLOSE (iunit)  
			  OPEN (iunit, file=GenFileName, status='old', form='binary', action='read')
			  REWIND (iunit)
			  read(iunit) RRDate, RR(:,irr)
			  do while(.not.EqMonthofYear(RRDate,StartDate)) 
				read(iunit) RRDate, RR(:,irr)
			  enddo
		  elseif (RRFlag(irr).eq.1) then ! time-series
			  iunit = 9000 + irr
			  CLOSE (iunit)  	  
			  OPEN (iunit, file=GenFileName, status='old', form='binary', action='read')
			  REWIND (iunit)
			  read(iunit) RRDate, RR(:,irr)
			  do while(.not.EqMonths(RRDate,StartDate)) 
				read(iunit,IOSTAT=istat) RRDate, RR(:,irr)
				if (istat.eq.-1) then
					write(*,*) 'InitModel: StartDate before FirstRRDate' 
					REWIND (iunit)
					read(iunit) NextRRDatetmp(irr), RRtmp(:,irr)
					write(*,"('Reverting to Climatology for ',a,' prior to ',3i5)") &
					          RRName(irr,1) , NextRRDatetmp(irr)%day,NextRRDatetmp(irr)%month, &
							  NextRRDatetmp(irr)%year
					RRFlag(irr)=2
					exit ! end of file
				endif
			  enddo
			  ! use climatology if StartDate before FirstRRDate
			  if (istat.eq.-1) then
				  iunit = 900+irr
				  GenFileName = trim(DataInRootDir) //  trim(DomName) // '/'  &
					  // 'met/9999010199991231_' // trim(RRFileID(irr))//'_clim.bin'
				  INQUIRE (file=GenFileName, exist=GenFileExistFlag)
				  if (GenFileExistFlag) then
				    CLOSE (iunit)  
					OPEN (iunit, file=GenFileName, status='old', form='binary', action='read')
					REWIND (iunit)
					read(iunit) RRDate, RR(:,irr)	  
					do while(.not.EqMonthofYear(RRDate,StartDate)) 
						read(iunit) RRDate, RR(:,irr)
					enddo
				  else
					write(*,"('File not found:',/,a)") GenFileName
					STOP 'InitModel: StartDate before FirstRRDate and no climatology' 
				  endif
			  endif ! istat eq -1
		  endif ! RRFlag eq 1
		else
		  write(*,"('File not found:',/,a)") GenFileName
		  STOP 'InitModel: RRFile not found'
		end if  ! GenFileExistFlag
		NextRRDate(irr) = RRDate
		NextRR(:,irr)=RR(:,irr)
    end if  ! RRFlag>0

	if (DiagsFlag == 1) then
	!write(*,"('First RR,(',i6,')',3f10.2)") irr,RR(:,irr)
	end if

end do  ! irr

if (ispin==1) then

! * OPEN AND INITIALISE OBSERVATION FILES
!   * NDVI
if (.not.allocated(NextANDVI)) allocate (NextANDVI(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(NDVIFileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 2
CLOSE (40)                                      ! ensure file closed before opening
OPEN (40, file=GenFileName, form='binary', status='old', action='read')
REWIND (40)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,1000                                     ! read NDVI records before StartDate
  read (40,end=2) NextNDVIDate, NextANDVI
  if (DayDifference(StartDate, NextNDVIDate) <= 0) exit ! found first record at/after StartDate
end do
goto 1
2 NextNDVIDate = dmydate(1,1,2100)              ! EOF or no data: set next NDVI in far future
NextANDVI = ErrVal                              !   and data to null
1 continue
if (DiagsFlag == 1) then
  write(*,"('Init NDVIDate, ANDVI:       ',3i6,3f10.2)") NextNDVIDate,NextANDVI(1:1)
end if


!   * LST and LSTTime
if (.not.allocated(NextALST))     allocate (NextALST(np))
if (.not.allocated(NextALSTTime)) allocate (NextALSTTime(np))
if (.not.allocated(NextALSTAngle)) allocate (NextALSTAngle(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(LSTFileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 4
GenFileName2 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(LSTTimeFileName)
INQUIRE (file=GenFileName2, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 4
GenFileName3 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(LSTAngleFileName)
INQUIRE (file=GenFileName3, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 4
CLOSE (41)                                      ! ensure file closed before opening
CLOSE (42)                                      ! ensure file closed before opening
CLOSE (43)                                      ! ensure file closed before opening
OPEN (41, file=GenFileName,  form='binary', status='old', action='read')    ! LST
OPEN (42, file=GenFileName2, form='binary', status='old', action='read')    ! LSTTime
OPEN (43, file=GenFileName3, form='binary', status='old', action='read')    ! LSTAngle
REWIND (41)
REWIND (42)
REWIND (43)
if (DiagsFlag == 1) then
  write(*,"('Open:',a)") trim(GenFileName)
  write(*,"('Open:',a)") trim(GenFileName2)
  write(*,"('Open:',a)") trim(GenFileName3)
end if
do i=1,100000                                     ! read LST and LSTime and LSTangle records before StartDate
  read (41,end=4) NextLSTDate,     NextALST
  read (42,end=4) NextLSTTimeDate, NextALSTTime
  read (43,end=4) NextLSTAngleDate, NextALSTAngle
  NextALST     = NextALST - 273.16              ! degK to degC
  NextALSTTime = NextALSTTime*24.0 + SolarUTCOffset  ! UTC to solar time
  if (DayDifference(StartDate, NextLSTDate) <= 0) exit  ! found first record at/after StartDate
end do
goto 3

4 NextLSTDate = dmydate(1,1,2100)               ! EOF or no data: set next NDVI in far future
NextALST     = ErrVal                           !   and data to null
NextALSTTime = ErrVal
NextALSTAngle = ErrVal
3 continue
if (DiagsFlag == 1) then
  write(*,"('First LSTDate, ALST:        ',3i6,3f10.2)") NextLSTDate,NextALST(1:1)
  write(*,"('First LSTTimeDate, ALSTTime:',3i6,3f10.2)") NextLSTTimeDate,NextALSTTime(1:1)
  write(*,"('First LSTTimeDate, ALSTAngle:',3i6,3f10.2)") NextLSTAngleDate,NextALSTAngle(1:1)
end if


!   * hourly fluxes and correponding time
if (.not.allocated(NextAphiRneth))      allocate (NextAphiRneth(np))
if (.not.allocated(NextAphiHh))         allocate (NextAphiHh(np))
if (.not.allocated(NextAphiEh))         allocate (NextAphiEh(np))
if (.not.allocated(NextAphiNEEh))       allocate (NextAphiNEEh(np))
if (.not.allocated(NextAphiGh))       allocate (NextAphiGh(np))
if (.not.allocated(NextATsoilh))       allocate (NextATsoilh(np))
if (.not.allocated(NextAphiHhTime))     allocate (NextAphiHhTime(np))



GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiRnethFileName)

INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 40
GenFileName2 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiHhFileName)
INQUIRE (file=GenFileName2, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 40
GenFileName3 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiEhFileName)
INQUIRE (file=GenFileName3, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 40
GenFileName4 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiNEEhFileName)
INQUIRE (file=GenFileName4, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 40
GenFileName5 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiGhFileName)
INQUIRE (file=GenFileName5, exist=GenFileExistFlag)
write(*,*) GenFileName
if (.not.(GenFileExistFlag)) goto 40
GenFileName6 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiHhTimeFileName)
INQUIRE (file=GenFileName6, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 40
GenFileName7 = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(TsoilhFileName)

INQUIRE (file=GenFileName7, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 40
CLOSE (551)                                      ! ensure file closed before opening
CLOSE (552)                                      ! ensure file closed before opening
CLOSE (553)                                      ! ensure file closed before opening
CLOSE (554)                                      ! ensure file closed before opening
CLOSE (555)                                      ! ensure file closed before opening
CLOSE (556)                                      ! ensure file closed before opening
CLOSE (557)                                      ! ensure file closed before opening
OPEN (551, file=GenFileName,  form='binary', status='old', action='read')    ! phiRneth
OPEN (552, file=GenFileName2,  form='binary', status='old', action='read')    ! phiHh
OPEN (553, file=GenFileName3,  form='binary', status='old', action='read')    ! phiEh
OPEN (554, file=GenFileName4,  form='binary', status='old', action='read')    ! phiNEEh
OPEN (555, file=GenFileName5,  form='binary', status='old', action='read')    ! phiGh
OPEN (556, file=GenFileName6, form='binary', status='old', action='read')    ! phiHhTime
OPEN (557, file=GenFileName7, form='binary', status='old', action='read')    ! Tsoilh
REWIND (551)
REWIND (552)
REWIND (553)
REWIND (554)
REWIND (555)
REWIND (556)
REWIND (557)
if (DiagsFlag == 1) then
  write(*,"('Open:',a)") trim(GenFileName)
  write(*,"('Open:',a)") trim(GenFileName2)
  write(*,"('Open:',a)") trim(GenFileName3)
  write(*,"('Open:',a)") trim(GenFileName4)
  write(*,"('Open:',a)") trim(GenFileName5)
  write(*,"('Open:',a)") trim(GenFileName6)
end if
do i=1,100000                                     ! read hourly flux records before StartDate
  read (551,end=40) NextAphiRnethDate,     NextAphiRneth
  read (552,end=40) NextAphiHhDate,        NextAphiHh
  read (553,end=40) NextAphiEhDate,        NextAphiEh
  read (554,end=40) NextAphiNEEhDate,      NextAphiNEEh
  read (555,end=40) NextAphiGhDate,      NextAphiGh
  read (557,end=40) NextATsoilhDate,      NextATsoilh
  read (556,end=40) NextAphiHhTimeDate,    NextAphiHhTime
  if (DayDifference(StartDate, NextAphiHhDate) <= 0) then
	exit  ! found first record at/after StartDate
  endif
end do
goto 30

!40 write(*,*) NextAphiRnethDate, NextAphiHhDate, NextAphiEhDate,NextAphiNEEhDate,NextAphiGhDate,NextATsoilhDate, NextAphiHhTimeDate
40 NextAphiHhDate = dmydate(1,1,2100)               ! EOF or no data: set next NDVI in far future
NextAphiHh     = ErrVal                           !   and data to null
NextAphiHhTime = ErrVal
30 continue
if (DiagsFlag == 1) then
  write(*,"('First AphiHhDate, AphiHh:        ',3i6,3f10.2)") NextAphiHhDate,NextAphiHh(1:1)
  write(*,"('First AphiHhTimeDate, AphiHhTime:',3i6,3f10.2)") NextAphiHhTimeDate,NextAphiHhTime(1:1)
end if



! * OPEN and Initialsise DAILy Flux files
!* phiRnet
if (.not.allocated(NextAphiRnet)) allocate (NextAphiRnet(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiRnetFileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 201
CLOSE (661)                                      ! ensure file closed before opening
OPEN (661, file=GenFileName, form='binary', status='old', action='read')
REWIND (661)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiRnet records before StartDate
  read (661,end=201) NextAphiRnetDate, NextAphiRnet
  if (DayDifference(StartDate, NextAphiRnetDate) <= 0) exit ! found first record at/after StartDate
end do
goto 101
201 NextAphiRnetDate = dmydate(1,1,2100)              ! EOF or no data: set next phiRnet in far future
NextAphiRnet = ErrVal                              !   and data to null
101 continue
if (DiagsFlag == 1) then
  write(*,"('Init phiRnetDate, AphiRnet:       ',3i6,3f10.2)") NextAphiRnetDate,NextAphiRnet(1:1)
end if

!* phiH
if (.not.allocated(NextAphiH)) allocate (NextAphiH(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiHFileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 202
CLOSE (662)                                      ! ensure file closed before opening
OPEN (662, file=GenFileName, form='binary', status='old', action='read')
REWIND (662)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiH records before StartDate
  read (662,end=202) NextAphiHDate, NextAphiH
  if (DayDifference(StartDate, NextAphiHDate) <= 0) exit ! found first record at/after StartDate
end do
goto 102
202 NextAphiHDate = dmydate(1,1,2100)              ! EOF or no data: set next phiH in far future
NextAphiH = ErrVal                              !   and data to null
102 continue
if (DiagsFlag == 1) then
  write(*,"('Init phiHDate, AphiH:       ',3i6,3f10.2)") NextAphiHDate,NextAphiH(1:1)
end if

!* phiE
if (.not.allocated(NextAphiE)) allocate (NextAphiE(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiEFileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 203
CLOSE (663)                                      ! ensure file closed before opening
OPEN (663, file=GenFileName, form='binary', status='old', action='read')
REWIND (663)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiE records before StartDate
  read (663,end=203) NextAphiEDate, NextAphiE
  if (DayDifference(StartDate, NextAphiEDate) <= 0) exit ! found first record at/after StartDate
end do
goto 103
203 NextAphiEDate = dmydate(1,1,2100)              ! EOF or no data: set next phiE in far future
NextAphiE = ErrVal                              !   and data to null
103 continue
if (DiagsFlag == 1) then
  write(*,"('Init phiEDate, AphiE:       ',3i6,3f10.2)") NextAphiEDate,NextAphiE(1:1)
end if

!* phiNEE
if (.not.allocated(NextAphiNEE)) allocate (NextAphiNEE(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(phiNEEFileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 204
CLOSE (664)                                      ! ensure file closed before opening
OPEN (664, file=GenFileName, form='binary', status='old', action='read')
REWIND (664)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read NDVI records before StartDate
  read (664,end=204) NextAphiNEEDate, NextAphiNEE
  if (DayDifference(StartDate, NextAphiNEEDate) <= 0) exit ! found first record at/after StartDate
end do
goto 104
204 NextAphiNEEDate = dmydate(1,1,2100)              ! EOF or no data: set next NDVI in far future
NextAphiNEE = ErrVal                              !   and data to null
104 continue
if (DiagsFlag == 1) then
  write(*,"('Init phiNEEDate, AphiNEE:       ',3i6,3f10.2)") NextAphiNEEDate,NextAphiNEE(1:1)
end if

! * OPEN and Initialsise DAILY soil moisture files
!* sm0008
if (.not.allocated(NextAsm0008)) allocate (NextAsm0008(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(sm0008FileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 701
CLOSE (771)                                      ! ensure file closed before opening
OPEN (771, file=GenFileName, form='binary', status='old', action='read')
REWIND (771)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiRnet records before StartDate
  read (771,end=701) NextAsm0008Date, NextAsm0008
  if (DayDifference(StartDate, NextAsm0008Date) <= 0) exit ! found first record at/after StartDate
end do
goto 711
701 NextAsm0008Date = dmydate(1,1,2100)              ! EOF or no data: set next Asm0008 in far future
NextAsm0008 = ErrVal                              !   and data to null
711 continue
if (DiagsFlag == 1) then
  write(*,"('Init Asm0008Date, Asm0008:       ',3i6,3f10.2)") NextAsm0008Date,NextAsm0008(1:1)
end if

!* sm0090
if (.not.allocated(NextAsm0090)) allocate (NextAsm0090(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(sm0090FileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 702
CLOSE (772)                                      ! ensure file closed before opening
OPEN (772, file=GenFileName, form='binary', status='old', action='read')
REWIND (772)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiRnet records before StartDate
  read (772,end=702) NextAsm0090Date, NextAsm0090
  if (DayDifference(StartDate, NextAsm0090Date) <= 0) exit ! found first record at/after StartDate
end do
goto 712
702 NextAsm0090Date = dmydate(1,1,2100)              ! EOF or no data: set next Asm0090 in far future
NextAsm0090 = ErrVal                              !   and data to null
712 continue
if (DiagsFlag == 1) then
  write(*,"('Init Asm0090Date, Asm0090:       ',3i6,3f10.2)") NextAsm0090Date,NextAsm0090(1:1)
end if

!* sm0030
if (.not.allocated(NextAsm0030)) allocate (NextAsm0030(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(sm0030FileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 703
CLOSE (773)                                      ! ensure file closed before opening
OPEN (773, file=GenFileName, form='binary', status='old', action='read')
REWIND (773)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiRnet records before StartDate
  read (773,end=703) NextAsm0030Date, NextAsm0030
  if (DayDifference(StartDate, NextAsm0030Date) <= 0) exit ! found first record at/after StartDate
end do
goto 713
703 NextAsm0030Date = dmydate(1,1,2100)              ! EOF or no data: set next Asm0030 in far future
NextAsm0030 = ErrVal                              !   and data to null
713 continue
if (DiagsFlag == 1) then
  write(*,"('Init Asm0030Date, Asm0030:       ',3i6,3f10.2)") NextAsm0030Date,NextAsm0030(1:1)
end if

!* sm3060
if (.not.allocated(NextAsm3060)) allocate (NextAsm3060(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(sm3060FileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 704
CLOSE (774)                                      ! ensure file closed before opening
OPEN (774, file=GenFileName, form='binary', status='old', action='read')
REWIND (774)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiRnet records before StartDate
  read (774,end=704) NextAsm3060Date, NextAsm3060
  if (DayDifference(StartDate, NextAsm3060Date) <= 0) exit ! found first record at/after StartDate
end do
goto 714
704 NextAsm3060Date = dmydate(1,1,2100)              ! EOF or no data: set next Asm3060 in far future
NextAsm3060 = ErrVal                              !   and data to null
714 continue
if (DiagsFlag == 1) then
  write(*,"('Init Asm3060Date, Asm3060:       ',3i6,3f10.2)") NextAsm3060Date,NextAsm3060(1:1)
end if

!* sm6090
if (.not.allocated(NextAsm6090)) allocate (NextAsm6090(np))
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // trim(sm6090FileName)
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (.not.(GenFileExistFlag)) goto 705
CLOSE (775)                                      ! ensure file closed before opening
OPEN (775, file=GenFileName, form='binary', status='old', action='read')
REWIND (775)
if (DiagsFlag == 1) then
  write(*,"('Initialise: ',a)") trim(GenFileName)
end if
do i=1,5000                                     ! read phiRnet records before StartDate
  read (775,end=705) NextAsm6090Date, NextAsm6090
  if (DayDifference(StartDate, NextAsm6090Date) <= 0) exit ! found first record at/after StartDate
end do
goto 715
705 NextAsm6090Date = dmydate(1,1,2100)              ! EOF or no data: set next Asm6090 in far future
NextAsm6090 = ErrVal                              !   and data to null
715 continue
if (DiagsFlag == 1) then
  write(*,"('Init Asm6090Date, Asm6090:       ',3i6,3f10.2)") NextAsm6090Date,NextAsm6090(1:1)
end if



!********************************************************************************

! * OPEN AND INITIALISE CATCHMENT DISCHARGE DATA FILES 
DisCMDataExist(1:nc) = .true.           ! assume all monthly data files exist
DisCDDataExist(1:nc) = .false.           ! assume all daily data files exist
!   * Step 1: use FindUniqueRegions to detect catchments (done in ReadControlFile)
!   * Step 2: open and position MONTHLY discharge data files, catchment by catchment
do ic=1,nc
  if (DisCMDataExist(ic)) then
    write(CatchName,"(i6)") iCatchID(ic)
    GenFileName = trim(DataInRootDir) // 'DischargeData.201111/'  &
                  // trim(CatchName) // 'monrun.dat'                ! sep-2005 data (monthly)
!    GenFileName = trim(DataInRootDir) // 'DischargeData.200610/'  &
!                  // trim(CatchName) // 'monrun_new.dat'            ! oct-2006 data (monthly)
    INQUIRE (file=GenFileName, exist=GenFileExistFlag)
    if (.not.(GenFileExistFlag)) goto 6
    OPEN (unit=2205+ic, file=GenFileName, status='old', action='read')
    REWIND (unit=2205+ic)
    read (2205+ic,"()",end=6)            ! skip header
    DischargeDate = dmydate(0,0,0)      ! initialise DischargeDate
    do while (DischargeDate .lt. StartDate)
                                        ! skip records prior to StartDate
      read (2205+ic,205,end=6) ObsYear, ObsMonth
      205 format (i6,i2,f12.2)                                      ! sep-2005 data (monthly)
!      read (205+ic,*,end=6) ObsYear, ObsMonth                       ! oct-2006 data (monthly)
      DischargeDate = dmydate(28, ObsMonth, ObsYear)
      DischargeDate = dmydate(DaysInMonth(DischargeDate), ObsMonth, ObsYear)
    end do
    backspace(2205+ic)                   ! ready to read first record after StartDate
    goto 5
    6 DisCMDataExist(ic) = .false.      ! no data (no file, or EOF before StartDate)
    5 continue
  end if
end do
!   * Step 3: open and position DAILY discharge data files, catchment by catchment
do ic=1,nc
  if (DisCDDataExist(ic)) then
    write(CatchName,"(i6)") iCatchID(ic)
    GenFileName = trim(DataInRootDir) // 'DischargeData.200509/'  &
                  // trim(CatchName) // 'dayrun.dat'                ! sep-2005 data (daily)
!    GenFileName = trim(DataInRootDir) // 'DischargeData.200610/'  &
!                  // trim(CatchName) // 'dayrun.dat'                ! oct-2006 data (daily)
    INQUIRE (file=GenFileName, exist=GenFileExistFlag)
    if (.not.(GenFileExistFlag)) goto 8
    OPEN (unit=1205+ic, file=GenFileName, status='old', action='read')
    REWIND (unit=1205+ic)
    read (1205+ic,"()",end=8)           ! skip header
    DischargeDate = dmydate(0,0,0)      ! initialise DischargeDate
    do while (DischargeDate .lt. StartDate)
                                        ! skip records prior to StartDate
      read (1205+ic,1205,end=8) ObsYear, ObsMonth, ObsDay
      1205 format (i6,2i2,f10.3)                                    ! sep-2005 data (daily)
 !    read (1205+ic,*,end=8) Obsday, ObsMonth, ObsYear              ! oct-2006 data (daily)
 !     DischargeDate = dmydate(ObsDay, ObsMonth, ObsYear)
    end do
    backspace(1205+ic)                  ! ready to read first record after StartDate
    goto 7
    8 DisCDDataExist(ic) = .false.      ! no data (no file, or EOF before StartDate)
    7 continue
  end if
end do
!   * Step 4: read catchment properties
if (DiagsFlag == 1) then
  write(*,"(/,'Catchment properties (ic, iCatchID, iCatchIDDum, TZRunCat, TZLchCat)')")
end if
GenFileName = trim(DataInRootDir) // 'DischargeData.200509/'      &
              // '00 CatchmentInformationEdited+Params.csv'
OPEN (unit=1200, file = GenFileName, status='old', action='read')
do ic=1,nc
  REWIND (unit=1200)
  do i=1,7                              ! skip headers
    read (1200,"()")
  end do
  do
    read (1200,*,end=9) dum(1), dumchar, CatchIDdum,            &
                        dum(1:2), dumchar, dum(1:14), TZRunCat(ic), TZLchCat(ic)
    if ( nint(CatchIDDum) == iCatchID(ic) ) goto 10
  end do
  9  write(*,"('iCatchID =',i8)") iCatchID(ic)
 ! STOP "InitModel: Catchment properties not found"         ! vh 24/12/09: commented out "STOP"
   write(*,*) "InitModel: Catchment properties not found"
  10 continue
  if (DiagsFlag == 1) then
    write(*,*) ic, iCatchID(ic), CatchIDDum, TZRunCat(ic), TZLchCat(ic)
  end if
end do
CLOSE (unit=1200)

! * INITIALISE ZRunCD, ZLchCD (low-pass filtered FWRunC, FWLchC in obs model)
ZRunCD = 0.0
ZLchCD = 0.0



	! * INITIALISE ACCUMULATORS
	!   * Initialise time averages of spatial averages of np-arrays (iTSOPFlag=0)
	do AvType = 1,3     ! AvType = (1,2,3) gives (month,ann,run averages)
						! arrays: XXSpStat(1,nxx,5), XXSpTAv(1,nxx,nAvType)
	  CALL TimeAverage (0, MMSpStat(:,:,1),  MMSpTAv(:,:,AvType),  MMSpTAvCount(:,:,AvType))
	  CALL TimeAverage (0, RRSpStat(:,:,1),  RRSpTAv(:,:,AvType),  RRSpTAvCount(:,:,AvType))
	  CALL TimeAverage (0, XXSpStat(:,:,1),  XXSpTAv(:,:,AvType),  XXSpTAvCount(:,:,AvType))
	  CALL TimeAverage (0, FFSpStat(:,:,1),  FFSpTAv(:,:,AvType),  FFSpTAvCount(:,:,AvType))
	  CALL TimeAverage (0, DDSpStat(:,:,1),  DDSpTAv(:,:,AvType),  DDSpTAvCount(:,:,AvType))
	  CALL TimeAverage (0, ZZPSpStat(:,:,1), ZZPSpTAv(:,:,AvType), ZZPSpTAvCount(:,:,AvType))
	  CALL TimeAverage (0, AAPSpStat(:,:,1), AAPSpTAv(:,:,AvType), AAPSpTAvCount(:,:,AvType))
	end do

	!   * Initialise nc-array time averages for each catchment (iTSOPFlag=0) 
	!     * NOTE: for catchment ZZ, AA, do whole-run time avs as well as month, year,
	!       for output to run-average discharge file
	do AvType = 1,3     ! AvType = (1,2,3) gives (month,ann,run averages)
						! arrays: ZZC1(nc,nzzC), ZZCTAv(nc,nzzC,nAvType)
	  CALL TimeAverage (0, ZZC1, ZZCTAv(:,:,AvType), ZZCTAvCount(:,:,AvType))
	  CALL TimeAverage (0, AAC1, AACTAv(:,:,AvType), AACTAvCount(:,:,AvType))
	end do

	!   * Initialise np-array time averages for maps (iTSOPFlag=0) 
	do AvType = 1,3     ! AvType = (1,2,3) gives (month,ann,run averages)
						! arrays: XXTAv(np,nxx,nAvType), XXTAvCount(np,nxx,nAvType)
	  CALL TimeAverage (0, MM,   MMTAv(:,:,AvType),  MMTAvCount(:,:,AvType))
	  CALL TimeAverage (0, RR,   RRTAv(:,:,AvType),  RRTAvCount(:,:,AvType))
	  CALL TimeAverage (0, XX1,  XXTAv(:,:,AvType),  XXTAvCount(:,:,AvType))
	  CALL TimeAverage (0, FF,   FFTAv(:,:,AvType),  FFTAvCount(:,:,AvType))
	  CALL TimeAverage (0, DD,   DDTAv(:,:,AvType),  DDTAvCount(:,:,AvType))
	  CALL TimeAverage (0, ZZP1, ZZPTAv(:,:,AvType), ZZPTAvCount(:,:,AvType))
	  CALL TimeAverage (0, AAP1, AAPTAv(:,:,AvType), AAPTAvCount(:,:,AvType))
	  CALL TimeAverageh (0, ZZPh1, ZZPhTAv(:,:,:,AvType), ZZPhTAvCount(:,:,:,AvType))
	  CALL TimeAverageh (0, AAPh1, AAPhTAv(:,:,:,AvType), AAPhTAvCount(:,:,:,AvType))
	end do

	! * OPEN AND INITIALISE OUTPUT FILES
	!   * Open time series output files and write headers
	!     * arrays: XX(np,nxx), XXSpStat(1,nxx,5) (5 for mean, sd, max, min, count)
	!     * TS output files are opened only when needed according to DomTSOutputFlag
	CALL SpatialStats (VV, VVSpStat)        ! find spatial averages of VV
	CALL SpatialStats (XX0, XXinitSpStat)   ! and of XX0 (= XXinit)
	do AvType=3,0,-1                        ! daily, monthly, annual, whole-run TS
	  if (4-AvType > DomTSOutputFlag) exit  ! DomTSOutputFlag=0 for sensitivity tests, so no TS output
	  if (ipDiag == 0) then                 ! TS for domain spatial averages
		CALL TimeSeriesOutput (0, AvType, ipDiag,         & 
		   UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1),      &
		   TTime, MM(1,:),RR(1,:), XX1(1,:), FF(1,:), DD(1,:), ZZP1(1,:), AAP1(1,:))
											! last line of args are dummies for initialisation call
	  else                                  ! TS for point ip
		CALL TimeSeriesOutput (0, AvType, ipDiag,         &
		   UU, VV(ipDiag,:), XX0(ipDiag,:),               &
		   TTime, MM(1,:),RR(1,:), XX1(1,:), FF(1,:), DD(1,:), ZZP1(1,:), AAP1(1,:))
											! last line of args are dummies for initialisation call
	  end if
	end do

	!   * Open catchment time series output files and write headers for day, mth, ann TS
	!     (Note that whole-run TS is initialised in subroutine CatchmentRunAvgOutput)
	do AvType=2,0,-1                        ! daily, monthly, annual TS
	  if (3-AvType>CatchTSOutputFlag) exit  ! CatchTSOutputFlag=0 for sensitivity tests, so no TS output
	  CALL CatchmentTSOutput                                                      &
			 (0, AvType, UU, VVSpStat(1,:,1), XXinitSpStat(1,:,1), iCatchID,      &
			  1, 2, TTime, ZZC1, AAC1)
	end do

	! OPEN DD TS binary files
	if (MapTSOutputFlag.eq.4) then
	 AvType = 0
	 CALL TimeSeriesOutput_bin (0,DDTSOPFlag, AvType,'day', TTime,DD,DDName(:,1),1)
	 CALL TimeSeriesOutput_bin (0,XXTSOPFlag, AvType,'day', TTime,XX1,XXName(:,1),2)
	 CALL TimeSeriesOutput_bin (0,FFTSOPFlag, AvType,'day', TTime,FF,FFName(:,1),3)
	endif
	if (MapTSOutputFlag >= 3) then 
		AvType = 1
		CALL TimeSeriesOutput_bin (0,DDTSOPFlag, AvType,'mth', TTime,DD,DDName(:,1),1)
		CALL TimeSeriesOutput_bin (0,XXTSOPFlag, AvType,'mth', TTime,XX1,XXName(:,1),2)
		CALL TimeSeriesOutput_bin (0,FFTSOPFlag, AvType,'mth', TTime,FF,FFName(:,1),3)
	endif
	if (MapTSOutputFlag >= 2) then 
		AvType = 2
		CALL TimeSeriesOutput_bin (0,DDTSOPFlag, AvType,'ann', TTime,DD,DDName(:,1),1)
		CALL TimeSeriesOutput_bin (0,XXTSOPFlag, AvType,'ann', TTime,XX1,XXName(:,1),2)
		CALL TimeSeriesOutput_bin (0,FFTSOPFlag, AvType,'ann', TTime,FF,FFName(:,1),3)
	endif
	if (MapTSOutputFlag >= 1) then 
		AvType = 3
		CALL TimeSeriesOutput_bin (0,DDTSOPFlag, AvType,'run', TTime,DD,DDName(:,1),1)
		CALL TimeSeriesOutput_bin (0,XXTSOPFlag, AvType,'run', TTime,XX1,XXName(:,1),2)
		CALL TimeSeriesOutput_bin (0,FFTSOPFlag, AvType,'run', TTime,FF,FFName(:,1),3)
	endif
		

	!   * Initialise PEST output files
	if (PESTOutputFlag >= 2) then
	  CALL WritePESTOutput (Time0,TTDate,AAP1,AAC1,ZZP1,ZZC1,FFTAv,VV,0,AAPH1,ZZPH1)  
	end if

end if ! (if ispin==1))

!if (DiagsFlag == 1) then
!  write(*,"(/,'InitModel complete, start run (y/n) ? ',$)") 
!  read (*,*) YN
!  if (YN /= 'y') STOP "Exit InitModel"
!end if

!open(unit=7,file='test.out',status='unknown')

END SUBROUTINE InitModel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE InitCABLE(UU)
Implicit none
real(sp),intent(in)::  UU(:)      ! UU(nuu) = spatially const params
INTEGER:: k, j,i, is
real:: totdepth                                    
REAL(sp):: totroot(np) 
 LOGICAL    :: vegparmnew   ! using new format input file 
 LOGICAL    :: spinup ! should the model spinup to soil state equilibrium?
 LOGICAL    :: spinConv ! has spinup converged? 
INTEGER(i_d) :: iday, imonth, iyear
CHARACTER(LEN=13) :: timestring                  ! will have a fixed value "seconds since"   



	CALL PointAll (TargetUU=UU, TargetXX=XX0,TargetRR=RR)
	mp = np
	


	ALLOCATE(latitude(np),longitude(np))  
	latitude = LatDeg ! set latitude
    longitude = LongDeg ! set longitude
	smoy=1
    spinConv = .FALSE.
	vegparmnew = .TRUE.  ! using new format
    spinup = .FALSE.  ! do we spin up the model?
	logn = 888      ! log file number - declared in input module
	model_structure_flag%canopy = 'canopy_vh' ! 'default','spatial' 'canopy_vh'
    model_structure_flag%soil = 'sli' ! 'default','soilsnow','sli','spatial'
    model_structure_flag%photosynthesis = 'default' ! 'default','hawkesbury','spatial'	
    model_structure_flag%sli_litter = 'resistance' ! 'default','on','off','resistance'
    model_structure_flag%sli_isotope = 'off' ! 'default', 'off','HDO','H218O','spatial'
    model_structure_flag%sli_coupled = 'coupled' ! 'coupled','uncoupled' energy and water
	check%ranges     = .TRUE.  ! variable ranges, input and output
    check%energy_bal = .TRUE.  ! energy balance
    check%mass_bal   = .TRUE.  ! water/mass balance

	! filenames for default soil and veg params
	filename%veg = trim(DataInRootDir) //'surface_data/def_veg_params_igbp.txt'
	filename%soil = trim(DataInRootDir) //'surface_data/def_soil_params.txt'
	filename%type =trim(DataInRootDir) //'surface_data/gridinfo1992_igbpGab.txt'
	filename%log = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // 'CABLElog.txt'
	filename%inits   = trim(DataInRootDir) //'surface_data/Mk3LsurfClimatology.nc'
	filename%LAI     = trim(DataInRootDir) //'surface_data/LAI_Monthly_Global.nc'

   IF(vegparmnew) THEN
          nvegt = 17
    ELSE
          nvegt = 13
   ENDIF

    nsoilt = 9


	OPEN(logn,FILE=filename%log)
	
	if (CableParamFlag.eq.2) then  ! Use global defaults
		CALL load_parameters(met,air,ssoil,veg,bgc,soil,canopy, &
       & rough,rad,sum_flux,bal,model_structure,logn,vegparmnew)

	else


		mp_patch = np*ntile   ! grassy and woody tile for each pixel (1=grassy, 2= woody)
		max_vegpatches = ntile
		if (CableParamFlag.eq.1) then
			tile_area(:,1) = 1.-fWoody
			tile_area(:,2) = fWoody
		else
			tile_area(:,1)=1.0
		endif
		do j=1,ntile
			tile_index(:,j) = (/(i,i=j,mp_patch,ntile)/)
		enddo


		CALL set_model_structure(mp_patch,model_structure)
		! Allocate CABLE's main variables:
		CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssoil, &
            sum_flux,veg,model_structure,mp_patch)
	
		
		! set parameters the same for each tile
		do k=1,ntile
			soil%albsoil(k:mp_patch:ntile) = albsoil
			soil%bch(k:mp_patch:ntile) = B1mult*min(bch1,16.0)
			soil%silt(k:mp_patch:ntile) = silt1
			soil%clay(k:mp_patch:ntile) = clay1
			soil%sand = 1.0 - soil%silt - soil%clay
			soil%css(k:mp_patch:ntile) = css1
			soil%hyds(k:mp_patch:ntile) = Ksat1Mult*max(hyds1,1.0e-8) ! m s-1
			soil%rs20(k:mp_patch:ntile) = rs20
			soil%sfc(k:mp_patch:ntile) = sfc1
			soil%ssat(k:mp_patch:ntile) = max(0.4,WVolSat1)
			soil%sucs(k:mp_patch:ntile) = max(sucs1,-2.0) * Psie1Mult         
			soil%swilt(k:mp_patch:ntile) = min(0.2,swilt1) 
			soil%rhosoil(k:mp_patch:ntile) = rhosoil1
			soil%clitt(k:mp_patch:ntile) = clittequil 
			soil%zeta(k:mp_patch:ntile) = zeta
			soil%fsatmax(k:mp_patch:ntile) = fsatmax

			! B horizon parameters (for use in soil-litter)
			soil%bchB(k:mp_patch:ntile) = B2mult*min(bch2,16.0)
			soil%siltB(k:mp_patch:ntile) = silt2
			soil%clayB(k:mp_patch:ntile) = clay2
			soil%cssB(k:mp_patch:ntile) = css2
			soil%hydsB(k:mp_patch:ntile) = Ksat2Mult*max(hyds2,1.0e-8)  ! m s-1
			soil%ssatB(k:mp_patch:ntile) = max(0.4,WVolSat2)
			soil%sucsB(k:mp_patch:ntile) = max(sucs2,-2.0) * Psie2Mult           
			soil%swiltB(k:mp_patch:ntile) = min(0.2,swilt2)
			soil%rhosoilB(k:mp_patch:ntile) = rhosoil2
			soil%sfcB(k:mp_patch:ntile) = sfc2
			soil%depthA(k:mp_patch:ntile) = Zsoil1
			!soil%depthA(k:mp_patch:ntile) = 30.0 ! temp test
			soil%depthB(k:mp_patch:ntile)= Zsoil2

			! veg parameters
			veg%vegcf(k:mp_patch:ntile) = vegcf
			veg%iveg(k:mp_patch:ntile) = int(iveg)
			veg%canst1(k:mp_patch:ntile) = canst1
			veg%dleaf(k:mp_patch:ntile) = dleaf    
			veg%hc(k:mp_patch:ntile) = hc
			veg%xfang(k:mp_patch:ntile) = xfang
			veg%rp20(k:mp_patch:ntile) = rp20
			veg%rpcoef(k:mp_patch:ntile) = rpcoef
			veg%shelrb(k:mp_patch:ntile) = shelrb
			veg%frac4(k:mp_patch:ntile) = frac4
			if (k.eq.2) veg%frac4(k:mp_patch:ntile)=0.0
			veg%tminvj(k:mp_patch:ntile) = tminvj
			veg%tmaxvj(k:mp_patch:ntile) = tmaxvj
			veg%vbeta(k:mp_patch:ntile) = vbeta	
			veg%meth(k:mp_patch:ntile) = 1 ! canopy turbulence parameterisation method: 0 or 1
			veg%ejmax(k:mp_patch:ntile) = ratioJV ! update later by multiplying by vcmax
			veg%a1c3(k:mp_patch:ntile) = a1
			veg%d0c3(k:mp_patch:ntile) = ds0
			veg%rootbeta(k:mp_patch:ntile) = rootbeta
			veg%refl(k:mp_patch:ntile,1) = min(max(scattVIS/2., 0.05),0.3)
			veg%refl(k:mp_patch:ntile,2) = min(max(scattNIR/2., 0.2),0.9)
			veg%refl(k:mp_patch:ntile,3) = 0.02

			bgc%cplant(k:mp_patch:ntile,1) = cplant1
			bgc%cplant(k:mp_patch:ntile,2) = cplant2
			bgc%cplant(k:mp_patch:ntile,3) = cplant3
			bgc%csoil(k:mp_patch:ntile,1) = csoil1
			bgc%csoil(k:mp_patch:ntile,2) = csoil2

			bgc%ratecp(1) =ratecp1
			bgc%ratecp(2) =ratecp2
			bgc%ratecp(3) =ratecp3
			bgc%ratecs(1) =ratecp1
			bgc%ratecs(2) =ratecp2
			
			rough%za_uv(k:mp_patch:ntile) = za ! lowest atm. model layer/reference height
			rough%za_tq(k:mp_patch:ntile) = za ! lowest atm. model layer/reference height
			rad%latitude(k:mp_patch:ntile) = LatDeg
			rad%longitude(k:mp_patch:ntile) = LongDeg	
			canopy%cansto(k:mp_patch:ntile) = cansto   ! canopy water storage (mm or kg/m2)		

		enddo


		
		if (CableParamFlag.eq.1) then  ! separate veg params for woody and grassy veg tiles
			veg%dleaf(1:mp_patch:ntile) = dleaf_g
			veg%vcmax(1:mp_patch:ntile) = vcmax_g
			veg%hc(1:mp_patch:ntile) = hc_g
			veg%xfang(1:mp_patch:ntile) = xfang_g
			veg%rp20(1:mp_patch:ntile) = rp20_g
			veg%vbeta(1:mp_patch:ntile) = vbeta_g
			veg%rootbeta(1:mp_patch:ntile) = rootbeta_g
			veg%F10(1:mp_patch:ntile) = F10_g
			veg%ZR(1:mp_patch:ntile) = ZR_g
			veg%gamma(1:mp_patch:ntile) = 10.**(loggamma_g)

			veg%dleaf(2:mp_patch:ntile) = dleaf_w
			veg%vcmax(2:mp_patch:ntile) = vcmax_w
			veg%hc(2:mp_patch:ntile) = hc_w
			veg%xfang(2:mp_patch:ntile) = xfang_w
			veg%rp20(2:mp_patch:ntile) = rp20_w
			veg%vbeta(2:mp_patch:ntile) = vbeta_w
			veg%rootbeta(2:mp_patch:ntile) = rootbeta_w
			veg%F10(2:mp_patch:ntile) = F10_w
			veg%ZR(2:mp_patch:ntile) = ZR_w
			veg%gamma(2:mp_patch:ntile) = 10.**(loggamma_w)

		elseif (CableParamFlag.eq.0) then ! separate veg params for woody and grassy veg (single tile)
		 WHERE (veg%iveg==10)                  ! grassy veg
	        veg%dleaf = dleaf_g
			veg%vcmax = vcmax_g
			veg%hc = hc_g
			veg%xfang = xfang_g
			veg%rp20 = rp20_g
			veg%vbeta = vbeta_g
			veg%rootbeta = rootbeta_g
			veg%F10 = F10_g
			veg%ZR = ZR_g
			veg%gamma = 10.**(loggamma_g)
		  ELSEWHERE                            ! woody veg
			veg%dleaf = dleaf_w
			veg%vcmax = vcmax_w
			veg%hc = hc_w
			veg%xfang = xfang_w
			veg%rp20 = rp20_w
			veg%vbeta = vbeta_w
			veg%rootbeta = rootbeta_w
			veg%F10 = F10_w
			veg%ZR = ZR_w
			veg%gamma = 10.**(loggamma_w)
		  ENDWHERE

		endif

        soil%zse = (/.022, 0.058, 0.07, .15, 0.30,&
	             0.30, 0.30, 1.20, 3.0, 4.5/)   ! soil layer thicknesses
	

		! derived parameters
	!	veg%vcmax = vcmax ! special for LOCAT
		veg%ejmax = veg%ejmax*veg%vcmax
		veg%taul = veg%refl
		soil%swilt_vec = SPREAD(soil%swilt,2,ms)
		soil%ssat_vec = SPREAD(soil%ssat,2,ms)
		soil%sfc_vec = SPREAD(soil%sfc,2,ms)
		soil%cnsd  = soil%sand*0.3 + soil%clay*0.25 &
            + soil%silt*0.265 ! set dry soil thermal conductivity [W/m/K]
		soil%hsbh  = soil%hyds*ABS(soil%sucs)*soil%bch        !difsat*etasat
		soil%ibp2  = NINT(soil%bch)+2
		soil%i2bp3 = 2*NINT(soil%bch)+3
		soil%sandB = 1.0 - soil%siltB - soil%clayB
		soil%cnsdB  = soil%sandB*0.3 + soil%clayB*0.25 &
            + soil%siltB*0.265 ! set dry soil thermal conductivity [W/m/K]
		soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints:
		soil%zshh(ms+1)=0.5*soil%zse(ms)
		soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
		soil%isoilm = 2  !(set to 9 for ice, otherwise any integer is ok)
	

! root density distribution


	! calculate vegin%froot from using rootbeta and soil depth 
    ! (Jackson et al. 1996, Oceologica, 108:389-411)
    totroot(:) = 1.0-veg%rootbeta(:)**(sum(soil%zse)*100.0) 
    totdepth = 0.0    
    do is=1,ms 
       totdepth = totdepth + soil%zse(is)*100.0
       veg%froot(:,is) = min(1.0,1.0-veg%rootbeta(:)**totdepth) 
    enddo
    do is=ms,2,-1
       veg%froot(:,is) = veg%froot(:,is)-veg%froot(:,is-1)
    enddo


 ! Construct derived parameters and zero initialisations, regardless 
    ! of where parameters and other initialisations have loaded from:
    CALL derived_parameters(soil,sum_flux,bal,ssoil,model_structure)

  
    ssoil%ssdnn  = 140.0 ! overall snow density (kg/m3)
    ssoil%snowd  = 0.0   ! snow liquid water equivalent depth (mm or kg/m2)
    ssoil%osnowd = 0.0   ! snow depth prev timestep (mm or kg/m2)
    ssoil%snage  = 0.0   ! snow age
    ssoil%isflag = 0     ! snow layer scheme flag (0 = no or little snow, 1=snow)
    ! snmin  = 0.11  ! switch b/w 1 and 3 layer snow at this depth (m)
    ssoil%wbice  = 0.0   ! soil ice 
    ssoil%tggsn  = 273.1 ! snow temperature per layer (K)
    ssoil%ssdn   = 140.0 ! snow density per layer (kg/m3)
    ssoil%smass  = 0.0   ! snow mass per layer (kg/m^2)
    ssoil%runoff = 0.0   ! runoff total = subsurface + surface runoff
    ssoil%rnof1  = 0.0   ! surface runoff (mm/timestepsize)
    ssoil%rnof2  = 0.0   ! deep drainage (mm/timestepsize)
    ssoil%rtsoil = 100.0 ! turbulent resistance for soil
    canopy%ga     = 0.0   ! ground heat flux (W/m2)
    canopy%dgdtg  = 0.0   ! derivative of ground heat flux wrt soil temp
    canopy%fev    = 0.0   ! latent heat flux from vegetation (W/m2)
    canopy%fes    = 0.0   ! latent heat flux from soil (W/m2)
    canopy%fhs    = 0.0   ! sensible heat flux from soil (W/m2)
	canopy%fhv    = 0.0
	canopy%fh    = 0.0
	canopy%fe    = 0.0

    ssoil%albsoilsn(:,1) = 0.1 ! soil+snow albedo
    ssoil%albsoilsn(:,2) = 0.3
    ssoil%albsoilsn(:,3) = 0.05
	ssoil%owetfac(:) = 0.0  ! edit vh 26/05/09
	canopy%fns    = 0.0 ! edit vh 26/05/09

	 bal%precip_tot = 0.0
     bal%rnoff_tot = 0.0
     bal%evap_tot = 0.0
     bal%WbalSum = 0.0
     bal%EbalSum = 0.0
     bal%RadbalSum = 0.0
     bal%canopy_drybal = 0.0
     bal%canopy_wetbal = 0.0
   

	do k=1,ntile
		! Set initial soil temperature and moisture:
		ssoil%tgg(k:mp_patch:ntile,1) = Tsoil1 + 273.15 ! soil temperature, 6 layers (K)
		ssoil%tgg(k:mp_patch:ntile,2) = Tsoil2 + 273.15
		ssoil%tgg(k:mp_patch:ntile,3) = Tsoil3 + 273.15
		ssoil%tgg(k:mp_patch:ntile,4) = Tsoil4 + 273.15
		ssoil%tgg(k:mp_patch:ntile,5) = Tsoil5 + 273.15
		ssoil%tgg(k:mp_patch:ntile,6) = Tsoil6 + 273.15
		ssoil%tgg(k:mp_patch:ntile,7) = Tsoil7 + 273.15
		ssoil%tgg(k:mp_patch:ntile,8) = Tsoil8 + 273.15
		ssoil%tgg(k:mp_patch:ntile,9) = Tsoil9 + 273.15
		ssoil%tgg(k:mp_patch:ntile,10) = Tsoil10 + 273.15
		!ssoil%tgg(k:mp_patch:ntile,11) = Tsoil10 + 273.15

		ssoil%wb(k:mp_patch:ntile,1) = Wrel1*soil%ssat(k:mp_patch:ntile) ! volumetric soil moisture,6 layers
		ssoil%wb(k:mp_patch:ntile,2) = Wrel2*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,3) = Wrel3*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,4) = Wrel4*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,5) = Wrel5*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,6) = Wrel6*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,7) = Wrel7*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,8) = Wrel8*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,9) = Wrel9*soil%ssat(k:mp_patch:ntile)
		ssoil%wb(k:mp_patch:ntile,10) = Wrel10*soil%ssat(k:mp_patch:ntile)
!		ssoil%wb(k:mp_patch:ntile,11) = Wrel10*soil%ssat(k:mp_patch:ntile)
		ssoil%S(k:mp_patch:ntile,1) = Wrel1 ! soil moisture content, relative to saturated value,6 layers
		ssoil%S(k:mp_patch:ntile,2) = Wrel2
		ssoil%S(k:mp_patch:ntile,3) = Wrel3
		ssoil%S(k:mp_patch:ntile,4) = Wrel4
		ssoil%S(k:mp_patch:ntile,5) = Wrel5
		ssoil%S(k:mp_patch:ntile,6) = Wrel6
		ssoil%S(k:mp_patch:ntile,7) = Wrel7
		ssoil%S(k:mp_patch:ntile,8) = Wrel8
		ssoil%S(k:mp_patch:ntile,9) = Wrel9
		ssoil%S(k:mp_patch:ntile,10) = Wrel10
		!ssoil%S(k:mp_patch:ntile,11) = WRel10
		ssoil%Tsurface(k:mp_patch:ntile) = Tsoil1	
		ssoil%h0(k:mp_patch:ntile) = 0.0			
		ssoil%TL(k:mp_patch:ntile) = Tsoil1	
		ssoil%SL(k:mp_patch:ntile) = 0.5
		ssoil%Tsurface(k:mp_patch:ntile) = Tsoil1
		ssoil%tss(k:mp_patch:ntile) = Tsoil1 + 273.15

		ssoil%albsoilsn(k:mp_patch:ntile,1) = max(albsoilVIS,0.02) ! soil+snow albedo
		ssoil%albsoilsn(k:mp_patch:ntile,2) = max(albsoilNIR,0.02)
	enddo


 endif

		  if (HourlyOutputFlag.eq.2) then
		  filename%out = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
					trim(RunID) // '/OutputTS/' // 'out_cable.nc'
		  
		    time_coord = 'LOC' ! set to LOC and the other option is GMT  
		  ! timeunits should have this string "seconds since 2002-01-01 00:01:00"
		  ! replacing the current date and time values 

		  TTDate   = StartDate + nint(Time0)	! current date as dmytype (use date arithmetic)
	      iday = TTDate%day
		  imonth = TTDate%month
		  iyear = TTDate%year

		  timestring = '                               '
		  WRITE(timeunits,*) timestring
		  timestring = 'seconds since'          
		  WRITE(timeunits(1:13),'(A13)')  timestring
		  WRITE(timeunits(15:18),'(I4)') TTDate%year
		  IF(iyear < 10) THEN  ! if the year is one digit concatenate with 0
			 timestring = '000' // timeunits(18:18)
			 WRITE(timeunits(15:18),'(A4)') timestring
		  ELSE IF(iyear > 10 .AND. iyear < 100) THEN
			 timestring = '00' // timeunits(17:18)
			 WRITE(timeunits(15:18),'(A4)') timestring
		  ELSE IF(iyear > 100 .AND. iyear < 1000) THEN
			 timestring = '0' // timeunits(16:18)
			 WRITE(timeunits(15:18),'(A4)') timestring
		  END IF
		  timestring = '-'
		  WRITE(timeunits(19:19),'(A1)')  timestring
		  WRITE(timeunits(20:21),'(I2)') imonth
		  IF(imonth < 10) THEN  ! if the month is one digit concatenate with 0
			 timestring = '0' // timeunits(21:21)
			 WRITE(timeunits(20:21),'(A2)') timestring
		  END IF
		  timestring = '-'
		  WRITE(timeunits(22:22),'(A1)')  timestring
		  timestring = '01'
		  WRITE(timeunits(23:24),'(I2)') iday
		  IF(iday < 10) THEN  ! if the day is one digit concatenate with 0
			 timestring = '0' // timeunits(24:24)
			 WRITE(timeunits(23:24),'(A2)') timestring
		  END IF

		  timestring = '00:00:00'
		  WRITE(timeunits(26:33),'(A8)')  timestring
		  
		  ! dimensions for netcdf output
		  max_vegpatches = ntile
		  xdimsize = mp
		  ydimsize = 1
		  output%grid='land'
		  output%patch = .FALSE.
		  output%averaging = 'all'
		  output%SWDown = .TRUE.
		  output%met = .TRUE. 
		  output%flux = .TRUE.  ! convective, runoff, NEE fluxes to output?
		  output%soil = .TRUE.  ! soil states written to output?
		  output%snow = .FALSE.  ! snow states written to output?
		  output%radiation = .TRUE.  ! net rad, albedo written to output?
		  output%carbon    = .False.  ! NEE, GPP, NPP, stores written to output?
		  output%GPP =.TRUE.
		  output%veg       = .TRUE.  ! vegetation states written to output?
		  output%params    = .TRUE.  ! input parameters used to produce run written to output?
		  output%balances  = .FALSE.  ! energy and water balances written to output?
		  output%Wbal = .TRUE.
		  output%EbalSum = .TRUE.
		  check%ranges     = .FALSE.  ! variable ranges, input and output
          check%energy_bal = .TRUE.  ! energy balance
          check%mass_bal   = .TRUE.  ! water/mass balance


		  ALLOCATE(landpt(mp))
		  ALLOCATE(lat_all(mp,1))
		  ALLOCATE(lon_all(mp,1))
		  lat_all(:,1)=LatDeg
		  lon_all(:,1) = LongDeg
		  landpt%nap = ntile
		  landpt%cstart = tile_index(:,1)
		  landpt%cend = tile_index(:,2)

		  CALL open_output_file(dels,soil,veg,bgc,rough)
	   endif

END SUBROUTINE InitCABLE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE WaterDynStep

!*******************************************************************************

SUBROUTINE InitStores (XX0)
!-------------------------------------------------------------------------------
! Initialise XX0 (because it has global scope in host program)
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Global variables
real(sp),intent(out):: XX0(:,:)     ! XX0(np,nxx) = state variables at T0
! Local variables (only needed here)
integer(i4b)  :: ixx
character(200):: GenFileName
!-------------------------------------------------------------------------------

! * INITIALISE XX0
!   * If XXInitFileFlag=0, set XX0 to spatially uniform value from ctl file
!   * If XXInitFileFlag=1, use spatial maps to set XX0:
!     maps are read from files (typically from a prior run) in \InitStores\
 do ixx=1,nxx
  if (XXInitFileFlag == 0) then
    XX0(:,ixx) = XXinitCtl(ixx)
  else
    GenFileName = trim(DataOutRootDir) // trim(DomName)              & 
	              // '/InitStores/' // trim(XXInitFileName(ixx))
    if (DiagsFlag ==1) then
      write(*,"('Open:',a)") trim(GenFileName)
    end if
	INQUIRE (file=GenFileName, exist=GenFileExistFlag)
    if (GenFileExistFlag) then
		OPEN (200, file=GenFileName, form='binary', status='old', action='read')
		REWIND (200)
		read (200) XX0(:,ixx)
		CLOSE (200)
	else
		XX0(:,ixx) = XXinitCtl(ixx)
	endif
  end if
end do
if (DiagsFlag == 1) then
  write(*,"(/,'Checks on XX0 maps:',/,                              &
  &  '    XXName               min(XX0)     max(XX0)   mean(XX0) count(XX0=Err)')")
  do ixx=1,nxx 
    write(*,"(a20,3f12.4,i10)")                                     &
            XXName(ixx,1), minval(XX0(:,ixx)), maxval(XX0(:,ixx)),  &
            sum(XX0(:,ixx))/size(XX0(:,ixx)), count(XX0(:,ixx) < ErrVal+1.0)
  end do
  do ixx=1,nxx
    write(*,"('ixx=',i1,'  XX0(1,ixx)=',4e10.3)") ixx, XX0(1,ixx)
  end do
end if

END SUBROUTINE InitStores
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE InitVV(VV)
!-------------------------------------------------------------------------------
! Initialise VV (because it has global scope in host program)
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Global variables
real(sp),intent(out):: VV(:,:)     ! VV(np,nvv) = variable parameters

!-------------------------------------------------------------------------------
! Initialise quantities with local scope to WaterDynStep:
! * Allocate arrays with dimension np or nc 
! * Initialise spatial arrays
! * Called from WaterDynStep on first call only
!-------------------------------------------------------------------------------
! Local variables (only needed here)
integer(i4b)  :: ivv, ixx, imm, ic, i, iunit
real(sp)      :: VVtest                  ! test element for reading VV
integer(i4b)  :: nVVtest(nvv)            ! test number of elements in VV data
logical(lgt)  :: nVVtestFlag(nvv)        ! flag for nVVtest = np
character(6)  :: CatchName
integer(i4b)  :: iday, iSkipDays
type(dmydate) :: FirstMetDate, MetDate, DischargeDate
integer(i4b)  :: ObsDay, ObsMonth, ObsYear
character(1)  :: YN
character(100):: dumchar
real(sp)      :: dum(20), CatchIDdum
!-------------------------------------------------------------------------------



! * READ AND ASSIGN SPATIALLY VARIABLE PARAMETERS (VV)
 do ivv = 1,nvv
  if (VVSpatialFlag(ivv) == 0) then
    VV(:,ivv) = VVConst(ivv)
  else if (VVSpatialFlag(ivv) == 1) then
    GenFileName = trim(DataInRootDir) //  trim(DomName) // '/'  &
                  // 'params/' // trim(VVFileName(ivv))
    INQUIRE (file=GenFileName, exist=GenFileExistFlag)
    if (GenFileExistFlag) then
      OPEN (unit=200, file=GenFileName, status='old', form='binary', action='read')
      REWIND (unit=200)
    else
      write(*,"('File not found:',/,a)") GenFileName
      STOP "InitModel: file not found"
    end if
    read (200) VV(:,ivv)
    CLOSE (unit=200)
  else
    write(*,"('InitModel: illegal VVSpatialFlag at ivv =',i4)")
    STOP
  end if
end do
!   * Write VV diagnostics
if (DiagsFlag == 1) then
  write(*,"(/,'Checks on VV maps:',/,                           &
  & '    VVName       VVSpatialFlag     min(VV)     max(VV)   mean(VV) count(VV=Err)')")
  do ivv=1,nvv 
    write(*,"(a20,i10,3f12.4,i10)")                             &
            VVName(ivv,1), VVSpatialFlag(ivv),                  &
                  minval(VV(:,ivv)), maxval(VV(:,ivv)),         &
                  sum(VV(:,ivv))/size(VV(:,ivv)),               &
                  count(VV(:,ivv) < ErrVal+1.0)
  end do
end if


END Subroutine InitVV

!*******************************************************************************

SUBROUTINE ReadControlFile
!-------------------------------------------------------------------------------
! * Open control file, read control data
! * All variables in ctl file are declared at head of WaterDynModule. They are
!   global to model and available to host program through USE WaterDynModule.
! * This routine is called from host program.
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
USE PointerModule
!#ifdef WIN32
USE DFLIB
!#endif


!DEC$ IF DEFINED (__INTEL_COMPILER)
! GETARG is intrinsic in Intel Fortran
!DEC$ ELSE
!use dflib, only: nargs
!DEC$ ENDIF 

implicit none
! Local variables
real(sp)    :: Dum
integer(i4b):: i, ixx, ivv, it
integer(i4b):: StartD_DailyOutput, StartM_DailyOutput, StartY_DailyOutput
integer(i4b):: EndD_DailyOutput, EndM_DailyOutput, EndY_DailyOutput
integer(i4b):: StartD_spin, StartM_spin, StartY_spin
integer(i4b):: EndD_spin, EndM_spin, EndY_spin
integer(i4b):: StartD, StartM, StartY
integer(i4b):: EndD, EndM, EndY
real(sp)     :: VVtest                  ! test element for reading VV
integer(i4b) :: nVVtest(nvv)            ! test number of elements in VV data
real(sp),allocatable:: CatchMapTemp(:)  ! temporary catchment ID np-array
logical(lgt) :: nVVtestFlag(nvv)        ! flag for nVVtest = np
integer(i2b), parameter :: Arg1 = 1  ! Specifies the first command-line argument (for GetArg)
!-------------------------------------------------------------------------------

if (Nargs().gt.1) then 
  CALL GetArg (Arg1, CtlFilename)
  write (*,*) 'Ctl File for this run:'
  write (*,*)  trim(CtlFilename)
else
	!write(*,"('Give ctl filename (in HomeDir; no .ctl extension:',$)")
	!read(*,*) CtlFileName
	!CtlFileName= 'Otway_13a'
	!CtlFileName= 'TumbKy_25a'
	!CtlFileName= 'Adelong_24'
	!CtlFileName= 'CD25ctlFile'
	CtlFileName= 'OzFlux_ext_0010_29a'
	!CtlFileName= 'Oz_0289_1965_29a'
	!CtlFileName= 'RegionPoint_29_VR_NPP_OzFlux_VV_RECCAP_25a'
endif
! * OPEN, READ AND CLOSE CTL FILE
OPEN (unit=11, file=trim(CtlFileName)//'.ctl', status='old', action='read')
REWIND (unit=11)
CALL ComSkp(11)
read (11,*) DomName                     ! domain name (eg Adelong)
read (11,*) RunID                       ! 3-char run ID (eg 05a)
read (11,*) Comment                     ! comment (second title line)
read (11,*) DataInRootDir               ! input  data root directory (with final /)
read (11,*) DataOutRootDir              ! output data root directory (with final /)
do i=1,nmm
  read (11,*) MetFileName(i)            ! Met filenames
end do
read (11,*) CO2FileName                ! CO2 filename
read (11,*) NDVIFileName                ! NDVI filename
read (11,*) LSTFileName                 ! LST filename
read (11,*) LSTTimeFileName             ! LSTTime filename
read (11,*) LSTAngleFileName            ! LSTAngle filename
read(11,*)  phiRnethFileName               ! phiHh filename
read(11,*)  phiHhFileName               ! phiHh filename
read(11,*)  phiEhFileName               ! phiEh filename
read(11,*)  phiNEEhFileName               ! phiNEEh filename
read(11,*)  phiGhFileName               ! phiGh filename
read(11,*)  TsoilhFileName               ! Tsoilh filename
read(11,*)  phiHhTimeFileName		    ! phiHhTime filename
read(11,*)  phiRnetFileName             ! phiRnet filename
read(11,*)  phiHFileName                ! phiH filename
read(11,*)  phiEFileName                ! phiE filename
read(11,*)  phiNEEFileName              ! phiNEE filename
read(11,*)  sm0008FileName              ! sm0008 filename
read(11,*)  sm0090FileName              ! sm0090 filename
read(11,*)  sm0030FileName              ! sm0030 filename
read(11,*)  sm3060FileName              ! sm3060 filename
read(11,*)  sm6090FileName              ! sm6090 filename
read(11,*)  LocalMetFileName            ! Local met filename
do i=1,nrr
  read (11,*) RRFileName(i)            ! remote sensing driver filenames
end do
read (11,*) ForwardModelFlag            ! flag to select forward model
read (11,*) CableParamFlag            ! flag to select best parameter estimates or CABLE defaults
read (11,*) SensTestFlag                ! flag for param sensitivity test
read (11,*) DiagsFlag                   ! diagnostics flag
read (11,*) ipDiag                      ! point (diagnostic) TS output flag
read (11,*) DomTSOutputFlag             ! domain-av time series output flag
read (11,*) MapTSOutputFlag             !  DD time series output flag
read (11,*) CatchTSOutputFlag           ! catchment time series output flag
read (11,*) MapOutputFlag               ! map output flag
 MapOutputFlag0 =  MapOutputFlag       ! save original value of map output flag
read (11,*) PESTOutputFlag              ! PEST output flag
read (11,*) WriteInitStoresFlag         ! flag to write stores at end of run
read (11,*) XXInitFileFlag              ! flag for XX initialisation from files
read (11,*) FracVegFlag                 ! flag for FracV (ext VegFPC, CLea)
read (11,*) HourlyOutputFlag            ! flag for hourly CABLE output to file
read (11,*) UseLocalMetDataFlag         ! Flag for using local met data if available
read (11,*) UseLAIFlag                  ! Flag for using LAI data if available
read (11,*) DelT                        ! CABLE time step (s)
read (11,*) StartD_DailyOutput,StartM_DailyOutput,StartY_DailyOutput        ! start date (hourly output)
read (11,*) EndD_DailyOutput,  EndM_DailyOutput,  EndY_DailyOutput          ! end date (hourly output)
read (11,*) StartD,StartM,StartY        ! start date
read (11,*) EndD,  EndM,  EndY          ! end date
read (11,*) nspin
nEnsemble = 1
CALL ComSkp(11)
if (nspin>0) then
	read (11,*) StartD_spin,StartM_spin,StartY_spin        ! start date
	read (11,*) EndD_spin,  EndM_spin,  EndY_spin         ! end date
endif
CALL ComSkp(11)
if (CableParamFlag.eq.1) then
	ntile = 2
else
	ntile = 1
endif

do i=1,ntt                              ! time variables (TT)
  read (11,*) TTstart(i), TTName(i,:)
end do
CALL ComSkp(11)
do i=1,nuu                              ! spatially uniform params (UU) 
  read (11,*) UUinitCtl(i), UUminCtl(i), UUmaxCtl(i),               &
              UUpinitCtl(i), UUqCtl(i), UUabsrel(i), UUscale(i),    &
              UUtransf(i), UUName(i,1:2), UUSensFlag(i),            &
              UUPESTFlag(i), UUKFFlag(i)           
end do
CALL ComSkp(11)
do i=1,nvv                              ! spatially variable params (VV)
  read (11,*) VVConst(i), VVName(i,1:2), VVSpatialFlag(i), VVFileName(i), VVSensFlag(i)
end do
CALL ComSkp(11)
do i=1,nmm                              ! met variables (MM)
  read (11,*) MMConst(i), MMName(i,1:2), MMMapOPFlag(i),             &
             MMFileID(i), MMmult(i), MMoffset(i)
end do
CALL ComSkp(11)
do i=1,nrr                              ! remote sensing variables (RR)
  read (11,*) RRConst(i), RRName(i,1:2), RRFlag(i), RRMapOPFlag(i) ,RRFileID(i)           
end do
CALL ComSkp(11)
do i=1,nxx                              ! initial state variables (XX): 
  read (11,*) XXinitCtl(i), XXpinitCtl(i), XXqCtl(i), XXabsrel(i),  &
              XXscale(i), XXtransf(i), XXName(i,1:2),               &
              XXMapOPFlag(i),XXTSOPFlag(i), XXInitFileName(i)
end do
CALL ComSkp(11)
do i=1,nff                              ! flux variables (FF)
  read (11,*) Dum, FFName(i,1:2), FFMapOPFlag(i),FFTSOPFlag(i), FFSensFlag(i)
end do
CALL ComSkp(11)
do i=1,ndd                              ! diagnostic variables (DD)
  read (11,*) Dum, DDName(i,1:2), DDMapOPFlag(i), DDTSOPFlag(i)
end do
CALL ComSkp(11)
do i=1,nzzP                             ! predicted point obs (ZZP)
  read (11,*) Dum, ZZPq(i), ZZPName(i,1:2), ZZPMapOPFlag(i)
end do
CALL ComSkp(11)
do i=1,naaP                             ! actual point obs (AAP)
  read (11,*) Dum, AAPr(i), AAPabsrel(i), AAPrlow(i), AAPscale(i),  & 
              AAPtransf(i), AAPflag(i), AAPdaymth(i), AAPName(i,1:2), AAPMapOPFlag(i)
end do
CALL ComSkp(11)
CALL ComSkp(11)
read(11,*) Ihour_output
CALL ComSkp(11)
do i=1,nzzPh                             ! predicted point obs (ZZPh)
  read (11,*) Dum, ZZPhq(i), ZZPhName(i,1:2), ZZPhMapOPFlag(i)
end do
CALL ComSkp(11)
do i=1,naaPh                             ! actual point obs (AAPh)
  read (11,*) Dum, AAPhr(i), AAPhabsrel(i), AAPhrlow(i), AAPhscale(i),  & 
              AAPhtransf(i), AAPhflag(i), AAPhdaymth(i), AAPhName(i,1:2), AAPhMapOPFlag(i)
end do
CALL ComSkp(11)
do i=1,nzzC                             ! predicted catchment obs (ZZC)
  read (11,*) Dum, ZZCq(i), ZZCName(i,1:2)
end do
CALL ComSkp(11)
do i=1,naaC                             ! actual catchment obs (AAC)
  read (11,*) Dum, AACr(i), AACabsrel(i), AACrlow(i), AACscale(i),  &
              AACtransf(i), AACflag(i), AACdaymth(i), AACName(i,1:2)
end do
CLOSE (unit=11)



! * CONSISTENCY CHECKS AND RESETS on dimensions and flags
!   * checks on nzz, naa
!if (nzzP /= naaP) STOP "Must have nzzP = naaP"
if (nzzC /= naaC) STOP "Must have nzzC = naaC"
!   * checks on flags
if (ForwardModelFlag  < 0 .or. ForwardModelFlag  > 3) STOP "Illegal ForwardModelFlag"
if (CableParamFlag  < 0 .or. CableParamFlag  > 2) STOP "Illegal CableParamFlag"
if (SensTestFlag      < 0 .or. SensTestFlag      > 5) STOP "Illegal SensTestFlag"
if (DiagsFlag         < 0 .or. DiagsFlag         > 1) STOP "Illegal DiagsFlag"
if (DomTSOutputFlag   < 0 .or. DomTSOutputFlag   > 4) STOP "Illegal DomTSOutputFlag"
if (MapTSOutputFlag   < 0 .or. MapTSOutputFlag   > 4) STOP "Illegal MapTSOutputFlag"
if (CatchTSOutputFlag < 0 .or. CatchTSOutputFlag > 4) STOP "Illegal CatchTSOutputFlag"
if (MapOutputFlag     < 0 .or. MapOutputFlag     > 4) STOP "Illegal MapOutputFlag"
if (PESTOutputFlag    < 0 .or. PESTOutputFlag    > 4) STOP "Illegal PESTOutputFlag"
if (XXInitFileFlag    < 0 .or. XXInitFileFlag    > 1) STOP "Illegal XXInitFileFlag"
if (any(VVSpatialFlag  < 0) .or. any(VVSpatialFlag  > 1)) STOP "Illegal VVSpatialFlag"
if (any(MMMapOPFlag    < 0) .or. any(MMMapOPFlag    > 4)) STOP "MMMapOPFlag"
if (any(XXMapOPFlag    < 0) .or. any(XXMapOPFlag    > 4)) STOP "XXMapOPFlag"
if (any(FFMapOPFlag    < 0) .or. any(FFMapOPFlag    > 4)) STOP "FFMapOPFlag"
if (any(DDMapOPFlag    < 0) .or. any(DDMapOPFlag    > 4)) STOP "DDMapOPFlag"
if (any(ZZPMapOPFlag   < 0) .or. any(ZZPMapOPFlag   > 4)) STOP "ZZPMapOPFlag"
if (any(AAPMapOPFlag   < 0) .or. any(AAPMapOPFlag   > 4)) STOP "AAPMapOPFlag"
if (any(ZZPhMapOPFlag   < 0) .or. any(ZZPhMapOPFlag   > 4)) STOP "ZZPhMapOPFlag"
if (any(AAPhMapOPFlag   < 0) .or. any(AAPhMapOPFlag   > 4)) STOP "AAPhMapOPFlag"
!   * checks on MMMult, MMOffset
if ( DiagsFlag ==1 .and. (any(MMmult /= 1.0) .or. any(MMoffset /= 0.0)) ) then
  write(*,"(/,'MM scaling (MM = MMmult*MM + MMoffset) is not identity')")
  write(*,"('MMName: = ',6(4x,a8))")   MMName(:,1)
  write(*,"('mult:   = ',6f12.4)") MMmult
  write(*,"('offset: = ',6f12.4)") MMoffset
!  write(*,"('Are you sure (y/n) ? ',$)")
!  read (*,*) YN
!  if (YN /= 'y') STOP "Reset MMmult, MMoffset"
end if

! * FIND ALLOCATABLE ARRAY DIMENSION np
!   * np = number of spatial points in domain (set as size of all VV arrays)
nVVtest(1:nvv) = nint(ErrVal)           ! initialise nVVtest = elements in VV(:,ivv)
CountVV: do ivv=1,nvv                   ! count elements in VV(:,ivv)
  if (.not.(VVSpatialFlag(ivv) == 0 .or. VVSpatialFlag(ivv) == 1))  &
    STOP "ReadControlFile: illegal VVSpatialFlag"
  if (VVSpatialFlag(ivv) == 0) cycle    ! spatial VV data not available
  GenFileName = trim(DataInRootDir) //  trim(DomName) // '/'        &
                // 'params/' // trim(VVFileName(ivv))
  INQUIRE (file=GenFileName, exist=GenFileExistFlag)
  if (GenFileExistFlag) then
    OPEN (unit=200, file=GenFileName, status='old', form='binary', action='read')
    REWIND (unit=200)
    do i = 1,1000000
      read (200,end=21) VVtest
      nVVtest(ivv) = i
    end do
    write(*,*) trim(GenFileName)
    STOP "ReadControlFile: VVtest count exceeded"
    21 continue
    CLOSE (unit=200)
  else
    write(*,*) trim(GenFileName)
    STOP "ReadControlFile: VV file not found - maybe wrong machine for ctl file!"
  end if
end do CountVV
!write(*,"('nVVtest =',/,10i7)") nVVtest
if ( all(nVVtest <= nint(ErrVal)) ) STOP "ReadControlFile: no spatial VV data available"
do ivv=1,nvv                            ! set np as first nVVtest /= ErrVal
  if (nVVtest(ivv) <= nint(ErrVal)) cycle
  np = nVVtest(ivv)
  exit
end do
do ivv=1,nvv
  if ( (nVVtest(ivv) > nint(ErrVal)) .and. (nVVtest(ivv) /= np) ) then
    STOP "ReadControlFile: not all available VV are of length np"
  end if
end do



! * FIND ALLOCATABLE ARRAY DIMENSION nc
!   * nc = number of catchments in domain
allocate (CatchMapTemp(np))             ! temporary catchment ID np-array
GenFileName = trim(DataInRootDir) //  trim(DomName) // '/'      &
              // 'params/' // trim(VVFileName(1))
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if (GenFileExistFlag) then
  OPEN (unit=200, file=GenFileName, status='old', form='binary', action='read')
  REWIND (unit=200)
else
  STOP "ReadControlFile: catchment map file not found"
end if
read (200) CatchMapTemp(:)
CLOSE (unit=200)
!   * find nc and put unique catchment IDs into iCatchID
CALL FindUniqueRegions (nint(CatchMapTemp(:)), nc, iCatchID)
deallocate (CatchMapTemp)

! * Write run information to screen (for forward and KF, but not PEST run)
if (DiagsFlag == 1 .or. nEnsemble > 1) then
  write(*,*)
  write(*,"('Domain:  ',a20)") DomName
  write(*,"('RunID:   ',a3 )") RunID
  write(*,"('Comment: ',a60)") Comment
  write(*,"('points     (np) =',i7)") np
  write(*,"('catchments (nc) =',i7)") nc
  write(*,"('iCatchID(1:nc)  =',/,(10i7))") iCatchID(1:nc)
end if

! * Assign StartDate_DailyOutput and EndDate_DailyOutput in date format
StartDate_DailyOutput = dmydate(StartD_DailyOutput, StartM_DailyOutput, StartY_DailyOutput)
EndDate_DailyOutput   = dmydate(EndD_DailyOutput, EndM_DailyOutput, EndY_DailyOutput)

! * Assign StartDate and EndDate in date format, and find nsteps, nmonths, ntime
StartDate = dmydate(StartD, StartM, StartY)
EndDate   = dmydate(EndD, EndM, EndY)
nSteps    = DayDifference (EndDate, StartDate)
nMonths = 0

StartDate_run = StartDate
EndDate_run = EndDate
nsteps_run = nsteps
nMonths_run = nMonths

do iT = 0,nsteps
 	! current date as dmytype (use date arithmetic)
	if (EndMonth(StartDate + iT)) nMonths = nMonths+1
enddo

! * Assign StartDate and EndDate in date format, and find nsteps, nmonths, ntime, for spin-up
StartDate_spin = dmydate(StartD_spin, StartM_spin, StartY_spin)
EndDate_spin   = dmydate(EndD_spin, EndM_spin, EndY_spin)
nSteps_spin    = DayDifference (EndDate_spin, StartDate_spin)
nMonths_spin = 0

do iT = 0,nsteps_spin
 	! current date as dmytype (use date arithmetic)
	if (EndMonth(StartDate_spin + iT)) nMonths_spin = nMonths_spin+1
enddo

! check for integral number of spin cycles
if ((nspin>0).and.(nsteps.ne.nsteps_spin)) then
	!if ((nsteps+1).ne.(nspin*(nsteps_spin+1)+ DayDifference(EndDate_run, EndDate_spin))) then
	if ((nsteps+1).lt.(nspin*(nsteps_spin+1))) then
	    write(*,*) nsteps+1
		write(*,*) nspin*(nsteps_spin+1)
		write(*,*) DayDifference (EndDate_run, EndDate_spin)
		STOP "total spin years greater than run length"
	endif
endif

ntime = int(24.*3600./delT)

END SUBROUTINE ReadControlFile

!*******************************************************************************

SUBROUTINE FindUniqueRegions (iRegMap, nReg, iRegID)
!-------------------------------------------------------------------------------
! Find unique regions in a (vector) map of integer region ID numbers, iRegMap(1:np).
! Return: nReg = number of unique regions, iRegID(1:nReg) = array of ID numbers.
! MRR, 09-nov-2005
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Global variables
integer(i4b),intent(in) :: iRegMap(:)
integer(i4b),intent(out):: nReg
integer(i4b),intent(out):: iRegID(:)
! Local variables
integer(i4b)  :: ip
character(200):: LocCatchName
logical(lgt)  :: DataExistFlag
!-------------------------------------------------------------------------------
nReg = 0                                    ! start with no valid regions
iRegID(:) = nint(ErrVal)                    ! initialise iRegID as all errors
do ip=1,np                                  ! set iRegID(1) to first valid region ID
  if (iRegMap(ip) <= nint(ErrVal)) cycle    ! error point: skip
  iRegID(1) = iRegMap(ip)
  exit
end do
if (size(iRegMap) == 1) then
  nReg = 1
  return
end if
do ip=1,np                                  ! test each element of iRegMap
  if (iRegMap(ip) <= nint(ErrVal)) cycle    ! error point: skip
  if (any(iRegMap(ip) == iRegID(1:nReg))  & ! found this region already, so skip
       .and. (nReg >= 1) ) cycle            !   (test only when nReg >= 1)
                                            ! test if there is a data file for catchment
  write(LocCatchName,"(i6)") iRegMap(ip)
  GenFileName = trim(DataInRootDir) // 'DischargeData.201111/'  &
                // trim(LocCatchName) // 'monrun.dat'
  INQUIRE (file=GenFileName, exist=DataExistFlag)
  if (.not.DataExistFlag) cycle             ! no monthly discharge data: cycle
  nReg = nReg + 1                           ! found a new region: augment nReg
  iRegID(nReg) = iRegMap(ip)                ! add new region to RegID
end do

END SUBROUTINE FindUniqueRegions


!*******************************************************************************

SUBROUTINE ReadMet (TTime, MM, hMM)
!-------------------------------------------------------------------------------
! * Read one scan of single-point meteorology from MetFile
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
implicit none
! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time from host program
real(sp),intent(out):: MM(:,:) , hMM(:,:,:)     ! MM(np,nmm) = met variables; hMM(np, nmm, ntime)
! Local variables
real(sp)     :: Dum(8)              ! dummy
real(sp)     :: ADisCDmm, ADisCMmm  ! catchment discharge (mm/d, mm/mth)
integer(i4b) :: itry, it, imm, ic, iunit
integer(i4b) :: MetDay, MetMonth, MetYear
integer(i4b) :: ObsDay, ObsMonth, ObsYear
type(dmydate):: TTDate              ! current date from host program
type(dmydate):: MetDate, NDVIDate, DischargeDate
type(dmydate):: LocalMetDate                                ! obs dates, as derived type in DateFunctions
integer(i4b)  :: EndLocalMetFlag=0, EndMetFlag=0
real(sp) :: rh(np), deltavph(np)
LOGICAL      :: isopen
integer(i4b) :: istat
!-------------------------------------------------------------------------------

! * set output arrays to null obs flag (values will be overwritten)
MM    = ErrVal
hMM = ErrVal

! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime, TargetMM=MM)

! * Construct date from TTime, in date format and as char string for filenames
TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
write(TTDateChar,"(i4.4,2i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay)

! * MET DATA
do imm = 1,nmm-2 ! read all met excpet vph
  iunit = 30 + imm
  INQUIRE(iunit,OPENED=isopen,IOSTAT=istat)
  if (isopen.and.(istat.eq.0)) then
	  read (iunit,end=2) MetDate, MM(:,imm)
	  if (MetDate .ne. TTDate) then
		write(*,"('imm, MetDate, TTDate:',7i6)") imm, MetDate, TTDate
		stop 'ReadMetObs: MetDate .ne. TTDate' 
	  end if
  endif
end do
! * Rescale met data as MM = MMmult*MM + MMoffset
forall (imm=1:nmm)
  MM(:,imm) = MMmult(imm)*MM(:,imm) + MMoffset(imm)
end forall

if (any(MM(:,3)>80)) then
write(*,*) 'Tmax too high'
endif

where (MM(:,3)>60.0)  ! temp fix for extreme high T in Tmax record
  MM(:,3) = 40.0
end where

do imm = nmm-1,nmm ! read or estimate vph
  iunit = 30 + imm
  INQUIRE(iunit,OPENED=isopen,IOSTAT=istat)
  if (isopen.and.(istat.eq.0)) then
	  read (iunit,end=2) MetDate, MM(:,imm)
	  if (MetDate .ne. TTDate) then
		write(*,"('imm, MetDate, TTDate:',7i6)") imm, MetDate, TTDate
		stop 'ReadMetObs: MetDate .ne. TTDate' 
	  end if
  else
	MM(:,imm)= esatf(MM(:,4))  ! set vp to saturated vp at Tmin
  end if
end do


where (MM(:,5)<2.0)  ! temp fix for -ve values in vph records
  MM(:,5) = 2.0
end where

where (MM(:,6)<2.0)  ! temp fix for -ve values in vph records
  MM(:,6) = 2.0
end where

where (MM(:,5)>140.0)  ! temp fix for v. large values in vph records
  MM(:,5) = 140.0
end where

where (MM(:,6)>140.0)  ! temp fix for v. large values in vph records
  MM(:,6) = 140.0
end where

goto 3

2 if (UseLocalMetDataFlag == 0) then
	stop 'ReadMetObs: end of daily met file'
  else
    EndMetFlag=1
	continue
  end if



! local Met data (optional)
3 if (UseLocalMetDataFlag >= 1) then
	iunit = 30 + nmm +1
	if (UseLocalMetDataFlag == 1) then
		read (iunit, end=4)  LocalMetDate, hMM(:,:,:)
	elseif (UseLocalMetDataFlag == 2) then
		read (iunit, end=4)  LocalMetDate, hMM(1,:,:)
		hMM = SPREAD(hMM(1,:,:),1,np)
	endif
	where (hmm(:,2,:)<ErrVal+1)
		hmm(:,2,:) = 335.97 * (((hmm(:,5,:)+ 273.16) / 293.0)**6)
	endwhere
	where ((hmm(:,7,:)>ErrVal+1).and.(hmm(:,7,:)<500))
		hmm(:,7,:)=hmm(:,7,:)*10.
	endwhere
	if (LocalMetDate .ne. TTDate) then
		write(*,"('LocalMetDate, TTDate:',7i6)") LocalMetDate, TTDate
		stop 'ReadMet: LocalMetDate .ne. TTDate' 
	end if
endif



goto 5

4 EndLocalMetFlag =1

5 if (EndLocalMetFlag==1 .and. EndMetFlag==1) then
		stop 'ReadMetObs: no met data available'
endif


RETURN

END SUBROUTINE ReadMet

!*******************************************************************************
SUBROUTINE ReadRR(TTime,RR)
!-------------------------------------------------------------------------------
! * Read one scan of spatial remote sensing driver from RRFile
! * Read one scan of atmospheric CO2 file
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
implicit none
! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time from host program
real(sp),intent(inout):: RR(:,:) 
! Local variables
real(sp)     :: tmpRR(size(RR,1))
integer(i4b) :: irr, iunit
type(dmydate):: TTDate              ! current date from host program
type(dmydate):: RRDate
integer(i4b) :: RRDay,RRMonth,RRyear
LOGICAL      :: isopen
integer(i4b) :: istat
real(i4b)     :: dum
!-------------------------------------------------------------------------------
! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime)
TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))


	

do irr= 1,nrr
! if time series data exist and time-series starts afer StartDate, check if time-series start-date has been reached
	if (EqDates(NextRRDatetmp(irr),TTDate).and.RRFlag(irr).eq.2) then
		RRFlag(irr)=1
		NextRRDate(irr) = NextRRDatetmp(irr)
		NextRR(:,irr) = RRtmp(:,irr)
		close(900+irr)
	endif
  
	if (EqDayofYear(NextRRDate(irr),TTDate).and.RRFlag(irr).eq.2) then !climatology
		iunit = 900 + irr
		INQUIRE(iunit,OPENED=isopen,IOSTAT=istat)
		RR(:,irr) = max(NextRR(:,irr),0.0)
		read(iunit,IOSTAT=istat)  NextRRDate(irr), NextRR(:,irr)
		if (isopen.and.(istat.eq.-1)) then
		 REWIND(iunit)
		 read(iunit,IOSTAT=istat)  NextRRDate(irr), NextRR(:,irr)
		endif
	!	write(*,"('reading climatology ',a,3i5)") trim(RRName(irr,1)) , NextRRDate(irr)%day,NextRRDate(irr)%month, &
	!						  NextRRDate(irr)%year
	elseif (EqDates(NextRRDate(irr),TTDate).and.RRFlag(irr).eq.1) then ! time series
		iunit = 9000 + irr
		INQUIRE(iunit,OPENED=isopen,IOSTAT=istat)
		RR(:,irr) = max(NextRR(:,irr),0.0)
		if (isopen.and.(istat.eq.0)) then
			read(iunit,IOSTAT=istat)  NextRRDate(irr), NextRR(:,irr)
		!	write(*,"('reading time series ',a,3i5)") trim(RRName(irr,1)) , NextRRDate(irr)%day,NextRRDate(irr)%month, &
		!					  NextRRDate(irr)%year
		endif
		if (isopen.and.(istat.eq.-1)) then
			write(*,"('end of time-series file ',a,' :revert to climatology')")  RRName(irr,1)

			! Initialise climatology file
			 iunit = 900 + irr	 
			 GenFileName = trim(DataInRootDir) //  trim(DomName) // '/'  &
					  // 'met/9999010199991231_' // trim(RRFileID(irr))//'_clim.bin'
			 INQUIRE (file=GenFileName, exist=GenFileExistFlag) 
			 if (GenFileExistFlag) then
				 OPEN (iunit, file=GenFileName, status='old', form='binary', action='read')
				 REWIND (iunit)
				 read(iunit) NextRRDate(irr), NextRR(:,irr)
				 do while(.not.EqMonthofYear(NextRRDate(irr),TTDate)) 
					read(iunit) NextRRDate(irr), NextRR(:,irr)
				 enddo
				 read(iunit,IOSTAT=istat) NextRRDate(irr), NextRR(:,irr)
				 if (isopen.and.(istat.eq.-1)) then
					REWIND(iunit)
					read(iunit,IOSTAT=istat)  NextRRDate(irr), NextRR(:,irr)
				endif
				!write(*,"('reading climatology ',a,4i5)") trim(RRName(irr,1)) , NextRRDate(irr)%day,NextRRDate(irr)%month, &
				!			  NextRRDate(irr)%year, RRFlag(irr)
				RRFlag(irr)=2
			 else
					write(*,"('File not found:',/,a)") GenFileName
					STOP ' EndDate after LastRRDate and no climatology' 
			 endif

		endif
	endif
enddo
	
		
RETURN
END SUBROUTINE ReadRR
!*******************************************************************************
SUBROUTINE ReadCO2(TTime,UU)
!-------------------------------------------------------------------------------
! * Read one scan of spatial remote sensing driver from RRFile
! * Read one scan of atmospheric CO2 file
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
implicit none
! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time from host program
real(sp),intent(inout):: UU(:) 
! Local variables
integer(i4b) :: irr, iunit
type(dmydate):: TTDate              ! current date from host program
LOGICAL      :: isopen
integer(i4b) :: istat
real(i4b)     :: dum
!-------------------------------------------------------------------------------
! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime,TargetUU=UU)
TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
!write(*,*) NextCO2Date,TTDate
if (EqDates(NextCO2Date,TTDate)) then 
	iunit = 90
	INQUIRE(iunit,OPENED=isopen,IOSTAT=istat)
	CO2A = NextCO2/1.e6
	read(iunit,IOSTAT=istat) NextCO2Date, NextCO2
	write(*,*) NextCO2Date, NextCO2
endif
		
RETURN
END SUBROUTINE ReadCO2

!*******************************************************************************
SUBROUTINE ReadObs (TTime, AAP,AAPh,AAC)
!-------------------------------------------------------------------------------

! * When available, read one scan of actual and ancillary observables
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
implicit none
! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time from host program
real(sp),intent(out):: AAP(:,:)     ! AAP(np,naaP) = actual point observations
real(sp),intent(out):: AAPh(:,:,:)     ! AAPh(np,naaP,ntime) = actual point observations
real(sp),intent(out):: AAC(:,:)     ! AAC(nc,naaC) = actual catchment observations
! Local variables
real(sp)     :: Dum(8)              ! dummy
real(sp)     :: ADisCDmm, ADisCMmm  ! catchment discharge (mm/d, mm/mth)
integer(i4b) :: itry, it, imm, ic, iunit, k, hour(np)
integer(i4b) :: MetDay, MetMonth, MetYear
integer(i4b) :: ObsDay, ObsMonth, ObsYear
type(dmydate):: TTDate              ! current date from host program
type(dmydate):: MetDate, NDVIDate, DischargeDate
                                    ! obs dates, as derived type in DateFunctions
!-------------------------------------------------------------------------------

! * set output arrays to null obs flag (values will be overwritten)
AAP   = ErrVal
AAPh = ErrVal
AAC   = ErrVal
! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime, TargetAAP=AAP, TargetAAPh=AAPh, TargetAAC=AAC)

! * Construct date from TTime, in date format and as char string for filenames
TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
write(TTDateChar,"(i4.4,2i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay)

 
! * REMOTELY SENSED DATA
!   * NDVI
if (DayDifference(NextNDVIDate, TTDate) > 0) then       ! Next NDVI is AFTER current date:
  ANDVI = ErrVal                                        !   so no data available
else if (DayDifference(NextNDVIDate, TTDate) == 0) then ! Next NDVI is AT current date
  ANDVI = NextANDVI                                     !   so data are available
  read (40,end=2) NextNDVIDate, NextANDVI               !   and read next (future) record
!  if (DiagsFlag == 1) then
!    write(*,"('NextNDVIDate,NextANDVI:',3i6,3f10.2)") NextNDVIDate,NextANDVI(1:3)
!  end if
else                                                    ! Next NDVI BEFORE current date
  STOP "ReadObs: NextNDVIDate before TTDate"         !   should not happen
end if
goto 1
2 NextNDVIDate = dmydate(1,1,2100)                      ! EOF: set next NDVI in far future
1 continue

!   * LST, LSTTime and LSTAngle (possibility of >1 per day)
  ALST     = ErrVal                                     !   so no data available
  ALSTTime = ErrVal
  ALSTAngle = ErrVal
  


do while (DayDifference(NextLSTDate, TTDate).eq.0)
   hour = floor(NextALSTTime)+1
  forall (k = 1:np, (hour(k)>=1))
	ALST(k,	hour(k)) = NextALST(k)
	ALSTTime(k,	hour(k)) = NextALSTTime(k)
	ALSTAngle(k,hour(k)) = NextALSTAngle(k)
  end forall
  read (41,end=4) NextLSTDate,     NextALST             !   and read next (future) record
  read (42,end=4) NextLSTTimeDate, NextALSTTime
  read (43,end=4) NextLSTTimeDate, NextALSTAngle
  where (NextALST < -999.+1.0)
    NextALST = ErrVal
  elsewhere
    NextALST = NextALST - 273.16                        ! degK to degC
  end where
  where (NextALSTTime < -999.+1.0)
    NextALSTTime = ErrVal
  elsewhere
    NextALSTTime = NextALSTTime*24. + SolarUTCOffset        ! UTC to solar time
  end where
  where (NextALSTAngle< -999.+1.0)
    NextALSTAngle = ErrVal
  end where
  if (NextLSTDate.ne.NextLSTTimeDate) then              ! check that LST dates match
    write(*,"('Error: NextLSTDate     =',3i6)") NextLSTDate
    write(*,"('       NextLSTTimeDate =',3i6)") NextLSTTimeDate
    STOP "ReadObs: mismatch between NextLSTDate, NextLSTTimeDate"
  end if
end do
goto 3



4 NextLSTDate     = dmydate(1,1,2100)                   ! EOF: set next LST in far future
  NextLSTTimeDate = dmydate(1,1,2100)
  NextLSTAngleDate = dmydate(1,1,2100)
3 continue

! phiRnet, phiH, phiE, phiNEE (daily daytime values)
! 1. phiRnet
if (DayDifference(NextAphiRnetDate, TTDate) > 0) then       ! Next phiRnet is AFTER current date:
  AphiRnet = ErrVal                                        !   so no data available
else if (DayDifference(NextAphiRnetDate, TTDate) == 0) then ! Next phiRnet is AT current date
  AphiRnet = NextAphiRnet                                   !   so data are available
  read (661,end=201) NextAphiRnetDate, NextAphiRnet               !   and read next (future) record
else                                                    ! Next phiRnet BEFORE current date
  STOP "ReadObs: NextAphiRnetDate before TTDate"         !   should not happen
end if
goto 101
201 NextAphiRnetDate = dmydate(1,1,2100)                      ! EOF: set next phiRnet in far future
101 continue

! 2. phiH
if (DayDifference(NextAphiHDate, TTDate) > 0) then       ! Next phiH is AFTER current date:
  AphiH = ErrVal                                        !   so no data available
else if (DayDifference(NextAphiHDate, TTDate) == 0) then ! Next phiH is AT current date
  AphiH = NextAphiH                                   !   so data are available
  read (662,end=202) NextAphiHDate, NextAphiH               !   and read next (future) record
else                                                    ! Next phiH BEFORE current date
  STOP "ReadObs: NextAphiHDate before TTDate"         !   should not happen
end if
goto 102
202 NextAphiHDate = dmydate(1,1,2100)                      ! EOF: set next phiH in far future
102 continue

! 3. phiE
if (DayDifference(NextAphiEDate, TTDate) > 0) then       ! Next phiE is AFTER current date:
  AphiE = ErrVal                                        !   so no data available
else if (DayDifference(NextAphiEDate, TTDate) == 0) then ! Next AphiE is AT current date
  AphiE = NextAphiE                                   !   so data are available
  read (663,end=203) NextAphiEDate, NextAphiE               !   and read next (future) record
else                                                    ! Next phiE BEFORE current date
  STOP "ReadObs: NextAphiEDate before TTDate"         !   should not happen
end if
goto 103
203 NextAphiEDate = dmydate(1,1,2100)                      ! EOF: set next phiE in far future
103 continue

! 4. phiNEE
if (DayDifference(NextAphiNEEDate, TTDate) > 0) then       ! Next phiNEE is AFTER current date:
  AphiNEE = ErrVal                                        !   so no data available
else if (DayDifference(NextAphiNEEDate, TTDate) == 0) then ! Next phiNEE is AT current date
  AphiNEE = NextAphiNEE                                   !   so data are available
  read (664,end=204) NextAphiNEEDate, NextAphiNEE               !   and read next (future) record

	else                                                    ! Next phiNEE BEFORE current date
  STOP "ReadObs: NextAphiNEEDate before TTDate"         !   should not happen
end if
goto 104
204 NextAphiNEEDate = dmydate(1,1,2100)                      ! EOF: set next phiNEE in far future
104 continue

! point NPP (whole of run)
GenFileName = trim(DataInRootDir) // trim(DomName) // '/obs/' // 'VR_NPP.bin'
INQUIRE (file=GenFileName, exist=GenFileExistFlag)
if(GenFileExistFlag) then
OPEN (900, file=GenFileName, form='binary', status='old', action='read')
	read(900,end=901) AphiNPP
endif
901 continue

! sm0008, sm0090, sm0030, sm3060, sm6090 (daily soil-moisture values)
! 1. sm0008

if (DayDifference(NextAsm0008Date, TTDate) > 0) then       ! Next sm0008 is AFTER current date:
  Asm0008 = ErrVal                                        !   so no data available
else if (DayDifference(NextAsm0008Date, TTDate) == 0) then ! Next sm0008 is AT current date
  Asm0008 = NextAsm0008                                   !   so data are available
  read (771,end=701) NextAsm0008Date, NextAsm0008               !   and read next (future) record
  where (NextAsm0008  <-99.+1.0)
    NextAsm0008  = ErrVal
  end where
else                                                    ! Next sm0008 BEFORE current date
  STOP "ReadObs: NextAsm0008Date before TTDate"         !   should not happen
end if
goto 711
701 NextAsm0008Date = dmydate(1,1,2100)                      ! EOF: set next sm0008 in far future
711 continue

! 2. sm0090
if (DayDifference(NextAsm0090Date, TTDate) > 0) then       ! Next sm0090 is AFTER current date:
  Asm0090 = ErrVal                                        !   so no data available
else if (DayDifference(NextAsm0090Date, TTDate) == 0) then ! Next sm0090 is AT current date
  Asm0090 = NextAsm0090                                   !   so data are available
  read (772,end=702) NextAsm0090Date, NextAsm0090               !   and read next (future) record
  where (NextAsm0090  <-99.+1.0)
    NextAsm0090  = ErrVal
  end where
else                                                    ! Next sm0090 BEFORE current date
  STOP "ReadObs: NextAsm0090Date before TTDate"         !   should not happen
end if
goto 712
702 NextAsm0090Date = dmydate(1,1,2100)                      ! EOF: set next sm0090 in far future
712 continue

! 3. sm0030
if (DayDifference(NextAsm0030Date, TTDate) > 0) then       ! Next sm0030 is AFTER current date:
  Asm0030 = ErrVal                                        !   so no data available
else if (DayDifference(NextAsm0030Date, TTDate) == 0) then ! Next sm0030 is AT current date
  Asm0030 = NextAsm0030                                   !   so data are available
  read (773,end=703) NextAsm0030Date, NextAsm0030               !   and read next (future) record
  where (NextAsm0030  <-99.+1.0)
    NextAsm0030  = ErrVal
  end where
else                                                    ! Next sm0030 BEFORE current date
  STOP "ReadObs: NextAsm0030Date before TTDate"         !   should not happen
end if
goto 713
703 NextAsm0030Date = dmydate(1,1,2100)                      ! EOF: set next sm0030 in far future
713 continue

! 4. sm3060
if (DayDifference(NextAsm3060Date, TTDate) > 0) then       ! Next sm3060 is AFTER current date:
  Asm3060 = ErrVal                                        !   so no data available
else if (DayDifference(NextAsm3060Date, TTDate) == 0) then ! Next sm3060 is AT current date
  Asm3060 = NextAsm3060                                   !   so data are available
  read (774,end=704) NextAsm3060Date, NextAsm3060               !   and read next (future) record
  where (NextAsm3060  <-99.+1.0)
    NextAsm3060  = ErrVal
  end where
else                                                    ! Next sm3060 BEFORE current date
  STOP "ReadObs: NextAsm3060Date before TTDate"         !   should not happen
end if
goto 714
704 NextAsm3060Date = dmydate(1,1,2100)                      ! EOF: set next sm3060 in far future
714 continue

! 5. sm6090
if (DayDifference(NextAsm6090Date, TTDate) > 0) then       ! Next sm6090 is AFTER current date:
  Asm6090 = ErrVal                                        !   so no data available
else if (DayDifference(NextAsm6090Date, TTDate) == 0) then ! Next sm6090 is AT current date
  Asm6090 = NextAsm6090                                   !   so data are available
  read (775,end=705) NextAsm6090Date, NextAsm6090               !   and read next (future) record
  where (NextAsm6090  <-99.+1.0)
    NextAsm6090  = ErrVal
  end where
else                                                    ! Next sm6090 BEFORE current date
  STOP "ReadObs: NextAsm6090Date before TTDate"         !   should not happen
end if
goto 715
705 NextAsm6090Date = dmydate(1,1,2100)                      ! EOF: set next sm6090 in far future
715 continue



!**************************************************************************


!   *Hourly fluxes phiRneth, phiHh, phiEh, phiNEEh  (possibility of >1 per day)
  AphiRneth    = ErrVal
  AphiHh       = ErrVal
  AphiEh       = ErrVal
  AphiNEEh     = ErrVal  
  AphiGh       = ErrVal   
  ATsoilh      = ErrVal                                
  AphiHhTime   = ErrVal

do while (DayDifference(NextAphiHhDate, TTDate).eq.0)
   
   hour = nint(NextAphiHhTime)


  forall (k = 1:np, (hour(k)>=1))
    AphiRneth(k,hour(k))  = NextAphiRneth(k)
    AphiHh(k,hour(k))     = NextAphiHh(k)
    AphiEh(k,hour(k))     = NextAphiEh(k)
	AphiNEEh(k,hour(k))   = NextAphiNEEh(k)
	AphiGh(k,hour(k))     = NextAphiGh(k)
	ATsoilh(k,hour(k))     = NextATsoilh(k)
	AphiHhTime(k,hour(k)) = NextAphiHhTime(k)
  end forall

  read (551,end=44) NextAphiRnethDate,    NextAphiRneth             !   and read next (future) record
  read (552,end=44) NextAphiHhDate,       NextAphiHh
  read (553,end=44) NextAphiEhDate,       NextAphiEh
  read (554,end=44) NextAphiNEEhDate,     NextAphiNEEh
  read (555,end=44) NextAphiGhDate,       NextAphiGh
  read (557,end=44) NextATsoilhDate,       NextATsoilh
  read (556,end=44) NextAphiHhTimeDate,   NextAphiHhTime

  where (NextAphiRneth  <ErrVal+1.0)
    NextaphiRneth  = ErrVal
  end where
  where (NextAphiHh  <ErrVal+1.0)
    NextaphiHh  = ErrVal
  end where
  where (NextAphiEh  <ErrVal+1.0)
    NextaphiEh  = ErrVal
  end where
  where (NextAphiNEEh  <ErrVal+1.0)
    NextaphiNEEh  = ErrVal
  end where
  where (NextAphiGh  <ErrVal+1.0)
    NextAphiGh  = ErrVal
  end where
  where (NextATsoilh  <ErrVal+1.0)
    NextATsoilh  = ErrVal
  end where
  where (NextAphiHhTime<ErrVal+1.0)
    NextAphiHhTime= ErrVal
  end where

  if (NextAphiHhDate.ne.NextAphiHhTimeDate) then              ! check that flux dates match
    write(*,"('Error: NextAphiHhDate     =',3i6)") NextAphiHhDate
    write(*,"('       NextphiHhATimeDate =',3i6)") NextAphiHhTimeDate
    STOP "ReadObs: mismatch between NextAphiHhDate , NextAphiHhTimeDate"
  end if
end do
goto 34



44 NextAphiHhDate     = dmydate(1,1,2100)                   ! EOF: set next flux in far future
   NextAphiHhTimeDate = dmydate(1,1,2100)

34 continue


! * MONTHLY DISCHARGE DATA
!   * file format (i6,i2,2f12.2): Year, Month, Discharge(mm/mth) (recorded, simulated)
!   * read discharge obs at end of month and when available
!   * ADisCM is mean monthly actual discharge [m/mth] on last day of month,
!     and ErrVal on all other days
do ic=1,nc
  CMDataExist: if (DisCMDataExist(ic)) then
    EndOfMonth: if (EndMonth(TTDate)) then
      read (2205+ic,205,end=6) ObsYear, ObsMonth, ADisCMmm
      205 format (i6,i2,f12.2)                                      ! sep-2005 data (monthly)
!      read (205+ic,*,end=6) ObsYear, ObsMonth, ADisCMmm             ! oct-2006 data (monthly)
      DischargeDate = dmydate(28, ObsMonth, ObsYear)
      DischargeDate = dmydate(DaysInMonth(DischargeDate), ObsMonth, ObsYear)
      if (ADisCMmm >= 0.0) then     ! rescale discharge data to m/mth
        ADisCM(ic) = ADisCMmm * 0.001
      else
        ADisCM(ic) = ErrVal         ! if obs < 0, set missing data flag
      end if
    else
      ADisCM(ic) = ErrVal           ! if not end month, set ADisCM as missing
    end if EndOfMonth
  end if CMDataExist
  goto 5
  6 DisCMDataExist(ic) =.false.     ! EOF: set data exist flag to false
  5 continue
end do

! * DAILY DISCHARGE DATA
!   * file format (i6,2i2,f10.3): Year, Month, Day, Discharge(mm/d)
do ic=1,nc
  CDDataExist: if (DisCDDataExist(ic)) then
    read (1205+ic,1205,end=8) ObsYear, ObsMonth, ObsDay, ADisCDmm
    1205 format (i6,2i2,f10.3)                                      ! sep-2005 data (daily)
!    read (1205+ic,*,end=8) ObsYear, ObsMonth, ObsDay, ADisCDmm      ! oct-2006 data (daily)
    DischargeDate = dmydate(ObsDay, ObsMonth, ObsYear)
      if (ADisCDmm >= 0.0) then     ! rescale discharge data to m/day
        ADisCD(ic) = ADisCDmm * 0.001
      else
        ADisCD(ic) = ErrVal         ! if obs < 0, set missing data flag
      end if
  
  end if CDDataExist
  goto 7
  8 DisCDDataExist(ic) =.false.     ! EOF: set data exist flag to false
  7 continue
end do

RETURN

END SUBROUTINE ReadObs
!*******************************************************************************


Subroutine GetSubdiurnalMet(Ttime,MM,MMprev,MMnext,VV,hMM)

USE TypeDef
USE Constants
USE DateFunctions
USE Utils
USE PointerModule
USE SubDiurnalMetModule

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

!*******************************************************************************
!*******************************************************************************


SUBROUTINE WaterDynFlux (TTime, XX, MM,RR,UU, VV,   FF, DD)
!-------------------------------------------------------------------------------
! * Calculate fluxes (FF), including dXXdt as first nxx entries of FF(np,nff)
! * Energy fluxes PhiXy [W/m2] are for entity X (= S,L,A,H,E) from surface 
!   component y (= v,s), averaged over daylight hours (DayltFrac*SecDay).
! * Mass fluxes FXppp [molX/m2/d] are for entity X and process ppp.
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
USE Utils
USE PointerModule
implicit none
! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time variables
real(sp),intent(in) :: XX(:,:)      ! XX(np,nxx) = state variables
real(sp),intent(in) :: MM(:,:)      ! MM(np,nmm) = met variables
real(sp),intent(in) :: RR(:,:)      ! RR(np,nrr) = met variables
real(sp),intent(in) :: UU(:)        ! UU(nuu)    = spatially uniform parameters
real(sp),intent(in) :: VV(:,:)      ! VV(np,nvv) = spatially variable parameters
real(sp),intent(out):: FF(:,:)      ! FF(np,nff) = fluxes and dXXdt
real(sp),intent(out):: DD(:,:)      ! DD(np,ndd) = derived quantities
! Local variables
real(sp):: DecDeg                   ! solar declination
real(sp):: SolarNorm(size(XX,1))    ! TOA solar frac (dummy only)
real(sp):: YearFrac                 ! year fraction                     [-]
real(sp):: Pmb(size(XX,1))          ! air pressure                      [mb]
real(sp):: TempA(size(XX,1))        ! average daytime air temp          [deg C]
real(sp):: TempS(size(XX,1))        ! average daytime surface temp      [deg C]
real(sp):: DefA(size(XX,1))         ! average daytime air deficit       [molW/molA]
real(sp):: DefS(size(XX,1))         ! average daytime surface deficit   [molW/molA]
real(sp):: CO2C(size(XX,1))         ! compensation point [CO2]          [molC/molA]
real(sp):: alfaWdef(size(XX,1))     ! WUE from deficit                  [molC/molW]
real(sp):: alfaW(size(XX,1))        ! WUE: actual                       [molC/molW]
real(sp):: PPc(size(XX,1))          ! PPc = Ga/(Ga+Grad) in CE          [-]
real(sp):: EpsA(size(XX,1))         ! dimensionless dQsat/dT            [-] 
real(sp):: PhiSd(size(XX,1))        ! daylight down solar flux          [W/m2]
real(sp):: PhiLd(size(XX,1))        ! daylight down thermal flux        [W/m2]
real(sp):: PhiAi(size(XX,1))        ! daylight isothermal avail energy  [W/m2]
real(sp):: PhiEq(size(XX,1))        ! equil latent heat flux            [W/m2]
real(sp):: FWTraE(size(XX,1)),  FWTraW(size(XX,1))  ! E-lim, W-lim FWTra    [m/d]
real(sp):: FWTr1W(size(XX,1)),  FWTr2W(size(XX,1))  ! layer W-lim FWTra     [m/d]
real(sp):: WRel1B(size(XX,1)),  WRel2B(size(XX,1))  ! bounded WRel
real(sp):: WRel1XS(size(XX,1)), WRel2XS(size(XX,1)) ! excess WRel for saturation checks
real(sp):: FQ(size(XX,1))                           ! incident PAR flux [molQ/m2/d]
real(sp):: WRelAv(size(XX,1))                       ! for leaf allocation coefficient
real(sp):: FracV(size(XX,1)), FracVExt1(size(XX,1)) ! vegetation fraction cover
real(sp):: FracVmax                                 ! maximum FracV
type(dmydate):: TTDate                              ! current date
! Parameters taken as constant in WaterDynFlux
real(sp),parameter:: Emis = 1.0                     ! emissivity        [-]
real(sp),parameter:: Pmb0 = 1000.0                  ! air pressure      [mb]
real(sp),parameter:: WRelLow = 0.01                 ! lowest WRel (to stop WRel < 0)
!-------------------------------------------------------------------------------

! * Reference output arrays (values will be overwritten)
FF = ErrVal
DD = ErrVal
! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime, TargetXX=XX, TargetMM=MM,TargetRR=RR,     &
               TargetUU=UU, TargetVV=VV,   TargetFF=FF, TargetDD=DD)
! * Find TTDate, YearFrac, DayLtFrac
TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
YearFrac = mod(YearDay(TTDate)/365.25, 1.0)         ! mod handles 31 Dec in leap year
CALL AstronDailySolar (YearFrac, LatDeg,   DecDeg, DayLtFrac, SolarNorm)

! * Find LAI and veg cover fraction, related by CoeffBeer. Pre-apply rLAIMult to LAI.
FracVmax  = 1.0 - exp(-CoeffBeer*rLAImax)
!   * FracV from external VegFPC
FracVExt1 = coeffPAR * FAPAR(:,nint(TTMonth))       ! externally specified FracV
where (FracVExt1 < FracVmax)
  rLAIExt = -log(1.0-FracVExt1) / CoeffBeer
elsewhere
  rLAIExt = rLAImax
end where
FracVExt = 1.0 - exp(-CoeffBeer*rLAIExt)            ! and calculate consistent FracV
!   * FracV from internal CLea (with CLea0 = RhoCLeaf*LeafThick)
rLAICLea  = min(max(CLea/CLea0, 0.0), rLAImax)      ! constrain to 0 <= rLAI <= rLAImax
FracVCLea = 1.0 - exp(-CoeffBeer*rLAICLea)          ! and calculate consistent FracV
!   * select FracV for computation
if (FracVegFlag <= 0) then
  FracV = FracVExt
else
  FracV = FracVCLea
end if

! * Find average daytime air temp, pressure, density, sat deficit
TempA  = 0.75*TempMax + 0.25*TempMin                ! average daytime temp
Pmb    = Pmb0 * exp( -(Grav*Altitude*RMA)   &       ! air pressure from altimeter eq
         / (Rgas*(273.16+TempA)) )                  ! (McIntosh and Thom p. 10)
RhoA   = (100.0*Pmb)/(Rgas*(273.16+TempA))          ! air density [molA/m3]
DefA   = (Esatf(TempMax) - Esatf(TempMin)) / Pmb    ! sat deficit [molW/molA]
DefA   = max(DefA, 0.001)                           ! prohibit supersaturation                      

! -------------
! ENERGY FLUXES
! -------------
! * Find irradiances, available energy, equilibrium latent heat flux
PPc    = Gaero / ( Gaero + 4.0*SBoltz*((TempA+273.16)**3)/(RhoA*Capp) )
                                                    ! PPc = Ga/(Ga+Gr)      [-]
EpsA   = Epsif(TempA, Pmb)                          ! Epsi at TempA         [-]
PhiSd  = SolarMJ * 1.0e6 / (DayltFrac*SecDay)       ! daylt down solar      [W/m2]
PhiLd  = 335.97 * (((TempA + 273.16) / 293.0)**6)   ! daylt down thermal    [W/m2]
                                                    !   (Swinbank formula)
PhiAi  = (1.0-Albedo)*PhiSd + Emis *    &           ! daylt iso-avail engy  [W/m2]
         (PhiLd - SBoltz*((TempA + 273.16)**4))     !   (veg + soil)
PhiEq  = PhiAi * (PPc*EpsA) / (PPc*EpsA + 1.0)      ! equil ltnt heat flux  [W/m2]
PhiEq  = max(PhiEq, 1.0)                            ! PhiEq > +1 W/m2, so non-negative
! * Find precip, PT evaporation
FWPrec = Precip                                     ! precipitation         [m/day]
FWPT   = CoeffPT * PhiEq * ((DayltFrac*SecDay)  & 
         / (RhoW*Rlat))                             ! PT (tot E-lim) evap   [m/day]

! ------------------
! EVAPORATIVE FLUXES
! ------------------
! * Find bounded initial WRel for calculation of FWSoil, FWTr1W, FWTr2W, FWLchA, FWLchB
WRel1B = min(max((WRel1-WRelLow)/(1.0-WRelLow), 0.0), 1.0)
WRel2B = min(max((WRel2-WRelLow)/(1.0-WRelLow), 0.0), 1.0)
! * Find transpiration from layers 1,2, soil evaporation, total evaporation
FWTraE = FracV * FWPT                               ! E-lim transpn         [m/day]
FWTr1W = FracV * 0.03 * WRel1B*ZSoil1*WVolSat1   ! W-lim transpn (1)     [m/day]
FWTr2W = FracV * 0.03 * WRel2B*ZSoil2*WVolSat2   ! W-lim transpn (2)     [m/day]
FWTraW = FWTr1W + FWTr2W                            ! total W-lim transpn   [m/day]
FWTraA = min(FWTraE, FWTraW) * FWTr1W/FWTraW        ! total transpn (1)     [m/day]
FWTraB = min(FWTraE, FWTraW) * FWTr2W/FWTraW        ! total transpn (2)     [m/day]
FWSoil = FWPT * (1.0-FracV) * WRel1B**(PwrFWSoil+1) ! soil evap             [m/day]
FWE    = FWTraA + FWTraB + FWSoil                   ! total evap            [m/day]
FWTra  = FWTraA + FWTraB                            ! total transpiration   [m/day]
! * Find latent and sensible heat fluxes
PhiE   = FWE * (RhoW*Rlat)/(DayltFrac*SecDay)       ! daylt latent flux     [W/m2]
PhiH   = PPc * (PhiAi - PhiE)                       ! daylt sens flux       [W/m2]
! * Find surface temperature and deficit
TempS  = TempA + PhiH / (RhoA*Capp*Gaero)               ! surface temp      [deg C]
DefS   = DefA + (EpsA*PhiH - PhiE) / (RhoA*Rlat*Gaero)  ! surface deficit   [molW/molA]

! --------------------------
! LEACHING AND RUNOFF FLUXES
! --------------------------
! * Find PRIOR water leaching fluxes out of layers 1,2
FWLchA   = HyConSat1 * WRel1B**(PwrFWLch+1)         ! leaching layer 1->2   [m/day]
FWLchB   = HyConSat2 * WRel2B**(PwrFWLch+1)         ! leaching ex layer 2   [m/day]
! * Prevent leaching from overdrawing water stores:
!   leaching flux can't be more than 0.5 of largest possible extraction in DelT
FWLchA   = min(FWLchA, 0.5*WRel1B*ZSoil1*WVolSat1/1.0 + FWPrec - FWTraA - FWSoil)
FWLchB   = min(FWLchB, 0.5*WRel2B*ZSoil2*WVolSat2/1.0 + FWLchA - FWTraB)
! * Apply saturation test in layer 2 to set FWLchA
WRel2XS  = WRel2 + 1.0 * (FWLchA - FWLchB - FWTraB) / (ZSoil2*WVolSat2)
where (WRel2XS > 1.0)                               ! layer 2 excess: reduce FWLchA
  FWLchA = FWLchA - (WRel2XS - 1.0) * (ZSoil2*WVolSat2)/1.0
end where
! * Apply saturation test in layer 1 to set FWRun
WRel1XS  = WRel1 + 1.0 * (FWPrec - FWLchA - FWTraA - FWSoil) / (ZSoil1*WVolSat1)
where (WRel1XS > 1.0)                               ! layer 1 excess: runoff occurs
  FWRun  = (WRel1XS - 1.0) * (ZSoil1*WVolSat1)/1.0
elsewhere                                           ! no layer 1 excess: no runoff
  FWRun  = 0.0
end where
FWDis    = FWRun + FWLchB                           ! total point discharge [m/day]

! -----------------------------------------------------------
! PLANT GROWTH = NET PRIMARY PRODUCTION = FCGro [molC/m2/day]  
! -----------------------------------------------------------
FQ     = SolarMJ * QMJSolar                         ! incident PAR flux [molQ/m2/d]
FCGroL = alfaQ * FracV * FQ                         ! Light-lim NPP     [molC/m2/day]
CO2C   = 37.0e-6 * (1.37**((TempS-25.0)/10.0))      ! CO2 comp pt       [molC/molA]
                                                    ! (von Caemmerer 2000)
alfaWdef = alfaWmul * (CO2A-CO2C) / (1.6 * max(DefS,0.001)) 
                                                    ! WUE from def      [molC/molW]
if (alfaWpri > 0.0) then                            ! use one of alfaWpri, alfaWdef
  alfaW = alfaWpri
else
  alfaW = alfaWdef
end if
FCGroW = alfaW * rhoW * (FWTraA + FWTraB)           ! Water-lim NPP     [molC/m2/day]
FCGro  = (FCGroL * FCGroW) / (FCGroL + FCGroW)      ! NPP = FCGro       [molC/m2/day]
                                                    ! Leaf allocation coefficient
where ((WRel1 .ge. 0.0) .and. (WRelA0 .ge. 0.0) .and. (FWTraW .gt. 0.0))
  WRelAv   = (WRel1*FWTr1W + WRel2*FWTr2W) / FWTraW
  AllocLea = sqrt(WRelAv) / (sqrt(WRelAv) + sqrt(WRelA0))
elsewhere
  AllocLea = 0.5
endwhere 

! --------------
! RATE EQUATIONS
! --------------
dwcoladt = (FWPrec - FWRun - FWLchA - FWTraA - FWSoil) / (ZSoil1*WVolSat1)
dwcolbdt = (FWLchA - FWLchB - FWTraB)                  / (ZSoil2*WVolSat2)
dCLeadt  = AllocLea*FCGro - RateCLea*CLea
ImBalA   = (ZSoil1*WVolSat1)*dwcoladt - (FWPrec - FWRun - FWLchA - FWTraA - FWSoil) 
ImBalB   = (ZSoil2*WVolSat2)*dwcolbdt - (FWLchA - FWLchB - FWTraB)

END SUBROUTINE WaterDynFlux

!********************##########################################################

SUBROUTINE WaterDynObs (TTime, XX, MM, UU, VV, AAP,AAPh, DD, FF,   ZZP,ZZPh,ZZC)
!-------------------------------------------------------------------------------
! * Observation models for NDVI, TempSt (at overpass time) and runoff
!-------------------------------------------------------------------------------
USE TypeDef
USE Constants
USE Utils
USE PointerModule
implicit none
! Global variables
real(sp),intent(in) :: TTime(:)     ! TTime(ntt) = time variables
real(sp),intent(in) :: XX(:,:)      ! XX(np,nxx) = state variables
real(sp),intent(in) :: MM(:,:)      ! MM(np,nmm) = met variables
real(sp),intent(in) :: UU(:)        ! UU(nuu)    = spatially uniform parameters
real(sp),intent(in) :: VV(:,:)      ! VV(np,nvv) = spatially variable parameters
real(sp),intent(inout) :: AAP(:,:)     ! AAP(np,naaP) = actual point obs (for ancillaries)
real(sp),intent(inout) :: AAPh(:,:,:)     ! AAP(np,naaP,ntime) = actual point obs (for ancillaries)
real(sp),intent(in)::  FF(:,:)      ! FF(np,nff) = fluxes and dXXdt
real(sp),intent(inout):: ZZP(:,:)     ! ZZP(np,nzzP) = predicted point obs
real(sp),intent(inout):: ZZPh(:,:,:)     ! ZZP(np,nzzP,ntime) = predicted point obs (hourly)
real(sp),intent(out):: ZZC(:,:)     ! ZZC(nc,nzzC) = predicted catchment obs
real(sp),intent(inout) :: DD(:,:)        ! DD(ndd)    = derived quantities

! Local variables
real(sp):: YearFrac                 ! year fraction                     [-]
real(sp):: TDDawn(size(XX,1))       ! time as dayfrac                   [-]
real(sp):: TDDusk(size(XX,1))       ! time as dayfrac                   [-]
real(sp):: TDTmax(size(XX,1))       ! time as dayfrac                   [-]
real(sp):: ALSTTime2(size(XX,1))    ! complete ALSTTime                 [-]
!real(sp):: FF(size(XX,1),nff)       ! fluxes and dXXdt
!real(sp):: DD(size(XX,1),ndd)       ! derived quantities
real(sp):: PhiHt(size(XX,1))        ! sens heat flux at overpass time   [W/m2]
!real(sp):: TempAt(size(XX,1))       ! air temp at overpass time         [deg C] 
real(sp):: FWRunC, FWLchC           ! Total point runoff, leaching over catchment
integer(i4b) :: iCatchMap(size(CatchMap,1))
integer(i4b) :: ic
type(dmydate):: TTDate              ! current date
real(sp),allocatable,save:: ZDisCMAcc(:)    ! accumulator for ZDisCM
integer(i4b),        save:: InitAccFlag = 0 ! initialisation flag

!-------------------------------------------------------------------------------

! * Allocate and initialise ZDisCMAcc (must use allocatable array to permit SAVE attribute)
if (InitAccFlag == 0) then
  allocate(ZDisCMAcc(nc))
  ZDisCMAcc = 0.0
  InitAccFlag = 1
end if



! * Reference output arrays (values will be overwritten)
!ZZP = ErrVal
!ZZPh = ErrVal
ZZC = ErrVal

! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime, TargetXX=XX, TargetMM=MM,         &
               TargetUU=UU, TargetVV=VV, TargetFF=FF, TargetDD=DD,  & 
               TargetAAP=AAP,TargetAAPh=AAPh, TargetZZP=ZZP, TargetZZPh=ZZPh,TargetZZC=ZZC)

! * Find fluxes and derived quantities with current stores, met, params
!CALL WaterDynFlux (TTime, XX, MM, UU, VV,   FF, DD)

! * Find TTDate, YearFrac
TTDate   = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
YearFrac = mod(YearDay(TTDate)/365.25, 1.0)     ! mod handles 31 Dec in leap year

! --------------------------
! OBSERVATION MODEL FOR NDVI
! --------------------------
ZNDVI = ErrVal
if (FracVegFlag <= 0) then
  ZNDVI    = cN0 + cn1*FracVExt + cn2*(FracVExt**2)
  else
  ZNDVI    = cN0 + cn1*FracVCLea + cn2*(FracVCLea**2)
end if

! --------------------------------------------------------------------
! OBSERVATION MODEL FOR SURFACE TEMPERATURE AT SATELLITE OVERPASS TIME
! --------------------------------------------------------------------

    
!   * All times TD are dayfracs (0 to 1, midnight=0)
!   * Overpass time is TDPass = ALSTTime/24.0   (don't define a separate symbol)
!   * Noon is TDNoon = 0.5                      (don't define a separate symbol)
!where (ALSTTime >= 0.0 .and. ALSTTime <= 24.0)  ! compute ZLST for all points
!	ALSTTime2 = ALSTTime
!elsewhere
!	ALSTTime2 = 12.0                              ! by setting ALSTTime to default if NA
!end where

!if (ForwardModelFlag.eq.5) then  ! temp fix vh
    where (ALST > ErrVal+1)
	  ZLST = TrEff                                 ! ZLST calculated from hourly CABLE output
	  DelTsTaZ = ZLST - TempAt
      DelTsTaA = ALST - TempAt
	elsewhere
	  DelTsTaA = ErrVal
	  ZLST = ErrVal
	  DelTsTaZ = ErrVal
	end where

	where (AphiRneth< ErrVal+1)
	  ZphiRneth = ErrVal
	end where
	where (AphiHh< ErrVal+1)
	  ZphiHh = ErrVal
	 end where
	where (AphiEh< ErrVal+1)
	  ZphiEh = ErrVal
	 end where
	where (AphiNEEh< ErrVal+1)
	  ZphiNEEh = ErrVal
	 end where
	where (AphiGh< ErrVal+1)
	  ZphiGh = ErrVal
	 end where
	where (ATsoilh< ErrVal+1)
	  ZTsoilh = ErrVal
	 end where
!end if
!else	
!	TDDawn = 0.5*(1.0-DayltFrac)            ! time of dawn
!	TDDusk = 1.0 - TDDawn                   ! time of dusk
!	TDTmax = TDDawn + DayltFrac*TimeTxFrac  ! time of Tmax
!	!  TDPass = ALSTTime2/24.0              ! overpass time (no separate symbol)
!	!  TDNoon = 0.5                         ! noon (no separate symbol)
!	PhiHt  = max( 2.0 * PhiH * cos(Pi*(ALSTTime2/24.0 - 0.5) / (TDDusk - TDDawn)),  &
!				  0.0 )                     ! sens heat flux at overpass time
!	TempAt = 0.5 * (TempMax + TempMin) +  & ! air temp at overpass time
!			 0.5 * (TempMax - TempMin) *  &
!			 cos(Pi*(TDTmax - ALSTTime2/24.0) / (2.0*(TDTmax - TDDawn)))
!	ZLST     = TempAt + PhiHt / (RhoA * Capp * Gaero)
!endif
	!ZLSTTime = ALSTTime2


! ---------------------------------------------------------
! OBSERVATION MODEL FOR DAILY AND MONTHLY CATCHMENT OUTFLOW
! ---------------------------------------------------------

! * Observation model for daily catchment outflow ZDisCD(nc)
iCatchMap(:) = nint(CatchMap(:))        ! integer catchment map (length nc)
do ic=1,nc
  FWRunC = sum( FWRun,  mask = ( (iCatchMap==iCatchID(ic)) ) ) /    &
                 count( mask = ( (iCatchMap==iCatchID(ic)) ) )
                                        ! total runoff over catchment [m/day]
  FWLchC = sum( FWLchB, mask = ( (iCatchMap==iCatchID(ic)) ) ) /    &
                 count( mask = ( (iCatchMap==iCatchID(ic)) ) )
                                        ! total leaching over catchment [m/day]
  if (TZRunDef > 0) then                ! low-pass filtered runoff
    ZRunCD(ic) = (1.0 - 1.0/TZRunDef) * ZRunCD(ic) + (1.0/TZRunDef) * FWRunC
  else
    ZRunCD(ic) = (1.0 - 1.0/TZRunCat(ic)) * ZRunCD(ic) + (1.0/TZRunCat(ic)) * FWRunC
  end if
  if (TZLchDef > 0) then                ! low-pass filtered leaching
    ZLchCD(ic) = (1.0 - 1.0/TZLchDef) * ZLchCD(ic) + (1.0/TZLchDef) * FWLchC
  else
    ZLchCD(ic) = (1.0 - 1.0/TZLchCat(ic)) * ZLchCD(ic) + (1.0/TZLchCat(ic)) * FWLchC
  end if
  ZDisCD(ic) = ZRunCD(ic) + ZLchCD(ic)  ! daily catchment outflow [m/day]
end do

! * Observation model for monthly catchment outflow ZDisCM(nc)
if (nEnsemble ==1) then
!   * Model 1 (use with PEST): 
!     * ZDisCM = monthly accumulation of ZDisCD on last day of each month
!     * ZDisCM = ErrVal on all other days
ZDisCMAcc = ZDisCMAcc + ZDisCD          ! accumulate monthly point outflow
if (EndMonth(TTDate)) then              ! end of month, so:
  ZDisCM  = ZDisCMAcc                   ! * set ZDisCM to accumulated daily outflow
  ZDisCMAcc = 0.0                       ! * reinitialise accumulator
else
  ZDisCM  = ErrVal                      ! not end of month, return error
end if
else if (nEnsemble >= 2) then
!   * Model 2 (use with EnKF): ZDisCM is discharge in mm/mth, defined on all days
ZDisCM = ZDisCD * DaysInMonth(TTDate)   ! monthly catchment outflow [m/mth]
else
  STOP "WaterDynObs: illegal nEnsemble"
end if

END SUBROUTINE WaterDynObs

!###############################################################################

SUBROUTINE SpatialStats (QQ, QQStats)
!-------------------------------------------------------------------------------
! * Calculates spatial statistics of array QQ(np,nqq) over np dimension, 
!   omitting null values
!   * QQStats(1,iqq,1) = mean of QQ(:,iqq)
!   * QQStats(1,nqq,2) = SD of QQ(:,iqq)
!   * QQStats(1,nqq,3) = max value of  of QQ(:,iqq)
!   * QQStats(1,nqq,4) = min value of QQ(:,iqq)
!   * QQStats(1,nqq,5) = count of non-null values in QQ(:,iqq)
! * First dimension of QQStats (1) is a quasi-spatial dimension so that QQ(1,nqq,1)
!   has shape compatible with QQ(np,nqq), for later call to TimeAverage
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Global variables
real(sp),intent(in) :: QQ(:,:)          ! QQ(np,nqq) = array to be averaged
real(sp),intent(out):: QQStats(:,:,:)   ! QQStats(1,iqq,istat)
                                        ! = stats (over np) of QQ(np,iqq)
! Local variables
real(sp)    :: QQM1, QQM2, QQMax, QQMin, QQCount
integer(i4b):: npp, nqq, ip, iqq
!-------------------------------------------------------------------------------

npp = size(QQ,1)
nqq = size(QQ,2)
if (npp == 1) then                      ! single spatial point (: does iqq=1,nqq)
  QQStats(1,:,1) = QQ(1,:)
  QQStats(1,:,2) = 0.0
  QQStats(1,:,3) = QQ(1,:)
  QQStats(1,:,4) = QQ(1,:)
  QQStats(1,:,5) = 1.0
else                                    ! multiple spatial points
  do iqq=1,nqq
    QQCount = 0.0
    QQM1    = 0.0
    QQM2    = 0.0
    QQMax   = -1.0e20
    QQMin   = +1.0e20
    do ip=1,npp
      if (QQ(ip,iqq) < ErrVal+1.0) cycle 
      QQCount = QQCount + 1
      QQM1    = QQM1    + QQ(ip,iqq)
      QQM2    = QQM2    + QQ(ip,iqq)**2
      QQMax   = max(QQMax, QQ(ip,iqq))
      QQMin   = min(QQMin, QQ(ip,iqq))
    end do
    if (nint(QQCount) >= 1) then
      QQStats(1,iqq,1) = QQM1 / QQCount
      QQStats(1,iqq,2) = (QQM2 - QQM1**2) / QQCount
      QQStats(1,iqq,3) = QQMax
      QQStats(1,iqq,4) = QQMin
    else 
      QQStats(1,iqq,1) = ErrVal
      QQStats(1,iqq,2) = ErrVal
      QQStats(1,iqq,3) = ErrVal
      QQStats(1,iqq,4) = ErrVal
    end if
    QQStats(1,iqq,5) = QQCount
  end do
end if
END SUBROUTINE SpatialStats

!###############################################################################

SUBROUTINE TimeAverage (AvFlag, QQ, QQTAv, QQCount)
!-------------------------------------------------------------------------------
! * Calculates cumulative time average of array QQ(np,nqq), omitting null values
! * AvFlag = (0,1) for (initialise, find cumulative av)
! * If no counts during use QQValidFlag to set QQTAv=error in any cell for which
!   data contains a null value during time averaging period
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Global variables
integer(i4b),intent(in):: AvFlag        ! (0,1) = (init, find cumulative av)
real(sp),intent(in)    :: QQ(:,:)       ! QQ(np,nqq) = array to be averaged
real(sp),intent(inout) :: QQTAv(:,:)    ! cumulative time average
real(sp),intent(inout) :: QQCount(:,:)  ! counters
!-------------------------------------------------------------------------------

if (AvFlag == 0) then
! * initialise
  QQCount = 0.0
  QQTAv   = 0.0
  if (size(QQ,1) /= size(QQTAv,1))       STOP "TimeAverage: array size mismatch 1A"
  if (size(QQ,1) /= size(QQCount,1))     STOP "TimeAverage: array size mismatch 1B"
  if (size(QQ,2) /= size(QQTAv,2))       STOP "TimeAverage: array size mismatch 2A"
  if (size(QQ,2) /= size(QQCount,2))     STOP "TimeAverage: array size mismatch 2B"
else if (AvFlag == 1) then
! * accumulate
  where (QQCount < 0.5)                 ! no data yet, start from QQTAv = 0 
    QQTAv    = 0.0
  end where
  where (QQ > ErrVal+1.0)               ! conditional because some QQ are errors
    QQCount = QQCount + 1.0
    QQTAv   = ( (QQCount-1.0)*QQTAv + QQ ) / QQCount
  end where
  where (QQCount < 0.5)                 ! still no data, return error
    QQTAv    = ErrVal
  end where
else
  stop "TimeAverage: illegal AvFlag"
end if

END SUBROUTINE TimeAverage
!###############################################################################

SUBROUTINE TimeAverageh (AvFlag, QQ, QQTAv, QQCount)
!-------------------------------------------------------------------------------
! * Calculates cumulative time average of array QQ(np,nqq), omitting null values
! * AvFlag = (0,1) for (initialise, find cumulative av)
! * If no counts during use QQValidFlag to set QQTAv=error in any cell for which
!   data contains a null value during time averaging period
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Global variables
integer(i4b),intent(in):: AvFlag        ! (0,1) = (init, find cumulative av)
real(sp),intent(in)    :: QQ(:,:,:)       ! QQ(np,nqq,ntime) = array to be averaged
real(sp),intent(inout) :: QQTAv(:,:,:)    ! cumulative time average
real(sp),intent(inout) :: QQCount(:,:,:)  ! counters
!-------------------------------------------------------------------------------

if (AvFlag == 0) then
! * initialise
  QQCount = 0.0
  QQTAv   = 0.0
  if (size(QQ,1) /= size(QQTAv,1))       STOP "TimeAverage: array size mismatch 1A"
  if (size(QQ,1) /= size(QQCount,1))     STOP "TimeAverage: array size mismatch 1B"
  if (size(QQ,2) /= size(QQTAv,2))       STOP "TimeAverage: array size mismatch 2A"
  if (size(QQ,2) /= size(QQCount,2))     STOP "TimeAverage: array size mismatch 2B"
else if (AvFlag == 1) then
! * accumulate
  where (QQCount < 0.5)                 ! no data yet, start from QQTAv = 0 
    QQTAv    = 0.0
  end where
  where (QQ > ErrVal+1.0)               ! conditional because some QQ are errors
    QQCount = QQCount + 1.0
    QQTAv   = ( (QQCount-1.0)*QQTAv + QQ ) / QQCount
  end where
  where (QQCount < 0.5)                 ! still no data, return error
    QQTAv    = ErrVal
  end where
else
  stop "TimeAverage: illegal AvFlag"
end if

END SUBROUTINE TimeAverageh


!###############################################################################

SUBROUTINE TimeSeriesOutput (iTSOPFlag, AvType, ipDiag, UU, VV, XXi,    &
                             TTime, MM,RR, XX, FF, DD, ZZP, AAP)
!-------------------------------------------------------------------------------
! Time series output of MM, XX, FF, DD, ZZP, AAP
! Output is for one spatial location only (may be a spatial average)
! iTSOPFlag=0: write headers
! iTSOPFlag=1: time series output at one time point
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
implicit none
! Global variables
integer(i4b),intent(in) :: iTSOPFlag    ! (0,1) = (init, output)
integer(i4b),intent(in) :: AvType       ! type of time average
integer(i4b),intent(in) :: ipDiag       ! point for diagnostic TS output (0 -> DomAv)
real(sp),intent(in) :: UU(:)            ! UU(nuu) = spatially uniform parameters
real(sp),intent(in) :: VV(:)            ! VV(nvv) = spatially variable parameters
real(sp),intent(in) :: XXi(:)           ! XXi(nxx) = state variables: initial values
real(sp),intent(in) :: TTime(:)         ! TTime(ntt) = time variables
real(sp),intent(in) :: MM(:)            ! MM(nmm) = met variables
real(sp),intent(in) :: RR(:)            ! RR(nrr) = remote-sensing driver variables
real(sp),intent(in) :: XX(:)            ! XX(nxx) = state variables
real(sp),intent(in) :: FF(:)            ! FF(nff) = fluxes and dXXdt
real(sp),intent(in) :: DD(:)            ! DD(ndd) = derived quantities
real(sp),intent(in) :: ZZP(:)           ! ZZP(nzzP) = predicted point observables
real(sp),intent(in) :: AAP(:)           ! AAP(naaP) = actual point observations
! Local variables
integer(i4b):: i, j, iunit              ! counters
character(200):: OutFile                ! output filename (temporary)
!-------------------------------------------------------------------------------

! * Test input flags for validity
if (AvType < 0 .or. AvType > nAvType) STOP "TimeSeriesOutput: illegal AvType"
if (ipDiag < 0 .or. ipDiag > np)      STOP "TimeSeriesOutput: illegal ipDiag"
! * set unit for output file
iunit = 150 + AvType

! * INITIALISATION
!   * open output file and write headers
if (iTSOPFlag == 0) then
  write(OutFile,"('DomAv',i1,'.csv')") AvType
  GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // trim(OutFile)
  OPEN (unit=iunit, file=trim(GenFileName))
  if (ipDiag == 0) then
    write(iunit,                                                        &
      "('WaterDyn: domain-averaged time series with AvType = ', i6)") AvType
  else
    write(iunit,                                                        &
      "('WaterDyn: daily point diagnostic TS at ip = ', i6)") ipDiag
  end if
  write(iunit,"('Domain:,,',  a)")     trim(DomName)
  write(iunit,"('Run ID:,,',  a)")     trim(RunID)
  write(iunit,"('Comment:,,', a)")     trim(Comment)
  write(iunit,"('UU: spatially uniform parameters')")
  write(iunit,"(100(a,    ','))") (trim(UUName(i,1)), i=1,nuu) 
  write(iunit,"(100(e13.6,','))") (UU(i),             i=1,nuu) 
  write(iunit,"('VV: spatially explicit parameters')")
  write(iunit,"(100(a,    ','))") (trim(VVName(i,1)), i=1,nvv) 
  write(iunit,"(100(e13.6,','))") (VV(i),             i=1,nvv) 
  write(iunit,"('XX: initial stores')")
  write(iunit,"(100(a,    ','))") (trim(XXName(i,1)), i=1,nxx) 
  write(iunit,"(100(e13.6,','))") (XXi(i),            i=1,nxx) 
  if (FracVegFlag <= 0) then
    write(iunit,"('FracVegFlag=0 => FracV from external VegFPC')")
  else
    write(iunit,"('FracVegFlag=0 => FracV from internal CLea')")
  end if
  do i=1,2                              ! skip lines
    write(iunit,"()")
  end do
  write(iunit,102) ('TT(',i,')',i=1,ntt-1),    &
                   ('MM(',i,')',i=1,nmm),      &
				   ('RR(',i,')',i=1,nrr),      &
                   ('XX(',i,')',i=1,nxx),      &
                   ('FF(',i,')',i=1,nff),      &
                   ('DD(',i,')',i=1,ndd),      &
                   ('ZZP(',i,')','AAP(',i,')',i=1,nzzP)
  102 format(200(a,i2,a,','))
  do j=2,1,-1
    write(iunit,101) (trim(TTName(i,j)), i=1,ntt-1),    &
                     (trim(MMName(i,j)), i=1,nmm),      &
					 (trim(RRName(i,j)), i=1,nrr),      &
                     (trim(XXName(i,j)), i=1,nxx),      &
                     (trim(FFName(i,j)), i=1,nff),      &
                     (trim(DDName(i,j)), i=1,ndd),      &
                     (trim(ZZPName(i,j)),trim(AAPName(i,j)), i=1,nzzP)
    101 format(200(a,','))
  end do

! * OUTPUT
else if (iTSOPFlag == 1) then
  write(iunit,100) TTime(1:ntt-1), MM,RR, XX, FF, DD, (ZZP(i),AAP(i),i=1,nzzP)
  100 format(200(e13.5,','))
else
  stop "TimeSeriesOutput: illegal iTSOPFlag"
end if

END SUBROUTINE TimeSeriesOutput

!###############################################################################
SUBROUTINE TimeSeriesOutput_bin (iTSOPFlag,TSOPFlag, AvType,AvName, TTime,QQ,QQName,QQid)
!-------------------------------------------------------------------------------
! Time series output of MM, XX, FF, DD, ZZP, AAP
! iTSOPFlag=0: write headers
! iTSOPFlag=1: time series output at one time point
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
USE PointerModule
implicit none
! Global variables
integer(i4b),intent(in) :: iTSOPFlag    ! (0,1) = (init, output)
integer(i4b), intent(in) :: TSOPFlag(:)    ! flag for map output
integer(i4b),intent(in) :: AvType       ! type of time average
integer(i4b),intent(in) :: QQid       ! integer for unse in iunit
real(sp),intent(in) :: TTime(:)         ! TTime(ntt) = time variables
real(sp),     intent(in) :: QQ(:,:)         ! QQ(np,nqq)  = array to be written
character(16),intent(in) :: QQname(:)       ! QQname(nqq) = name of QQ(:,iqq)
character(3), intent(in) :: AvName          ! (day,mth,ann,run)
character(8) :: StartDateChar         ! start date as character string (YYYYMMDD)
character(8) :: EndDateChar         ! end date as character string (YYYYMMDD)

! Local variables
integer(i4b):: i, j, iunit, iqq             ! counters

!-------------------------------------------------------------------------------


! * Construct date from TTime as char string
!CALL PointAll (TargetTTime=TTime)
!write(TTDateChar,"(i4.4,2i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay)

write(StartDateChar,"(i4.4,2i2.2)") (StartDate%year), &
                                     (StartDate%month), (StartDate%day)

write(EndDateChar,"(i4.4,2i2.2)") (EndDate%year), &
                                     (EndDate%month), (EndDate%day)

! * Test input flags for validity
if (AvType < 0 .or. AvType > nAvType) STOP "TimeSeriesOutput_bin: illegal AvType"


j=0

if (iTSOPFlag==0) then
		do iqq=1,size(QQ,2)
			if (TSOPFlag(iqq) < (4-AvType)) cycle
			j=j+1
			iunit = 350*(Avtype+1+QQID*100) + j
			!write(7,*) trim(QQname(iqq)),iTSOPFlag,Avtype, iunit
			GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
					trim(RunID) // '/OutputTS/'//StartDateChar//EndDateChar &
					// '_' //trim(AvName) // '_' // trim(QQname(iqq)) //        &
					'_TS.bin'
			OPEN(iunit, file=GenFileName, form='binary')
		enddo
else if (iTSOPFlag == 1) then
		do iqq=1,size(QQ,2)
			if (TSOPFlag(iqq)< (4-AvType)) cycle
			j=j+1
			iunit = 350*(Avtype+1+QQID*100) + j
			!write(7,*) trim(QQname(iqq)),iTSOPFlag,Avtype, iunit
			write(iunit) nint(TTDay), nint(TTMonth), nint(TTYear),QQ(:,iqq)
		enddo
else
		stop "TimeSeriesOutput: illegal iTSOPFlag"
end if

END SUBROUTINE TimeSeriesOutput_bin

!###############################################################################

SUBROUTINE CatchmentTSOutput (iTSOPFlag, AvType, UU, VV, XXi, iCatchID,     &
                              izzaaL, izzaaH, TTime, ZZC, AAC)
!-------------------------------------------------------------------------------
! Time series output of ZZC, AAC for all catchments, for selected observables
! (set by (izzaaL <= izz <= izzaaH) and (izzaaL <= iaa <= izzaaH))
! iTSOPFlag=0: write headers
! iTSOPFlag=1: time series output at one time point
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
implicit none
! Global variables
integer(i4b),intent(in):: iTSOPFlag     ! (0,1) = (init, output)
integer(i4b),intent(in):: AvType        ! type of time average
integer(i4b),intent(in):: iCatchID(:)   ! integer unique catchment ID numbers
integer(i4b),intent(in):: izzaaL,izzaaH ! selectors for observables to be output
real(sp),intent(in) :: UU(:)            ! UU(nuu) = spatially uniform parameters
real(sp),intent(in) :: VV(:)            ! VV(nvv) = spatially variable parameters
real(sp),intent(in) :: XXi(:)           ! XXi(nxx)= state variables: initial values
real(sp),intent(in) :: TTime(:)         ! TTime(ntt) = time variables
real(sp),intent(in) :: ZZC(:,:)         ! ZZC(nc,nzzC) = predicted catchment obs
real(sp),intent(in) :: AAC(:,:)         ! AAC(nc,naaC) = actual catchment obs
! Local variables
real(sp)      :: AACwrite(size(AAC,1),size(AAC,2))  ! array for writing output
integer(i4b)  :: i, j, iunit, ic        ! counters
character(200):: OutFile                ! output filename (temporary)
!-------------------------------------------------------------------------------

! * Test input flags for validity
if (AvType < 0 .or. AvType > nAvType) stop "CatchmentTSOutput: illegal AvType"
! * set unit for output file
iunit = 160 + AvType

! * INITIALISATION
!   * open output file and write headers
if (iTSOPFlag == 0) then
  write(OutFile,"('CatchAv',i1,'.csv')") AvType
  GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // trim(OutFile)
  OPEN (unit=iunit, file=trim(GenFileName))
  write(iunit,"('WaterDyn: catchment time series with AvType = ', i6)") AvType
  write(iunit,"('Domain:,,',  a)")     trim(DomName)
  write(iunit,"('Run ID:,,',  a)")     trim(RunID)
  write(iunit,"('Comment:,,', a)")     trim(Comment)
  write(iunit,"('UU: spatially uniform parameters')")
  write(iunit,"(100(a,    ','))") (trim(UUName(i,1)), i=1,nuu) 
  write(iunit,"(100(e13.6,','))") (UU(i),             i=1,nuu) 
  write(iunit,"('VV: spatially explicit parameters (spatial averages)')")
  write(iunit,"(100(a,    ','))") (trim(VVName(i,1)), i=1,nvv) 
  write(iunit,"(100(e13.6,','))") (VV(i),             i=1,nvv) 
  write(iunit,"('XX: initial stores (spatial averages)')")
  write(iunit,"(100(a,    ','))") (trim(XXName(i,1)), i=1,nxx) 
  write(iunit,"(100(e13.6,','))") (XXi(i),            i=1,nxx) 
  do i=1,3                              ! skip lines
    write(iunit,"()")
  end do
  write(iunit,102) ('TT(',i,')',i=1,ntt-1),             &
                   (('ZZC(',i,')','AAC(',i,')',i=izzaaL,izzaaH),ic=1,nc)
  102 format(400(a,i2,a,','))
  write(iunit,103) ((trim(ZZCName(i,2)),trim(AACName(i,2)),i=izzaaL,izzaaH),ic=1,nc)
  103 format(',,,,',400(a,','))
  write(iunit,101) ((iCatchID(ic),trim(ZZCName(i,1)),   &
                     iCatchID(ic),trim(AACName(i,1)),i=izzaaL,izzaaH),ic=1,nc)
  101 format(',,,,',400(i6,':',a,','))

! * OUTPUT
else if (iTSOPFlag == 1) then
  do ic=1,nc                            ! set AAC=error in catchments with no data
    if (DisCMDataExist(ic)) then
      AACwrite(ic,:) = AAC(ic,:)
    else
      AACwrite(ic,:) = ErrVal
    end if
  end do
  write(iunit,100) TTime(1:ntt-1), ((ZZC(ic,i), AACwrite(ic,i),i=izzaaL,izzaaH),ic=1,nc)
  100 format(400(e13.5,','))
else
  stop "TimeSeriesOutput: illegal iTSOPFlag"
end if

END SUBROUTINE CatchmentTSOutput

!###############################################################################

SUBROUTINE CatchmentRunAvgOutput (UU, VV, XXi, iCatchID, izzaaL, izzaaH, ZZC, AAC)
!-------------------------------------------------------------------------------
! Run-average output of ZZC, AAC for all catchments, for selected observables
! (set by (izzaaL <= izz <= izzaaH) and (izzaaL <= iaa <= izzaaH))
! This routine is called only once at end of run, with AvType=3 (whole-run av).
! Routine writes output with catchments as rows and (iCatchID, ZZ, AA) as columns,
! to fit data into xls spreadsheet.
!-------------------------------------------------------------------------------
USE TypeDef
USE Utils
implicit none
! Global variables
integer(i4b),intent(in):: iCatchID(:)   ! integer unique catchment ID numbers
integer(i4b),intent(in):: izzaaL,izzaaH ! selectors for observables to be output
real(sp),intent(in) :: UU(:)            ! UU(nuu) = spatially uniform parameters
real(sp),intent(in) :: VV(:)            ! VV(nvv) = spatially variable parameters
real(sp),intent(in) :: XXi(:)           ! XXi(nxx)= state variables: initial values
real(sp),intent(in) :: ZZC(:,:)         ! ZZC(nc,nzzC) = predicted catchment obs
real(sp),intent(in) :: AAC(:,:)         ! AAC(nc,naaC) = actual catchment obs
! Local variables
integer(i4b)  :: i, j, iunit, ic        ! counters
character(200):: OutFile                ! output filename (temporary)
!-------------------------------------------------------------------------------

! * set unit for output file
iunit = 170

! * INITIALISATION
!   * open output file and write headers
  write(OutFile,"('CatchAv3.csv')")
  GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //      &
                trim(RunID) // '/OutputTS/' // trim(OutFile)
  OPEN (unit=iunit, file=trim(GenFileName))
  write(iunit, "('WaterDyn: catchment whole-run averages with AvType = 3')")
  write(iunit,"('Domain:,,',  a)")     trim(DomName)
  write(iunit,"('Run ID:,,',  a)")     trim(RunID)
  write(iunit,"('Comment:,,', a)")     trim(Comment)
  write(iunit,"('UU: spatially uniform parameters')")
  write(iunit,"(50(a,    ','))") (trim(UUName(i,1)), i=1,nuu) 
  write(iunit,"(50(e13.6,','))") (UU(i),             i=1,nuu) 
  write(iunit,"('VV: spatially explicit parameters (spatial averages)')")
  write(iunit,"(50(a,    ','))") (trim(VVName(i,1)), i=1,nvv) 
  write(iunit,"(50(e13.6,','))") (VV(i),             i=1,nvv) 
  write(iunit,"('XX: initial stores (spatial averages)')")
  write(iunit,"(50(a,    ','))") (trim(XXName(i,1)), i=1,nxx) 
  write(iunit,"(50(e13.6,','))") (XXi(i),            i=1,nxx) 
  do i=1,4                              ! skip lines
    write(iunit,"()")
  end do
  write(iunit,101) (trim(ZZCName(i,2)),trim(AACName(i,2)), i=izzaaL,izzaaH)
  write(iunit,101) (trim(ZZCName(i,1)),trim(AACName(i,1)), i=izzaaL,izzaaH)
  101 format('iCatchID,',400(a,','))

! * OUTPUT
  do ic=1,nc
    write(iunit,100) iCatchID(ic),(ZZC(ic,i),AAC(ic,i), i=izzaaL,izzaaH)
    100 format(i10,',',400(e13.5,','))
  end do

! * close file
  CLOSE(iunit)

END SUBROUTINE CatchmentRunAvgOutput

!###############################################################################

SUBROUTINE MapOutput (QQ, QQname, MapOPFlag, AvName, TTime)
!-------------------------------------------------------------------------------
! * Map output of array QQ(np,nqq), as spatial np-vectors to nqq binary files
! * Output map only if MapOPFlag(iqq) = 1
! * Filename for spatial np-vector iqq is: AvName_QQname(iqq)_TTDateChar.bin
!-------------------------------------------------------------------------------
USE TypeDef
USE PointerModule
implicit none
! Global variables
real(sp),     intent(in) :: QQ(:,:)         ! QQ(np,nqq)  = array to be written
character(16),intent(in) :: QQname(:)       ! QQname(nqq) = name of QQ(:,iqq)
integer(i4b), intent(in) :: MapOPFlag(:)    ! flag for map output
character(3), intent(in) :: AvName          ! (day,mth,ann,run)
real(sp),      intent(in):: TTime(:)        ! TTime(ntt) = time variables
! Local variables
integer(i4b)  :: iqq
!-------------------------------------------------------------------------------

! * Construct date from TTime as char string
CALL PointAll (TargetTTime=TTime)
write(TTDateChar,"(i4.4,2i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay)

! * For each iqq: build filename, open file, write vector map, close file
do iqq=1,size(QQ,2)
  if (MapOPFlag(iqq) /= 1) cycle
  GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' // trim(RunID) //   &
                '/OutputMap/' // trim(AvName) // '_' // trim(QQname(iqq)) //        &
                '_' // trim(TTDateChar) // '.bin'
  OPEN(unit=170, file=GenFileName, form='binary')
  write(170) QQ(:,iqq)
  CLOSE (170)
end do

END SUBROUTINE MapOutput
!###############################################################################

SUBROUTINE MapOutput_tmp (QQ, QQname, MapOPFlag, AvName, TTime)
!-------------------------------------------------------------------------------
! * Map output of array QQ(np,nqq), as spatial np-vectors to nqq binary files
! * Output map only if MapOPFlag(iqq) = 1
! * temporary file, overwritten at each call
! * Filename for spatial np-vector iqq is: AvName_QQname(iqq)_tmp.bin
!-------------------------------------------------------------------------------
USE TypeDef
USE PointerModule
implicit none
! Global variables
real(sp),     intent(in) :: QQ(:,:)         ! QQ(np,nqq)  = array to be written
character(16),intent(in) :: QQname(:)       ! QQname(nqq) = name of QQ(:,iqq)
integer(i4b), intent(in) :: MapOPFlag(:)    ! flag for map output
character(3), intent(in) :: AvName          ! (day,mth,ann,run)
real(sp),      intent(in):: TTime(:)        ! TTime(ntt) = time variables
! Local variables
integer(i4b)  :: iqq
!-------------------------------------------------------------------------------

! * Construct date from TTime as char string
CALL PointAll (TargetTTime=TTime)
write(TTDateChar,"(i4.4,2i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay)

! * For each iqq: build filename, open file, write vector map, close file
do iqq=1,size(QQ,2)
  !if (MapOPFlag(iqq) /= 1) cycle
  GenFileName = trim(DataOutRootDir) // trim(DomName) // '/InitStores/' //   &
                 trim(AvName) // '_' // trim(QQname(iqq)) //        &
                '_' // 'tmp.bin'
  OPEN(unit=170, file=GenFileName, form='binary')
  write(170) QQ(:,iqq)
  CLOSE (170)
end do

END SUBROUTINE MapOutput_tmp
!###############################################################################

SUBROUTINE MapOutputh (QQ, QQname, MapOPFlag,Ihour_output, AvName, TTime)
!-------------------------------------------------------------------------------
! * Map output of array QQ(np,nqq), as spatial np-vectors to nqq binary files
! * Output map only if MapOPFlag(iqq) = 1
! * Filename for spatial np-vector iqq is: AvName_QQname(iqq)_TTDateChar.bin
!-------------------------------------------------------------------------------
USE TypeDef
USE PointerModule
implicit none
! Global variables
real(sp),     intent(in) :: QQ(:,:,:)         ! QQ(np,nqq,ntime)  = array to be written
character(16),intent(in) :: QQname(:)       ! QQname(nqq) = name of QQ(:,iqq)
integer(i4b), intent(in) :: MapOPFlag(:)    ! flag for map output
integer(i4b),  intent(in):: Ihour_output(24)    ! specifies which hours to output for AAPh and ZZPh
character(3), intent(in) :: AvName          ! (day,mth,ann,run)
real(sp),      intent(in):: TTime(:)        ! TTime(ntt) = time variables
! Local variables
integer(i4b)  :: iqq, itt
!-------------------------------------------------------------------------------

! extract hourly data to be written out


! * Construct date from TTime as char string
CALL PointAll (TargetTTime=TTime)
write(TTDateChar,"(i4.4,2i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay)

! * For each iqq: build filename, open file, write vector map, close file
do iqq=1,size(QQ,2)
  if (MapOPFlag(iqq) /= 1) cycle
  do itt = 1,size(QQ,3)
   if (Ihour_output(itt).eq.1) then
    if (ALL(QQ(:,iqq,itt)<(ErrVal+1))) cycle
	write(TTDateTimeChar,"(i4.4,3i2.2)") nint(TTYear), nint(TTMonth), nint(TTDay),itt
    GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' // trim(RunID) //   &
                '/OutputMap/' // trim(AvName) // '_' // trim(QQname(iqq)) //        &
                '_' // trim(TTDateTimeChar) // '.bin'
	OPEN(unit=170, file=GenFileName, form='binary')
	write(170) QQ(:,iqq,itt)
	CLOSE (170)
    endif
  enddo
end do

END SUBROUTINE MapOutputh

!###############################################################################

SUBROUTINE WritePESTOutput (Time0,TTDate,AAP,AAC,ZZP,ZZC,FFTAv,VV,JobFlag,AAPh,ZZPh,AAPTAv,ZZPTAv)
!-------------------------------------------------------------------------------
! ** CMT, March 2006 (revised around 19-feb-07)
! Write output for PEST (either to prepare for a PEST run, or during a PEST run)
! PESTOutputFlag=2: preparing for PEST, write actual obs to file 'ActualObs.dat' and 
!                   write corresponding instructions for ins file to 'Instructions.dat'
! PESTOutputFlag=3: preparing for PEST, write pseudo obs to file 'PseudoObs.dat' and 
!                   write corresponding instructions for ins file to 'Instructions.dat'
! PESTOutputFlag=4: during PEST run, write observables to file 'ModelObs.dat' 
!
! ** Updated 13 Mar 07 for allocatable arrays
! 
! ** Updated 28 Mar 07 to run entirely from one ctl file and a bat file.
!
! ** Revised 4 Apr 07 to include whole-of-run-mean obs choice
! Current valid choices of PESTOutputFlag are:
! PESTOutputFlag=1: whole-of-run obs, generates PEST files for an optimisation run
! PESTOutputFlag=2: actual obs, generates PEST files for an optimisation run
! PESTOutputFlag=4: single model run to be used during PEST run, write closen actual 
! observables to file 'ModelObs.dat' 
! PESTOutputFlag=5: single model run to be used during PEST run, write chosen 
! whole-of-run observables to file 'ModelObs.dat' 
! Point scale obs still disablesd

!-------------------------------------------------------------------------------
USE TypeDef
USE DateFunctions
USE PointerModule
implicit none

real(sp),intent(in) :: Time0           ! time at start of step
type(dmydate),intent(in):: TTDate      ! current date
real(sp),intent(in) :: ZZP(:,:)        ! ZZP(np,nzzP) = predicted point observables
real(sp),intent(in), optional :: ZZPh(:,:,:)        ! ZZP(np,nzzP,ntime) = predicted point observables (hourly)
real(sp),intent(in) :: ZZC(:,:)        ! ZZC(nc,nzzC) = predicted catchment observables
real(sp),intent(in) :: AAP(:,:)        ! AAP(np,naaP) = actual point observations
real(sp),intent(in) :: FFTAv(:,:,:)
real(sp),intent(in), optional :: AAPh(:,:,:)        ! AAPh(np,nAAPh,ntime) = actual point observables (hourly)
real(sp),intent(in), optional :: AAPTAv(:,:,:) ! (np,naaP,nAvType)
real(sp),intent(in), optional :: ZZPTAv(:,:,:) ! (np,nzzP,nAvType)
real(sp),intent(in) :: AAC(:,:)        ! AAC(nc,naaC) actual catchment observations
integer(i4b),intent(in) :: JobFlag     ! (0,1,2) = (open files, write output, close files; calc weights; overwrite soil moisture as anomalies)
real(sp),intent(in):: VV(:,:) 

real(sp) :: obs, err, wt, pred
integer(i4b) :: ia,ip,ic,it,k,j
integer(i4b) :: day, month, yr
character*34 str
character*2 itstr
character*120 :: str120
character,save :: charAZ
integer,save :: countobs
character :: charC = 'C'
character :: charP = 'P'
character :: charH = 'H'
character :: chhash = '#'
logical :: finished
integer(i4b) :: ios,idx,pcount
integer,save :: countobsh
real(sp),allocatable,save:: obsaveh(:,:,:)
real(sp),allocatable,save:: PredAveh(:,:,:)
integer(i4b),allocatable,save:: obscounth(:,:,:)
integer(i4b):: countobsgp
integer(i4b):: countobs2
integer(i4b)::	n_ip_ia(np,nAAP)
integer(i4b)::	n_ic_ia(nc,nAAC)
integer(i4b)::	n_ip_iah(np,nAAPh)
integer(i4b)::	nt(1000000) 
character*8:: obsgp_temp='00000000'
character*8:: obsgp(10000)='00000000'
character*80000 :: obsgp_str
integer(i4b)::obsgp_count(10000)=0
character*30:: dum
character*1:: dum1
real(sp):: RMSE(10000)=0.0, weight(10000) = 0.0, phi(10000) = 0.0
real(sp):: resid(1000000)=0.0
integer(i4b):: Iobsgp(1000000)
integer(i4b)::	countobsgpAAPh = 0
integer(i4b)::	countobsgpAAP(nAAP)
integer(i4b)::	countobsgpAAC = 0
integer(i4b):: P_Index
real(sp):: obsval, obssig
character(30):: obsstr, obsname
integer(i4b),allocatable:: A9count(:), A8count(:), A10count(:), A11count(:), A12count(:)
real(sp),allocatable:: A8mean(:), A9mean(:), A10mean(:), A11mean(:), A12mean(:)
integer(i4b),allocatable:: Z9count(:), Z8count(:), Z10count(:), Z11count(:), Z12count(:)
real(sp),allocatable:: Z8mean(:), Z9mean(:), Z10mean(:), Z11mean(:), Z12mean(:)
real(sp):: FWPrecC        ! Total precip over catchment
integer(i4b) :: iCatchMap(size(CatchMap,1))

if (.not.allocated(PredAveh)) allocate(PredAveh(np,naaph,ntime))
if (.not.allocated(ObsAveh)) allocate(ObsAveh(np,naaph,ntime))
if (.not.allocated(obscounth)) allocate(obscounth(np,naaph,ntime))
if (.not.allocated(A8count)) allocate(A8count(np))
if (.not.allocated(A9count)) allocate(A9count(np))
if (.not.allocated(A10count)) allocate(A10count(np))
if (.not.allocated(A11count)) allocate(A11count(np))
if (.not.allocated(A12count)) allocate(A12count(np))
if (.not.allocated(Z8count)) allocate(Z8count(np))
if (.not.allocated(Z9count)) allocate(Z9count(np))
if (.not.allocated(Z10count)) allocate(Z10count(np))
if (.not.allocated(Z11count)) allocate(Z11count(np))
if (.not.allocated(Z12count)) allocate(Z12count(np))
if (.not.allocated(A8mean)) allocate(A8mean(np))
if (.not.allocated(A9mean)) allocate(A9mean(np))
if (.not.allocated(A10mean)) allocate(A10mean(np))
if (.not.allocated(A11mean)) allocate(A11mean(np))
if (.not.allocated(A12mean)) allocate(A12mean(np))
if (.not.allocated(Z8mean)) allocate(Z8mean(np))
if (.not.allocated(Z9mean)) allocate(Z9mean(np))
if (.not.allocated(Z10mean)) allocate(Z10mean(np))
if (.not.allocated(Z11mean)) allocate(Z11mean(np))
if (.not.allocated(Z12mean)) allocate(Z12mean(np))
CALL PointAll (TargetVV=VV)

! Open files
if (JobFlag .eq. 0) then ! i.e. create and open files 
  if ((PESTOutputFlag .eq. 1) .or. (PESTOutputFlag .eq. 2)) then    
  ! Create PEST tpl file
    write(6,*) 'Create tpl file'
    open(unit=24,file='PESTCD25.tpl')
    write(24,'(a5)') 'ptf #'
    open(unit=11, file=trim(CtlFileName)//'.ctl', status='old')
    finished = .false.
    do while (.not.(finished))      
      read(11,'(a120)',IOSTAT=ios) str120
      if (ios .ne. 0) STOP "EOF in ctl file in WritePESTOutput"
      idx = index(str120,'PESTOutputFlag')
      if (idx .ne. 0) then
        if (PESTOutputFlag == 1) str120(1:1) = '5'
        if (PESTOutputFlag == 2) str120(1:1) = '4'
      end if  
      write(24,'(a)') trim(str120)
      if (str120(1:36) .eq. '! UU = spatially constant parameters') finished = .true.
    end do
    read(11,'(a120)',IOSTAT=ios) str120
    write(24,'(a)') trim(str120) 
    finished = .false.
    pcount = 0
    do while (.not.(finished))
      read(11,'(a120)',IOSTAT=ios) str120
      if (ios .ne. 0) STOP "EOF in ctl file in WritePESTOutput"  
      if (str120(1:1) .eq. '!') then 
         finished = .true.
		 write(24,'(a)') trim(str120)
      else 
        pcount = pcount + 1
        if (UUPESTFlag(pcount) .eq. 1) then
          idx = index(str120,',')
          write(24,'(a1,a16,a1,a100)') chhash,UUName(pcount,1),chhash,str120(idx:120) 
        else
          write(24,'(a)') trim(str120)
        end if  
      end if
    end do
    finished = .false.  
    do while (.not.(finished))
      read(11,'(a120)',IOSTAT=ios) str120
      if (ios .ne. 0) then
        finished = .true.
      else  
        write(24,'(a)') trim(str120)
      end if
    end do  
    close(11)
    close(24)
  ! Create end part of pst file (PESTpst2.txt)
    write(6,*) 'Create file PESTpst2.txt'
    open(unit=24,file='PESTpst2.txt')    
    write(24,'(a)') '* model command line'


!DEC$ IF DEFINED (__INTEL_COMPILER)
write(24,'(a)') '/home/cmar/hav014/bios_batch/PEST/runcable_hav014.sh'
!DEC$ ELSE
write(24,'(a)') 'c:/CableDyn/code/CableDyn25/release/CableDyn25 CD25ctlFile'
!DEC$ ENDIF 

    write(24,'(a)') '* model input/output'
    write(24,'(a)') 'PESTCD25.tpl CD25ctlFile.ctl'
    write(24,'(a)') 'PESTCD25.ins ModelObs_anomaly.dat'
    write(24,'(a)') '* prior information'
    close(24)
  ! Open files ActualObs.dat and ins file'  
    write(6,*) 'Writing actual obs to file ''ActualObs.dat''' 
    write(6,*) 'Creating ins file ''PESTCD25.ins'''   
    open(unit=21,file='ActualObs.dat',STATUS='REPLACE',ACTION='WRITE')
	!open(unit=221,file='ActualObs_temp.dat',STATUS='unknown',ACTION='WRITE')
    open(unit=22,file='PESTCD25.ins')
    write(22,'(a5)')'pif #'
    write(22,'(a8)')'#output#'
    countobs = 0
    charAZ = 'A'
  end if  ! ((PESTOutputFlag .eq. 1) .or. (PESTOutputFlag .eq. 2)) 
  if ((PESTOutputFlag == 2) .or.(PESTOutputFlag == 4) .or. (PESTOutputFlag == 5)) then
    open(unit=23,file='ModelObs.dat',STATUS='REPLACE',ACTION='WRITE')
    write(23,*) 'output'  
  end if 
end if ! (JobFlag .eq. 0) 

! Close files
if (JobFlag .eq. 2) then ! i.e. close files
  ! Write Whole-of-run-average observations or observables now  
  if ((PESTOutputFlag == 1) .or. (PESTOutputFlag == 5)) then  
   do ia = 1, naaC  
    if (AACflag(ia) .gt. 0) then 
     do ic = 1, nc 
        select case (PESTOutputFlag)
        case (1)     
          if (AACabsrel(ia) .eq. 'a') err = AACr(ia)
          if (AACabsrel(ia) .eq. 'r') err = AACr(ia)*AACTAv(ic,ia,3)
          write(21,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,AACTAv(ic,ia,3),1.0/err,AACName(ia,1)
		  write(221,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,AACTAv(ic,ia,3),1.0/err,AACName(ia,1)
          write(str,122) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia 
          write(22,*) str  ! writing instruction file
          countobs = countobs + 1
		  resid(countobs) = ZZCTAv(ic,ia,3)-AACTAv(ic,ia,3)
        case (5)     
          write(23,123) ZZCTAv(ic,ia,3)
        end select 
     end do
    end if 
   end do       
  end if
  ! Close files
  if ((PESTOutputFlag == 1) .or. (PESTOutputFlag == 2)) then
    write(6,*) 'Closing obs and instructions files'
    close(21)
	close(221)
    close(22)
  ! Create first part of pst file (PESTpst1.txt)  
    write(6,*) 'Create file PESTpst1.txt'
    open(unit=24,file='PESTpst1.txt')
    write(24,'(a)')'pcf'
    write(24,'(a)')'* control data'
    write(24,'(a)')'restart estimation'
    pcount = sum(UUPESTFlag)
	countobsgp = sum(AAPFlag) + sum(AACFlag) + sum(AAPhFlag)
    write(24,*) pcount,countobs,' 1 0', countobsgp
    write(24,'(a)') '1 1 single point 1 0 0'
    write(24,'(a)') '5.0 2.0 0.3 0.01 10'
    write(24,'(a)') '5 5 0.001'
    write(24,'(a)') '0.1'
    write(24,'(a)') '120 0.01 4 3 0.01 3'
    write(24,'(a)') '1 1 1'
    write(24,'(a)') '* parameter groups'
    write(24,'(a)') 'PP relative 0.01 0.0 always_2 2.0 parabolic'
    write(24,'(a)') '* parameter data'
    write(6,*) 'Calibrate parameters:-'
    do pcount = 1,nuu
      if (UUPESTFlag(pcount) .eq. 1) then
      write(6,*) UUName(pcount,1)
        write(24,'(a12,'' none factor '',f13.7,f13.7,f13.7,''  PP 1.0 0.0 1'')')    &
                                 UUName(pcount,1),UUinitCtl(pcount),                &
                                 UUminCtl(pcount),UUmaxCtl(pcount)
      end if
    end do
    write(24,'(a)') '* observation groups'
    if (AAPFlag(1) .eq. 1) write(24,*) 'ANDVI'
    if (AAPFlag(2) .eq. 1) write(24,*) 'ALST'
    if (AACFlag(1) .eq. 1) write(24,*) 'ADisCD'
    if (AACFlag(2) .eq. 1) write(24,*) 'ADisCM'
	if (AAPhFlag(6) .eq. 1) write(24,*) 'AphiRneth'
	if (AAPhFlag(7) .eq. 1) write(24,*) 'AphiHh'
	if (AAPhFlag(8) .eq. 1) write(24,*) 'AphiEh'
	if (AAPhFlag(9) .eq. 1) write(24,*) 'AphiNEEh'
	if (AAPFlag(8) .eq. 1) write(24,*) 'Asm0008'
	if (AAPFlag(9) .eq. 1) write(24,*) 'Asm0090'
    if (AAPFlag(10) .eq. 1) write(24,*) 'Asm0030'
	if (AAPFlag(11) .eq. 1) write(24,*) 'Asm3060'
	if (AAPFlag(12) .eq. 1) write(24,*) 'Asm6090'
	if (AAPFlag(2) .eq. 1) write(24,*) 'AphiRnet'
	if (AAPFlag(3) .eq. 1) write(24,*) 'AphiH'
	if (AAPFlag(4) .eq. 1) write(24,*) 'AphiE'
	if (AAPFlag(6) .eq. 1) write(24,*) 'AphiNPP'
	if (AAPFlag(7) .eq. 1) write(24,*) 'AphiGPP'
    write(24,'(a)') '* observation data'
    close(24)
	close(23)
  end if
  if (PESTOutputFlag == 4) then
    close(23)
  end if
end if

! overwrite PEST actualobs file with AA anomalies (soil moisture only)
if (JobFlag .eq. 2) then
	A8mean = 0.0
	A9mean = 0.0
	A10mean = 0.0
	A11mean = 0.0
	A12mean = 0.0
	Z8mean = 0.0
	Z9mean = 0.0
	Z10mean = 0.0
	Z11mean = 0.0
	Z12mean = 0.0
	A8count = 0
	A9count = 0
	A10count = 0
	A11count = 0
	A12count = 0
	Z8count = 0
	Z9count = 0
	Z10count = 0
	Z11count = 0
	Z12count = 0

	if (PESTOutputFlag .eq. 2) then
	
		write(6,*) 'Writing actual obs to file ''ActualObs_anomaly.dat'''    
		open(unit=21,file='ActualObs.dat')
		open(unit=221,file='ActualObs_anomaly.dat')
		do 
			read(21,*, end=91) dum
			if (INDEX(dum,'A08').gt.0.or.INDEX(dum,'A09').gt.0.or.INDEX(dum,'A10').gt.0.or.INDEX(dum,'A11').gt.0.or.INDEX(dum,'A12').gt.0) then 
				BACKSPACE(21)
				read(21,121) day,month,yr,charC,ip,charAZ,ia,obs,wt,AAPName(ia,1)
				if (ia.eq.08) then
					A8mean(ip) = A8mean(ip) + obs
					A8count(ip) = A8count(ip)+1
				endif
				if (ia.eq.09) then
					A9mean(ip) = A9mean(ip) + obs
					A9count(ip) = A9count(ip)+1
				endif
				if (ia.eq.10) then
					A10mean(ip) = A10mean(ip) + obs
					A10count(ip) = A10count(ip)+1
				endif
				if (ia.eq.11) then
					A11mean(ip) = A11mean(ip) + obs
					A11count(ip) = A11count(ip)+1
				endif
				if (ia.eq.12) then
					A12mean(ip) = A12mean(ip) + obs
					A12count(ip) = A12count(ip)+1
				endif
			endif
		enddo

		91 close(21)

		WHERE (A8count.gt.0)
			A8mean = A8mean/A8count
		ENDWHERE
		WHERE (A9count.gt.0)
			A9mean = A9mean/A9count
		ENDWHERE
		WHERE (A10count.gt.0)
			A10mean = A10mean/A10count
		ENDWHERE
		WHERE (A11count.gt.0)
			A11mean = A11mean/A11count
		ENDWHERE
		WHERE (A12count.gt.0)
			A12mean = A12mean/A12count
		ENDWHERE

		open(unit=21,file='ActualObs.dat')
		do 

			read(21,*, end=991) dum
			BACKSPACE(21)
			if (INDEX(dum,'AH').gt.0) then  ! hourly point observables
				read(21,1221) day,month,yr,charP,ip,charAZ,charH,ia,itstr,obs,wt,AAPhName(ia,1)
				write(221,1221) day,month,yr,charP,ip,charAZ,charH,ia,itstr,obs,wt,AAPhName(ia,1)
			elseif (INDEX(dum,'C').gt.0) then    ! catchment observables
				read(21,121) day,month,yr,charC,ic,charAZ,ia,obs,wt,AACName(ia,1)	
				write(221,121) day,month,yr,charC,ic,charAZ,ia,obs,wt,AACName(ia,1)
			else           ! point observables
				
				read(21,121) day,month,yr,charC,ip,charAZ,ia,obs,wt,AAPName(ia,1)
				if (ia==08) obs = obs - A8mean(ip)
				if (ia ==09) obs = obs - A9mean(ip)	
				if (ia ==10) obs = obs - A10mean(ip)	
				if (ia ==11) obs = obs - A11mean(ip)	
				if (ia ==12) obs = obs - A12mean(ip)		
				write(221,121) day,month,yr,charC,ip,charAZ,ia,obs,wt,AAPName(ia,1)
			endif
		enddo

		991 close(21)
		close(221)

	endif

	if (PESTOutputFlag == 4.or.PESTOutputFlag == 2) then
	! overwrite PEST modelobs file with ZZ anomalies (soil moisture only)
		write(6,*) 'Writing model obs to file ''ModelObs_anomaly.dat'''
		open(unit=23,file='ModelObs.dat')
 		open(unit=22,file='PESTCD25.ins',ACTION='READ')
		read(23,*) dum
		read(22,*) dum
		read(22,*) dum
		do
			read(23,*,end=9991) obsval
			read(22,*,end=9991) dum,dum
			if (INDEX(dum,'A08').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				Z8mean(ip) = Z8mean(ip) + obsval
				Z8count(ip) = Z8count(ip) + 1
			endif
			if (INDEX(dum,'A09').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				Z9mean(ip) = Z9mean(ip) + obsval
				Z9count(ip) = Z9count(ip) + 1
			endif
			if (INDEX(dum,'A10').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				Z10mean(ip) = Z10mean(ip) + obsval
				Z10count(ip) = Z10count(ip) + 1
			endif
			if (INDEX(dum,'A11').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				Z11mean(ip) = Z11mean(ip) + obsval
				Z11count(ip) = Z11count(ip) + 1
			endif
			if (INDEX(dum,'A12').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				Z12mean(ip) = Z12mean(ip) + obsval
				Z12count(ip) = Z12count(ip) + 1
			endif
		enddo
	9991	close(23)
		close(22)

		WHERE (Z8count.gt.0)
			Z8mean = Z8mean/Z8count
		ENDWHERE
		WHERE (Z9count.gt.0)
			Z9mean = Z9mean/Z9count
		ENDWHERE
		WHERE (Z10count.gt.0)
			Z10mean = Z10mean/Z10count
		ENDWHERE
		WHERE (Z11count.gt.0)
			Z11mean = Z11mean/Z11count
		ENDWHERE
		WHERE (Z12count.gt.0)
			Z12mean = Z12mean/Z12count
		ENDWHERE

		open(unit=23,file='ModelObs.dat')
		open(unit=223,file='ModelObs_anomaly.dat') 
		open(unit=224,file='ModelObs_anomaly_ASCII.dat') 
		open(unit=22,file='PESTCD25.ins',ACTION='READ')
		read(23,*) dum
		write(223,*) dum
		read(22,*) dum
		read(22,*) dum
		do
			read(23,*,end=9992) obsval
			read(22,*,end=9992) dum,dum
			if (INDEX(dum,'A08').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				obsval = obsval - Z8mean(ip)
			endif
			if (INDEX(dum,'A09').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				obsval = obsval - Z9mean(ip)
			endif
			if (INDEX(dum,'A10').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				obsval = obsval - Z10mean(ip)
			endif
			if (INDEX(dum,'A11').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				obsval = obsval - Z11mean(ip)
			endif
			if (INDEX(dum,'A12').gt.0) then
				P_Index = INDEX(dum,'P')
				read(dum(P_Index+1:P_Index+3),*) ip
				obsval = obsval - Z12mean(ip)
			endif
			write(223,"(e13.5)") obsval
			write(224,"(e13.5)") obsval
		enddo
	9992	close(23)
		close(22)
		close(223)
		close(224)

	endif
endif



! Overwrite ActualObs.dat with appropriate weights
if (JobFlag .eq. 2.and.PESTOutputFlag .eq. 2) then
	write(6,*) 'reading ActualObs_anomaly.dat'
	open(unit=221,file='ActualObs_anomaly.dat',STATUS='OLD',ACTION='READ')
	countobsgp=0
	countobs2=0
	countobsgpAAPh = 0
	countobsgpAAP = 0
	countobsgpAAC = 0
	n_ip_ia = 0
	n_ic_ia = 0
	n_ip_iah = 0
	nt = 0
	do 
		read(221,*, end=99) dum
		countobs2=countobs2+1
		BACKSPACE(221)
		if (INDEX(dum,'AH').gt.0) then  ! hourly point observables
			read(221,1221) day,month,yr,charP,ip,charAZ,charH,ia,itstr,obs,wt,AAPhName(ia,1)
			write(obsgp_temp,126) charP,charH,ip,ia
			n_ip_iah(ip,ia) = n_ip_iah(ip,ia)+1
			nt(countobs2) = n_ip_iah(ip,ia)+1 ! size of (ip,ia) set to which countobs_th obs belongs

			if (INDEX(obsgp_str,obsgp_temp).gt.0) then
				obsgp_count((INDEX(obsgp_str,obsgp_temp)-1)/8+1) = obsgp_count((INDEX(obsgp_str,obsgp_temp)-1)/8+1) + 1
				!RMSE((INDEX(obsgp_str,obsgp_temp)-1)/8+1) = RMSE((INDEX(obsgp_str,obsgp_temp)-1)/8+1) + resid(countobs2)**2
				Iobsgp(countobs2) = (INDEX(obsgp_str,obsgp_temp)-1)/8+1
			else
				countobsgp = countobsgp + 1
				countobsgpAAPh = countobsgpAAPh +1
				write(obsgp_str((countobsgp-1)*8+1:(countobsgp-1)*8+8),126) charP,charH,ip,ia
				obsgp_count(countobsgp) = 1
				!RMSE(countobsgp) = resid(countobs2)**2
				write(obsgp(countobsgp),126) charP,charH,ip,ia
				Iobsgp(countobs2) = countobsgp
			endif
		elseif (INDEX(dum,'C').gt.0) then    ! catchment observables
			read(221,121) day,month,yr,charC,ic,charAZ,ia,obs,wt,AACName(ia,1)	
			write(obsgp_temp,126) charC,'0',111,ia
			write(obsgp(countobs2),126) charC,'0',ic,ia
			n_ic_ia(ic,ia) = n_ic_ia(ic,ia)+1
			nt(countobs2) = n_ic_ia(ic,ia)+1 ! size of (ip,ia) set to which countobs_th obs belongs
			if (INDEX(obsgp_str,obsgp_temp).gt.0) then
				obsgp_count((INDEX(obsgp_str,obsgp_temp)-1)/8+1) = obsgp_count((INDEX(obsgp_str,obsgp_temp)-1)/8+1) + 1
				!RMSE((INDEX(obsgp_str,obsgp_temp)-1)/8+1) = RMSE((INDEX(obsgp_str,obsgp_temp)-1)/8+1) + resid(countobs2)**2
				Iobsgp(countobs2) = (INDEX(obsgp_str,obsgp_temp)-1)/8+1
			else
				countobsgp = countobsgp + 1
				countobsgpAAC = countobsgpAAC +1
				write(obsgp_str((countobsgp-1)*8+1:(countobsgp-1)*8+8),126) charC,'0',111,ia
				obsgp_count(countobsgp) = 1
				!RMSE(countobsgp) = resid(countobs2)**2
				
				Iobsgp(countobs2) = countobsgp
			endif

		else           ! point observables
			read(221,121) day,month,yr,charC,ip,charAZ,ia,obs,wt,AAPName(ia,1)	
			write(obsgp(countobs2),126) charC,'0',ip,ia	
			if (ia==08) resid(countobs2) = resid(countobs2) - (Z8mean(ip)-A8mean(ip))
			if (ia==09) resid(countobs2) = resid(countobs2) - (Z9mean(ip)-A9mean(ip))
			if (ia==10) resid(countobs2) = resid(countobs2) - (Z10mean(ip)-A10mean(ip))
			if (ia==11) resid(countobs2) = resid(countobs2) - (Z11mean(ip)-A11mean(ip))
			if (ia==12) resid(countobs2) = resid(countobs2) - (Z12mean(ip)-A12mean(ip))
			write(obsgp_temp,126) charC,'0',111,ia
			n_ip_ia(ip,ia) = n_ip_ia(ip,ia)+1
			nt(countobs2) = n_ip_ia(ip,ia)  ! size of (ip,ia) set to which countobs_th obs belongs
			if (INDEX(obsgp_str,obsgp_temp).gt.0) then
				obsgp_count((INDEX(obsgp_str,obsgp_temp)-1)/8+1) = obsgp_count((INDEX(obsgp_str,obsgp_temp)-1)/8+1) + 1
				!RMSE((INDEX(obsgp_str,obsgp_temp)-1)/8+1) = RMSE((INDEX(obsgp_str,obsgp_temp)-1)/8+1) + resid(countobs2)**2
				Iobsgp(countobs2) = (INDEX(obsgp_str,obsgp_temp)-1)/8+1
			else
				countobsgp = countobsgp + 1
				countobsgpAAP(ia) = countobsgpAAP(ia) +1
				write(obsgp_str((countobsgp-1)*8+1:(countobsgp-1)*8+8),126) charC,'0',111,ia
				obsgp_count(countobsgp) = 1
				!RMSE(countobsgp) = resid(countobs2)**2
				
				Iobsgp(countobs2) = countobsgp
			endif
		endif
	enddo

	99 close(221)


	DO k=1,countobs2
		j = Iobsgp(k)
		if (INDEX(obsgp(k),'PH').gt.0) then 
			read(obsgp(k),126) charP,charH,ip,ia
			RMSE(j) = RMSE(j)+resid(k)**2/real(n_ip_iah(ip,ia))
		elseif (INDEX(obsgp(k),'C').gt.0) then
			read(obsgp(k),126) charC,dum1,ic,ia
			RMSE(j) = RMSE(j)+resid(k)**2/real(n_ic_ia(ic,ia))
		elseif (INDEX(obsgp(k),'P').gt.0) then
			read(obsgp(k),126) charC,dum1,ip,ia
		    RMSE(j) = RMSE(j)+resid(k)**2/real(n_ip_ia(ip,ia))
		endif
	ENDDO

	DO k=1,countobs2
		j = Iobsgp(k)
		if (INDEX(obsgp(k),'PH').gt.0) then 
			read(obsgp(k),126) charP,charH,ip,ia
			weight(k) = (1/SQRT(real(n_ip_iah(ip,ia))))*SQRT(1./(REAL(countobsgp)*RMSE(j)))
		elseif (INDEX(obsgp(k),'C').gt.0) then
			read(obsgp(k),126) charC,dum1,ic,ia
			weight(k) = (1/SQRT(real(n_ic_ia(ic,ia))))*SQRT(1./(REAL(countobsgp)*RMSE(j)))
		elseif (INDEX(obsgp(k),'P').gt.0) then
			read(obsgp(k),126) charC,dum1,ip,ia
			weight(k) = (1/SQRT(real(n_ip_ia(ip,ia))))*SQRT(1./(REAL(countobsgp)*RMSE(j)))
		endif
	ENDDO

	write(6,*) 'writing ActualObs.Dat with new weights'

	open(unit=21,file='ActualObs.dat',STATUS='REPLACE',ACTION='WRITE')
	open(unit=223,file='ActualObs_ASCII.dat',STATUS='REPLACE',ACTION='WRITE')
	open(unit=221,file='ActualObs_anomaly.dat',STATUS='OLD',ACTION='READ')
	countobsgp=0
	countobs2=0
	countobsgpAAPh = 0
	countobsgpAAP = 0
	countobsgpAAC = 0
	do 
		read(221,*, end=999) dum
		countobs2=countobs2+1
		BACKSPACE(221)
		if (INDEX(dum,'AH').gt.0) then  ! hourly point observables
			read(221,1221) day,month,yr,charP,ip,charAZ,charH,ia,itstr,obs,wt,AAPhName(ia,1)
			write(21,1221) day,month,yr,charP,ip,charAZ,charH,ia,itstr,obs,weight((countobs2)),AAPhName(ia,1)
		elseif (INDEX(dum,'C').gt.0) then    ! catchment observables
			read(221,121) day,month,yr,charC,ic,charAZ,ia,obs,wt,AACName(ia,1)	
			write(21,121) day,month,yr,charC,ic,charAZ,ia,obs,weight((countobs2)),AACName(ia,1)
		else           ! point observables
			read(221,121) day,month,yr,charC,ic,charAZ,ia,obs,wt,AAPName(ia,1)		
			write(21,121) day,month,yr,charC,ic,charAZ,ia,obs,weight((countobs2)),AAPName(ia,1)
			write(223,*) day,month,yr,ic,ia,obs
		endif
	enddo

	999 close(21)
	close(223)
	close(221)
		

endif

! Write Output
if (JobFlag .eq. 1) then ! write obs (actual or model) at point then catchment scale
  yr = mod(TTDate%year,100)
  
  ! Write obs each step for PESTOutputFlag = 2 or 4  
if ((PESTOutputFlag .eq. 2) .or. (PESTOutputFlag .eq. 4)) then

! Point scale

! Write point obs (eg daily soil moisture)
 if (TTDate.ge.StartDate+365) then
   do ia = 1, naaP  
    if ((AAPflag(ia) .gt. 0.).and.(INDEX(AAPdaymth(ia),'d').gt.0)) then 
     do ip = 1, np
      if (AAP(ip,ia) .gt. -30.) then 
        select case (PESTOutputFlag)
        case (2)     
          if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
          if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*AAP(ip,ia),AAPrlow(ia))
			  write(21,121) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,AAP(ip,ia),1.0/err,AAPName(ia,1)
			  write(221,121) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,AAP(ip,ia),1.0/err,AAPName(ia,1)
			  write(str,122) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia 
			  write(22,*) str 
			  countobs = countobs + 1
			  resid(countobs) = ZZP(ip,ia) - AAP(ip,ia)
			  write(23,123) ZZP(ip,ia)
			  
        case (4)    
          write(23,123) ZZP(ip,ia)
        end select
      end if  
     end do
    end if 
   end do
endif ! (TTDate.ge.StartDate+365) 

! Write monthly mean of point obs
    do ia = 1, naaP        
     if ((AAPflag(ia) .gt. 0).and.(INDEX(AAPdaymth(ia),'m').gt.0)) then  
      do ip = 1, np
		if (AAPTAv(ip,ia,1).gt. -998.and.(EndMonth(TTDate))) then
         if (PESTOutputFlag .eq. 2) obs = AAPTAv(ip,ia,1)
         if ((PESTOutputFlag .eq. 3) .or. (PESTOutputFlag .eq. 4)) obs = ZZPTAv(ip,ia,1)         
         !ObsCount(ip,ia) = ObsCount(ip,ia) + 1    
         select case (PESTOutputFlag)
         case (2) 
           if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
           if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*AAP(ip,ia),AAPrlow(ia))         
           write(21,121) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,obs,1.0/err,AAPName(ia,1)
           write(str,122) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia 
           write(22,*) str 
           countobs = countobs + 1
		   resid(countobs) = ZZPTAv(ip,ia,1) - AAPTAv(ip,ia,1)
		   write(23,123) ZZPTAv(ip,ia,1)
		   !write(*,*)  ip, ia, ZZPTAv(ip,ia,1) , AAPTAv(ip,ia,1)
         case (3)
           if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
           if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*ZZP(ip,ia),AAPrlow(ia))         
           write(21,121) TTDate%day,TTDate%month,yr,charP,ic,charAZ,ia,obs,1.0/err,ZZPName(ia,1)
           write(str,122) TTDate%day,TTDate%month,yr,charP,ic,charAZ,ia 
           write(22,*) str  
           countobs = countobs + 1         
         case (4)                  
           write(23,123) obs
		   !write(*,*)  ip, ia, ZZPTAv(ip,ia,1) , AAPTAv(ip,ia,1)
         end select
 
         !ObsCount(ip,ia) = 0         
       end if   !(AAPTAv(ip,ia,1).gt. -998.and.(EndMonth(TTDate))
      end do !np
     end if  ! ((AAPflag(ia) .gt. 0).and.(AACdaymth.eq.'m'))
    end do  ! ia 

! Write annual mean of point obs
    do ia = 1, naaP        
     if ((AAPflag(ia) .gt. 0).and.(INDEX(AAPdaymth(ia),'y').gt.0)) then  
      do ip = 1, np
		if (AAPTAv(ip,ia,2).gt. -998.and.(EndYear(TTDate))) then
         if (PESTOutputFlag .eq. 2) obs = AAPTAv(ip,ia,2)
         if ((PESTOutputFlag .eq. 3) .or. (PESTOutputFlag .eq. 4)) obs = ZZPTAv(ip,ia,2)         
         !ObsCount(ip,ia) = ObsCount(ip,ia) + 1    
         select case (PESTOutputFlag)
         case (2) 
           if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
           if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*AAP(ip,ia),AAPrlow(ia))         
           write(21,121) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,obs,1.0/err,AAPName(ia,1)
           write(str,122) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia 
           write(22,*) str 
           countobs = countobs + 1
		   resid(countobs) = ZZPTAv(ip,ia,2) - AAPTAv(ip,ia,2)
		   write(23,123) ZZPTAv(ip,ia,2)
		   !write(*,*) TTDate, ip, ia, ZZP(ip,ia)
         case (3)
           if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
           if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*ZZP(ip,ia),AAPrlow(ia))         
           write(21,121) TTDate%day,TTDate%month,yr,charP,ic,charAZ,ia,obs,1.0/err,ZZPName(ia,2)
           write(str,122) TTDate%day,TTDate%month,yr,charP,ic,charAZ,ia 
           write(22,*) str  
           countobs = countobs + 1         
         case (4)                  
           write(23,123) obs
         end select
 
         !ObsCount(ip,ia) = 0         
       end if   !(AAPTAv(ip,ia,1).gt. -998.and.(EndYear(TTDate))
      end do !np
     end if  ! ((AAPflag(ia) .gt. 0).and.(AACdaymth.eq.'y'))
    end do  ! ia 

! Write whole of run mean of point obs
    do ia = 1, naaP        
     if ((AAPflag(ia) .gt. 0).and.(INDEX(AAPdaymth(ia),'r').gt.0)) then  
      do ip = 1, np
		if (AAPTAv(ip,ia,3).gt. -998.and.((TTDate .ge. EndDate))) then
         if (PESTOutputFlag .eq. 2) obs = AAPTAv(ip,ia,3)
         if ((PESTOutputFlag .eq. 3) .or. (PESTOutputFlag .eq. 4)) obs = ZZPTAv(ip,ia,3)         
         !ObsCount(ip,ia) = ObsCount(ip,ia) + 1    
         select case (PESTOutputFlag)
         case (2) 
           if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
           if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*AAP(ip,ia),AAPrlow(ia))         
           write(21,121) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,obs,1.0/err,AAPName(ia,1)
           write(str,122) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia 
           write(22,*) str 
           countobs = countobs + 1
		   resid(countobs) = ZZPTAv(ip,ia,3) - AAPTAv(ip,ia,3)
		   write(23,123) ZZPTAv(ip,ia,3)
		   !write(*,*)  ip, ia,  ZZPTAv(ip,ia,3) , AAPTAv(ip,ia,3)
         case (3)
           if (AAPabsrel(ia) .eq. 'a') err = AAPr(ia)
           if (AAPabsrel(ia) .eq. 'r') err = max(AAPr(ia)*ZZP(ip,ia),AAPrlow(ia))         
           write(21,121) TTDate%day,TTDate%month,yr,charP,ic,charAZ,ia,obs,1.0/err,ZZPName(ia,2)
           write(str,122) TTDate%day,TTDate%month,yr,charP,ic,charAZ,ia 
           write(22,*) str  
           countobs = countobs + 1         
         case (4)                  
           write(23,123) obs
		   !write(*,*)  ip, ia,  ZZPTAv(ip,ia,3) , AAPTAv(ip,ia,3)
         end select
 
         !ObsCount(ip,ia) = 0         
       end if   !(AAPTAv(ip,ia,1).gt. -998.and.((TTDate .ge. EndDate))
      end do !np
     end if  ! ((AAPflag(ia) .gt. 0).and.(AACdaymth.eq.'r'))
    end do  ! ia 
  
! Point scale

! Write monthly mean of hourly point obs
	 if  (TTDate%day.eq.1) then
		PredAveh = 0.0
		ObsAveh= 0.0
		ObsCounth = 0  		
	 endif

 if (TTDate.ge.StartDate+365) then
    do ia = 1, naaPh      
     if (AAPhflag(ia) .gt. 0) then  
      do ip = 1, np
		do it = 1,ntime
		   if (AAPh(ip,ia,it) .gt. -998) then

			 if (PESTOutputFlag .eq. 2) then
				ObsAveh(ip,ia,it) = ObsAveh(ip,ia,it) + AAPh(ip,ia,it)
				PredAveh(ip,ia,it) = PredAveh(ip,ia,it) + ZZPh(ip,ia,it)

			 endif
			 if ((PESTOutputFlag .eq. 3) .or. (PESTOutputFlag .eq. 4)) ObsAveh(ip,ia,it) = ObsAveh(ip,ia,it) + ZZPh(ip,ia,it)           
			 ObsCounth(ip,ia,it) = ObsCounth(ip,ia,it) + 1  
		   end if
		   if ((EndMonth(TTDate)) .and. (ObsCounth(ip,ia,it) .gt. 0)) then  
			obs = ObsAveh(ip,ia,it)/ObsCounth(ip,ia,it)
			 select case (PESTOutputFlag)
			 case (2) 
			   pred = PredAveh(ip,ia,it)/ObsCounth(ip,ia,it)
			   if (AAPhabsrel(ia) .eq. 'a') err = AAPhr(ia)
			   if (AAPhabsrel(ia) .eq. 'r') err = max(AAPhr(ia)*AAPh(ip,ia,it),AAPhrlow(ia))     
			   if ((it-1).lt.10) then
				write(itstr,"(1(a1,i1))") '0',it-1
			   else
				write(itstr,"(1(i2))") it-1
			   endif
			   
				   write(21,1221) TTDate%day,TTDate%month,yr,charP,ip,charAZ,charH,ia,itstr,obs,1.0/err,AAPhName(ia,1)
				   write(221,1221) TTDate%day,TTDate%month,yr,charP,ip,charAZ,charH,ia,itstr,obs,1.0/err,AAPhName(ia,1)

				   write(str,1222) TTDate%day,TTDate%month,yr,charP,ip,charAZ,charH,ia,itstr
				   write(22,*) str 
				   countobs = countobs + 1  
				   resid(countobs) = pred - obs
				   write(23,123) obs
		    
			 case (4) 
                 
			   write(23,123) obs

			 end select
			 PredAveh(ip,ia,it) = 0.0
			 ObsAveh(ip,ia,it) = 0.0
			 ObsCounth(ip,ia,it) = 0         
		   end if  
	   end do 
      end do
     end if 
    end do 
endif

! Write hourly point obs

!if (present(AAPh)) then
!   do ia = 1, naaPh  
!    if (AAPhflag(ia) .gt. 0) then 
!     do ip = 1, np
!		do it = 1,ntime
!			if (AAPh(ip,ia,it) .gt. -998) then 
!			select case (PESTOutputFlag)
!			case (2)     
!				if (AAPhabsrel(ia) .eq. 'a') err = AAPhr(ia)
!				if (AAPhabsrel(ia) .eq. 'r') err = max(AAPhr(ia)*AAPh(ip,ia,it),AAPhrlow(ia))
!				if (TTDate.ge.StartDate+365) then
!					write(21,124) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,it,AAPh(ip,ia,it),1.0/err,AAPhName(ia,1)
!					write(str,1222) TTDate%day,TTDate%month,yr,charP,ip,charAZ,ia,it
!					write(22,*) str 
!					countobs = countobs + 1
!				endif
!			case (4) 
!			if (TTDate.ge.StartDate+365) then    
!				write(23,123) ZZPh(ip,ia,it)
!			endif
!			end select
!			end if  
!		end do
!	 end do
!    end if 
!   end do
!endif





! Catchments (annual)
 if (EndYear(TTDate)) then
   do ia = 1, naaC  
    if (AACflag(ia) .gt. 0.and.(INDEX(AACdaymth(ia),'y').gt.0)) then 
     do ic = 1, nc
      if (AACTav(ic,ia,2) .gt. -998) then 
        select case (PESTOutputFlag)
        case (2)     
          if (AACabsrel(ia) .eq. 'a') err = AACr(ia)
          if (AACabsrel(ia) .eq. 'r') err = max(AACr(ia)*AACTav(ic,ia,2),AACrlow(ia))
          write(21,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,AACTav(ic,ia,2),1.0/err,AACName(ia,1)
		  write(221,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,AACTav(ic,ia,2),1.0/err,AACName(ia,1)
          write(str,122) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia 
          write(22,*) str 
          countobs = countobs + 1
		  resid(countobs) = log(ZZCTav(ic,ia,2)) - log(AACTav(ic,ia,2))
		  write(23,123) ZZCTAv(ic,ia,2)
        case (4)     
          write(23,123) ZZCTAv(ic,ia,2)
        end select
      end if  
     end do
    end if 
   end do    
 endif   

iCatchMap(:) = nint(CatchMap(:))        ! integer catchment map (length nc)
 ! Catchments (whole of run) : output mean Precip - mean discharge
 if ((TTDate .ge. EndDate)) then
   do ia = 1, naaC  
    if (AACflag(ia) .gt. 0.and.(INDEX(AACdaymth(ia),'r').gt.0)) then 
     do ic = 1, nc
	 
	  FWPrecC = sum( FFTAV(:,3,3),  mask = ( (iCatchMap==iCatchID(ic)) ) ) /    &
                 count( mask = ( (iCatchMap==iCatchID(ic)) ) )

      if (AACTav(ic,ia,3) .gt. -998) then 
        select case (PESTOutputFlag)
        case (2)     
          if (AACabsrel(ia) .eq. 'a') err = AACr(ia)
          if (AACabsrel(ia) .eq. 'r') err = max(AACr(ia)*AACTav(ic,ia,3),AACrlow(ia))
          !write(21,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,log(max(AACTav(ic,ia,3),1.e-10)),1.0/err,AACName(ia,1)
		  !write(221,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,log(max(AACTav(ic,ia,3),1.e-10)),1.0/err,AACName(ia,1)
		  write(21,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia, FWPrecC-AACTav(ic,ia,3)/30.4375,1.0/err,AACName(ia,1)
		  write(221,121) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia,FWPrecC-AACTav(ic,ia,3)/30.4375,1.0/err,AACName(ia,1)
          write(str,122) TTDate%day,TTDate%month,yr,charC,ic,charAZ,ia 
          write(22,*) str 
          countobs = countobs + 1
		  !resid(countobs) = log(max(ZZCTav(ic,ia,3),1.e-10)) - log(max(AACTav(ic,ia,3),1.e-10))
		  !write(23,123) log(max(ZZCTAv(ic,ia,3),1.e-10))
		  resid(countobs) = (FWPrecC-ZZCTav(ic,ia,3)/30.4375)-(FWPrecC-AACTav(ic,ia,3)/30.4375)
		  write(23,123) FWPrecC-AACTav(ic,ia,3)/30.4375
		  !write(*,*) (FWPrecC-ZZCTav(ic,ia,3)/30.4375), (FWPrecC-AACTav(ic,ia,3)/30.4375)
        case (4)     
          !write(23,123) log(max(ZZCTAv(ic,ia,3),1.e-10))
		  write(23,123) FWPrecC-ZZCTav(ic,ia,3)/30.4375
		  !write(*,*) (FWPrecC-ZZCTav(ic,ia,3)/30.4375), (FWPrecC-AACTav(ic,ia,3)/30.4375)
        end select
      end if  
     end do
    end if 
   end do    
 endif   
  
  121 format(i2.2,i2.2,i2.2,a1,i3.3,a1,i2.2,'  ',e13.5,'    ',e13.5,'    ',a8) ! for pst file
  122 format('l1 [',i2.2,i2.2,i2.2,a1,i3.3,a1,i2.2,']1:14')      ! for ins file
  1221 format(i2.2,i2.2,i2.2,a1,i3.3,a1,a1,i1,a2'  ',e13.5,'    ',e13.5,'    ',a9) ! for pst file (hourly obs)
  1222 format('l1 [',i2.2,i2.2,i2.2,a1,i3.3,a1,a1,i1,a2']1:14')      ! for ins file (hourly obs)
  124 format(i2.2,i2.2,i2.2,a1,i3.3,a1,i1,i3.3'  ',e13.5,'    ',f8.2,'    ',a8)  ! for pst file
  125 format('l1 [',i2.2,i2.2,i2.2,a1,i3.3,a1,i1,i3.3']1:17')      ! for ins file
  123 format(e13.5)  
  126 format(a1,a1,i3.3,i3.3) ! for obsgp string
  
 end if 
  
end if  ! jobflag .eq. 1 (i.e. write output)

END SUBROUTINE WritePESTOutput

!###############################################################################

END MODULE WaterDynModule

!###############################################################################
!###############################################################################

PROGRAM WaterDynM
!-------------------------------------------------------------------------------
USE TypeDef
USE DateFunctions
USE WaterDynModule
USE PointerModule

USE define_dimensions                ! for CABLE
USE main_cable_module				     ! for CABLE
USE checks_module					 ! for CABLE
USE Physical_constants				 ! for CABLE
USE parameter_module				 ! for CABLE
USE canopy_module					 ! for CABLE
USE define_types					 ! for CABLE
USE constants				   		! for CABLE


implicit none
real(sp)            :: Time0        ! time (in host model)
real(sp)            :: Time0_met     ! time of met (in host model)
real(sp)            :: UU(nuu)      ! spatially uniform parameters (adjustable)
real(sp) ,allocatable:: VV(:,:)      ! spatially variable parameters (adjustable)
real(sp),allocatable:: XX0(:,:)     ! XX0(np,nxx) = state variables at T0
real(sp),allocatable:: XX1(:,:)     ! XX1(np,nxx) = state variables at T1
real(sp),allocatable:: ZZP1(:,:)    ! ZZP1(np,nzzP) = predicted point obs at T1
real(sp),allocatable:: ZZPh1(:,:,:)    ! ZZPh1(np,nzzP,ntime) = predicted point obs at T1 (hourly)
real(sp),allocatable:: ZZC1(:,:)    ! ZZC1(nc,nzzC) = predicted catchment obs at T1
real(sp),allocatable:: AAP1(:,:)    ! AAP1(np,naaP) = actual point obs at T1
real(sp),allocatable:: AAPh1(:,:,:)    ! AAP1(np,naaP,ntime) = actual point obs at T1 (hourly)
real(sp),allocatable:: AAC1(:,:)    ! AAC1(nc,naaC) = actual catchment obs at T1
real(sp)            :: TTime(ntt)   ! model time variables
type(dmydate):: TTDate			    ! current date

integer(i4b):: it, is                  ! counters
! Variables specific to WDsens:
real(sp)            :: UURef(nuu)       ! spatially uniform parameters (adjustable)
real(sp),allocatable:: VVRef(:,:)       ! spatially variable parameters (adjustable)
real(sp),allocatable:: Resp(:)          ! response of model (value of a given flux) with perturbed UU
real(sp),allocatable:: RespUmin(:)          ! response of model (value of a given flux) with min UU 
real(sp),allocatable:: RespUmax(:)          ! response of model (value of a given flux) with max UU
real(sp),allocatable:: Sens(:,:)        ! sensitivity of param = P/Resp * d(Resp)/dP
real(sp),allocatable:: RespRef(:)       ! model response for default parameters
real(sp),allocatable:: RespRefMonth(:,:)       ! monthly model response for default parameters
real(sp),allocatable:: RespMonth(:,:)       ! monthly model response for perturbed parameters
real(sp),allocatable:: SensMonth(:,:,:)       ! monthly sensitivity of param = P/Resp * d(Resp)/dP
character(20),allocatable:: names(:)    ! response flux names (for writing to output)
integer(i4b):: ip,im,ir,iuu,ivv,iff,nuusens,nvvsens,nffsens,nuucount,nvvcount,nffcount,nMonthCount

real(sp),allocatable:: minsm0008(:), maxsm0008(:), minsm0090(:), maxsm0090(:)
real(sp):: obsval, obssig
character(30):: obsstr, obsname, dum
integer(i4b):: A9count=0, A8count=0
real(sp):: A8mean=0.0, A9mean=0.0
integer(i4b):: Z9count=0, Z8count=0
real(sp):: Z8mean=0.0, Z9mean=0.0
!-------------------------------------------------------------------------------

CALL CPU_TIME (cpuTimeBegin)

! FORWARD RUN
! * Initialise WaterDyn
CALL ReadControlFile

! * Allocate arrays declared in main program
if (.not.allocated(VV))  allocate (VV(np,nVV))    !  spatially variable parameters 
if (.not.allocated(VVref))  allocate (VVref(np,nVV))    !  spatially variable parameters 
if (.not.allocated(XX0))  allocate (XX0(np,nxx))    ! state variables at T0
if (.not.allocated(XX1))  allocate (XX1(np,nxx))    ! state variables at T1
if (.not.allocated(ZZP1)) allocate (ZZP1(np,nzzP))  ! predicted point obs at T1
if (.not.allocated(ZZPh1)) allocate (ZZPh1(np,nzzPh,ntime))  ! predicted point obs at T1
if (.not.allocated(ZZC1)) allocate (ZZC1(nc,nzzC))  ! predicted catchment obs at T1
if (.not.allocated(AAP1)) allocate (AAP1(np,naaP))  ! actual point obs at T1
if (.not.allocated(AAPh1)) allocate (AAPh1(np,naaPh,ntime))  ! actual point obs at T1
if (.not.allocated(AAC1)) allocate (AAC1(nc,naaC))  ! actual catchment obs at T1
if (.not.allocated(minsm0008)) allocate (minsm0008(np))
if (.not.allocated(maxsm0008)) allocate (maxsm0008(np))
if (.not.allocated(minsm0090)) allocate (minsm0090(np))
if (.not.allocated(maxsm0090)) allocate (maxsm0090(np))

! * Assign pointers to generic arrays
CALL PointAll (TargetTTime=TTime)



if (SensTestFlag ==3 ) then
	nffsens = 0
	do iff = 1, nff
	  if (FFSensFlag(iff) .eq. 1) nffsens = nffsens + 1
	end do 
	nvvsens = 0
	do ivv = 1, nvv
	  if (VVSensFlag(ivv) .eq. 1) nvvsens = nvvsens + 1
	end do
	  
	if (.not.allocated(RespRefMonth)) allocate (RespRefMonth(nMonths,nffsens))  ! monthly model response for default parameters
	if (.not.allocated(RespMonth)) allocate (RespMonth(nMonths,nffsens))  ! monthly model response for default parameters
	if (.not.allocated(SensMonth)) allocate (SensMonth(nMonths,nffsens,nvvsens))  ! monthly model response for default parameters
    if (.not.allocated(names)) allocate(names(nffsens))
endif
! * Initialise quantities with host-program scope (VV, UU, XX0)
CALL InitVV(VVref)
VV=VVref
UU = UUinitCtl
CALL InitStores (XX0)
call PointAll(TargetZZP=ZZP1)
ispin =1
InitFlagCO2 = 0 

if (nspin>0) then
! spin-up with repeated meteorology
it = -1 ! time-step counter for whole run
is = -1 ! time-step counter for current spin-cycle
	do ispin = 1,nspin
		if (ispin>1) InitFlagCO2 = 1 ! don't initialise CO2 file unless first spin 
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1             ! Flags last timestep in run
		nMonthCount = 0
		do is=0,nSteps_spin	 
		  if (is.eq.nSteps_spin) EndFlag = 0
		  if (nsteps_spin.eq.nsteps_run) then
			it=is
		  else
			it = it+1  ! increment time-step counter for whole run
		  endif
		  Time0 = real(it) * 1.0
		  Time0_met = real(is) * 1.0
		  StartDate=StartDate_spin
		  nsteps=nsteps_spin
		  CALL WaterDynStep (Time0,Time0_met, UU,VV, XX0,    XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		end do
	end do
	TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
	if ((EndDate.gt.TTDate)) then
	    InitFlagCO2 = 1
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1             ! Flags last timestep in run
		TTDate = dmydate(nint(TTDay), nint(TTMonth), nint(TTYear))
		StartDate=TTDate+1
		EndDate = EndDate_run
		nsteps = nsteps_run-nspin*nsteps_spin
		nmonths = nmonths_run - nspin*nmonths_spin
		do is=0,nSteps		  
		  it = it+1  ! increment time-step counter for whole run
		  if (it.eq.nSteps_run) EndFlag = 0
		  Time0 = real(it) * 1.0
		  Time0_met = real(is) * 1.0
		  CALL WaterDynStep (Time0,Time0_met, UU,VV, XX0,    XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		enddo
	endif
endif

!  model run (no spin-up)
if (nspin==0) then

	EndDate = EndDate_run
	StartDate=StartDate_run
	nsteps=nsteps_run
	nmonths = nmonths_run

	InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
	EndFlag = 1             ! Flags last timestep in run
	nMonthCount = 0
	do it=0,nSteps
	  Time0 = real(it) * 1.0
	  Time0_met=Time0
	  if (it.eq.nSteps) EndFlag = 0
	  CALL WaterDynStep (Time0,Time0_met, UU,VV, XX0,    XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
	  if (it.eq.1) then
		minsm0008 = Zsm0008
		maxsm0008 = Zsm0008
		minsm0090 = Zsm0090
		maxsm0090 = Zsm0090
	  else
		minsm0008 = min(Zsm0008,minsm0008)
		maxsm0008 = max(Zsm0008,maxsm0008)
		minsm0090 = min(Zsm0090,minsm0090)
		maxsm0090 = max(Zsm0090,maxsm0090)
	  endif

	  XX0 = XX1
	  if(SensTestFlag ==3 .and. EndMonth(StartDate + iT)) then
		nffcount = 0
		nMonthCount = nMonthCount+1
  		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			RespRefMonth(nMonthCount,nffcount) = FFSpTAvTemp(1,iff,1) 
			names(nffcount) = FFName(iff,1)
		  end if  
		end do
	   endif
	end do
endif

CALL CPU_TIME (cpuTimeNow)
if (DiagsFlag == 1) then
  write(*,"('Forward runtime (sec)     =',f12.2)") cpuTimeNow - cpuTimeBegin
end if

if (SensTestFlag == 0) write(*,*) "End BIOS2"

if (SensTestFlag == 0) STOP

! SENSITIVITY TESTS ABOUT REFERENCE CASE (only if SensTestFlag = 1 or 2)
! * reset output flags to turn off any further output
DiagsFlag         = 0
DomTSOutputFlag   = 0
CatchTSOutputFlag = 0
MapOutputFlag     = 0
PESTOutputFlag    = 0
nEnsemble         = 1   ! used in test for writing run info to screen in ReadControlFile

if (SensTestFlag.eq.1) then
	! * Count cases
	nuusens = 0
	do iuu = 1, nuu
	  if (UUSensFlag(iuu) .eq. 1) nuusens = nuusens + 1
	end do  
	nffsens = 0
	do iff = 1, nff
	  if (FFSensFlag(iff) .eq. 1) nffsens = nffsens + 1
	end do   
	write(*,"(/,'START PERTURBATION RUNS FOR SENSITIVITY TEST')")
	write(*,"('Parameters varied         =',i10)") nuusens
	write(*,"('Flux sensitivites         =',i10)") nffsens
	write(*,"('Model runs (inc ref case) =',i10)") nuusens+1

	! * allocate arrays for sensitivity tests
	allocate(Resp(nffsens))
	allocate(names(nffsens))
	allocate(RespRef(nffsens))  
	allocate(Sens(nffsens,nuusens))

	! * Save reference-case parameters and outputs, as UURef and RespRef(nffsens)
	UURef = UU              ! save reference UU
	nffcount = 0
	do iff = 1,nff  
	  if (FFSensFlag(iff) .eq. 1) then 
		nffcount = nffcount + 1
		RespRef(nffcount) = FFSpTAv(1,iff,3) 
		names(nffcount) = FFName(iff,1)
	  end if  
	end do
	write(*,"(/,'      Flux      RefValue')")
	write(*,"(a10,2x,e12.4)") (trim(names(iff)),RespRef(iff),iff=1,nffsens)
  
	! Open output files (in OutputSens directory)
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'RawFluxes.csv'
	OPEN (unit=28, file=trim(GenFileName))
	write(28,"('Domain:,,',  a)") trim(DomName)
	write(28,"('Run ID:,,',  a)") trim(RunID)
	write(28,"('Comment:,,', a)") trim(Comment)
	write(28,"(20(a10,','))") 'Fluxes,  ',names(:)
	write(28,"(a15,' ',20(e13.6,','))") 'Default,   , ',RespRef
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'Sensitivity.csv'
	OPEN (unit=29, file=trim(GenFileName))
	write(29,"('Domain:,,',  a)") trim(DomName)
	write(29,"('Run ID:,,',  a)") trim(RunID)
	write(29,"('Comment:,,', a)") trim(Comment)
	write(29,"(20(a10,','))") 'Fluxes ',names(:)

	! Loop over parameters varying one at a time
	do ip = 1, nuu
	  nuucount = 0
	  if (UUSensFlag(ip) .eq. 1) then
		write(*,"('Vary parameter ',a16)") trim(UUName(ip,1))
		nuucount = nuucount + 1
		UU = UUinitCtl          ! Initialise quantities with global scope (XX0, UU)
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1
		UU(ip) = UU(ip)*1.01    ! perturb reference UU
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  XX0 = XX1
		end do
		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			RespUmin(nffcount) = FFSpTAv(1,iff,3)        
			Sens(nffcount,nuucount) = UURef(ip)/RespRef(nffcount) *         &
			  (Resp(nffcount)-RespRef(nffcount))/(UU(ip)-UURef(ip)) 
		  end if  
		end do
		CALL CPU_TIME (cpuTimeNow)
		write(*,"('Done parameter ',a16,'      CPUtime (sec) = ',f9.2)")    &
		  trim(UUName(ip,1)), cpuTimeNow - cpuTimeBegin
	!    write(6,*) ' Resp =',Resp(:)
	!    write(6,*) ' Sens =',Sens(:,nuucount)
	!    write(6,*) ' '
		write(28,"(a10,', ',20(e13.6,','))") UUName(ip,1),UURef(ip),Resp(:)    
		write(29,"(a10,', ',20(e13.6,','))") UUName(ip,1),Sens(:,nuucount)
	  end if   
	end do 
elseif (SensTestFlag.eq.4) then
	! * Count cases
	nuusens = 0
	do iuu = 1, nuu
	  if (UUSensFlag(iuu) .eq. 1) nuusens = nuusens + 1
	end do  
	nffsens = 0
	do iff = 1, nff
	  if (FFSensFlag(iff) .eq. 1) nffsens = nffsens + 1
	end do   
	write(*,"(/,'START PERTURBATION RUNS FOR SENSITIVITY TEST 4')")
	write(*,"('Parameters varied         =',i10)") nuusens
	write(*,"('Flux sensitivites         =',i10)") nffsens
	write(*,"('Model runs (inc ref case) =',i10)") 2.*nuusens+1

	! * allocate arrays for sensitivity tests
	allocate(RespUmin(nffsens))
	allocate(RespUmax(nffsens))
	allocate(names(nffsens))
	allocate(RespRef(nffsens))  
	allocate(Sens(nffsens,nuusens))

	! * Save reference-case parameters and outputs, as UURef and RespRef(nffsens)
	UURef = UU              ! save reference UU
	nffcount = 0
	do iff = 1,nff  
	  if (FFSensFlag(iff) .eq. 1) then 
		nffcount = nffcount + 1
		RespRef(nffcount) = FFSpTAv(1,iff,3) 
		names(nffcount) = FFName(iff,1)
	  end if  
	end do
	write(*,"(/,'      Flux      RefValue')")
	write(*,"(a10,2x,e12.4)") (trim(names(iff)),RespRef(iff),iff=1,nffsens)
  
	! Open output files (in OutputSens directory)
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'RawFluxes.csv'
	OPEN (unit=28, file=trim(GenFileName))
	write(28,"('Domain:,,',  a)") trim(DomName)
	write(28,"('Run ID:,,',  a)") trim(RunID)
	write(28,"('Comment:,,', a)") trim(Comment)//' '//'SENSITIVITY TEST 4'
	write(28,"(20(a10,','))") 'Fluxes,  ',names(:)
	write(28,"(a15,' ',20(e13.6,','))") 'Default,   , ',RespRef
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'Sensitivity.csv'
	OPEN (unit=29, file=trim(GenFileName))
	write(29,"('Domain:,,',  a)") trim(DomName)
	write(29,"('Run ID:,,',  a)") trim(RunID)
	write(29,"('Comment:,,', a)") trim(Comment)
	write(29,"(20(a10,','))") 'Fluxes ',names(:)

	! Loop over parameters varying one at a time
	do ip = 1, nuu
	  nuucount = 0
	  if (UUSensFlag(ip) .eq. 1) then
		write(*,"('Vary parameter ',a16)") trim(UUName(ip,1))
		nuucount = nuucount + 1
		UU = UUinitCtl          ! Initialise quantities with global scope (XX0, UU)
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1
		UU(ip) = UUminctl(ip)   ! set UU to minimum value
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  XX0 = XX1
		end do
		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			RespUmin(nffcount) = FFSpTAv(1,iff,3)        
		  end if  
		end do

		UU = UUinitCtl          ! Initialise quantities with global scope (XX0, UU)
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1
		UU(ip) = UUmaxctl(ip)   ! set UU to minimum value
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  XX0 = XX1
		end do
		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			RespUmax(nffcount) = FFSpTAv(1,iff,3)        
		  end if  
		end do

		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1  
			Sens(nffcount,nuucount) = (RespUmax(nffcount)-RespUmin(nffcount))/RespRef(nffcount)
		  end if  
		end do

		CALL CPU_TIME (cpuTimeNow)
		write(*,"('Done parameter ',a16,'      CPUtime (sec) = ',f9.2)")    &
		  trim(UUName(ip,1)), cpuTimeNow - cpuTimeBegin
	!    write(6,*) ' Resp =',Resp(:)
	!    write(6,*) ' Sens =',Sens(:,nuucount)
	!    write(6,*) ' '
		write(28,"(a10,', ',20(e13.6,','))") UUName(ip,1),UURef(ip),Resp(:)    
		write(29,"(a10,', ',20(e13.6,','))") UUName(ip,1),Sens(:,nuucount)
	  end if   
	end do 

elseif (SensTestFlag.eq.2) then
	! * Count cases
	nvvsens = 0
	do ivv = 1, nvv
	  if (VVSensFlag(ivv) .eq. 1) nvvsens = nvvsens + 1
	end do  
	nffsens = 0
	do iff = 1, nff
	  if (FFSensFlag(iff) .eq. 1) nffsens = nffsens + 1
	end do   
	write(*,"(/,'START PERTURBATION RUNS FOR SENSITIVITY TEST')")
	write(*,"('Parameters varied         =',i10)") nvvsens
	write(*,"('Flux sensitivites         =',i10)") nffsens
	write(*,"('Model runs (inc ref case) =',i10)") nvvsens+1

	! * allocate arrays for sensitivity tests
	allocate(Resp(nffsens))
	allocate(names(nffsens))
	allocate(RespRef(nffsens))  
	allocate(Sens(nffsens,nvvsens))

	! * Save reference-case parameters and outputs, as UURef and RespRef(nffsens)
	
	nffcount = 0
	do iff = 1,nff  
	  if (FFSensFlag(iff) .eq. 1) then 
		nffcount = nffcount + 1
		RespRef(nffcount) = FFSpTAv(1,iff,3) 
		names(nffcount) = FFName(iff,1)
	  end if  
	end do
	write(*,"(/,'      Flux      RefValue')")
	write(*,"(a10,2x,e12.4)") (trim(names(iff)),RespRef(iff),iff=1,nffsens)
  
	! Open output files (in OutputSens directory)
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '\Run' //    &
				  trim(RunID) // '/OutputSens/' // 'RawFluxes.csv'
	OPEN (unit=28, file=trim(GenFileName))
	write(28,"('Domain:,,',  a)") trim(DomName)
	write(28,"('Run ID:,,',  a)") trim(RunID)
	write(28,"('Comment:,,', a)") trim(Comment)
	write(28,"(20(a10,','))") 'Fluxes,  ',names(:)
	write(28,"(a15,' ',20(e13.6,','))") 'Default,   , ',RespRef
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'Sensitivity.csv'
	OPEN (unit=29, file=trim(GenFileName))
	write(29,"('Domain:,,',  a)") trim(DomName)
	write(29,"('Run ID:,,',  a)") trim(RunID)
	write(29,"('Comment:,,', a)") trim(Comment)
	write(29,"(20(a10,','))") 'Fluxes ',names(:)

	! Loop over parameters varying one at a time
	do ip = 1, nvv
	  nvvcount = 0
	  if (VVSensFlag(ip) .eq. 1) then
		write(*,"('Vary parameter ',a16)") trim(VVName(ip,1))
		nvvcount = nvvcount + 1
		VV=VVref
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		Endflag=1
		VV(:,ip) = VV(:,ip)*1.01    ! perturb reference VV
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  XX0 = XX1
		end do
		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			Resp(nffcount) = FFSpTAv(1,iff,3)        
			Sens(nffcount,nvvcount) = VVRef(1,ip)/RespRef(nffcount) *         &
			  (Resp(nffcount)-RespRef(nffcount))/(VV(1,ip)-VVRef(1,ip)) 
		  end if  
		end do
		CALL CPU_TIME (cpuTimeNow)
		write(*,"('Done parameter ',a16,'      CPUtime (sec) = ',f9.2)")    &
		  trim(VVName(ip,1)), cpuTimeNow - cpuTimeBegin
	!    write(6,*) ' Resp =',Resp(:)
	!    write(6,*) ' Sens =',Sens(:,nuucount)
	!    write(6,*) ' '
		write(28,"(a10,', ',20(e13.6,','))") VVName(ip,1),VVRef(1,ip),Resp(:)    
		write(29,"(a10,', ',20(e13.6,','))") VVName(ip,1),Sens(:,nvvcount)
	  end if   
	end do 

elseif (SensTestFlag.eq.3) then
	! * Count cases
	nvvsens = 0
	do ivv = 1, nvv
	  if (VVSensFlag(ivv) .eq. 1) nvvsens = nvvsens + 1
	end do  
	nffsens = 0
	do iff = 1, nff
	  if (FFSensFlag(iff) .eq. 1) nffsens = nffsens + 1
	end do   
	write(*,"(/,'START PERTURBATION RUNS FOR SENSITIVITY TEST')")
	write(*,"('Parameters varied         =',i10)") nvvsens
	write(*,"('Flux sensitivites         =',i10)") nffsens
	write(*,"('Model runs (inc ref case) =',i10)") nvvsens+1



	! * Save reference-case parameters and outputs, as UURef and RespRef(nffsens)

	!write(*,"(/,'      Flux      RefValue')")
	!write(*,"(a10,2x,e12.4)") (trim(names(iff)),RespRef(iff),iff=1,nffsens)
  
	! Open output files (in OutputSens directory)
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'RawFluxes.csv'
	OPEN (unit=28, file=trim(GenFileName))
	write(28,"('Domain:,,',  a)") trim(DomName)
	write(28,"('Run ID:,,',  a)") trim(RunID)
	write(28,"('Comment:,,', a)") trim(Comment)
	write(28,"(200(a10,','))") ' ','nMonth',TTName(1:ntt-1,1),'Default',names(:)
!	write(28,"(a15,' ',20(e13.6,','))") 'Default,   , ',RespRef
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'Sensitivity.csv'
	OPEN (unit=29, file=trim(GenFileName))
	write(29,"('Domain:,,',  a)") trim(DomName)
	write(29,"('Run ID:,,',  a)") trim(RunID)
	write(29,"('Comment:,,', a)") trim(Comment)
	write(29,"(200(a10,','))") ' ','nMonth',TTName(1:ntt-1,1),names(:)

	! Loop over parameters varying one at a time
	do ip = 1, nvv
	  nvvcount = 0
	  if (VVSensFlag(ip) .eq. 1) then
		write(*,"('Vary parameter ',a16)") trim(VVName(ip,1))
		nvvcount = nvvcount + 1
		VV=VVref
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		Endflag=1
		nMonthCount=0
		VV(:,ip) = VV(:,ip)*1.01    ! perturb reference UU
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  XX0 = XX1
		  if(EndMonth(StartDate + iT)) then
			nffcount = 0
			nMonthCount = nMonthCount+1
			do iff = 1,nff  
				if (FFSensFlag(iff) .eq. 1) then 
					nffcount = nffcount + 1
					RespMonth(nMonthCount,nffcount) = FFSpTAvTemp(1,iff,1) 
					SensMonth(nMonthCount,nffcount,nvvcount) = VVRef(1,ip)/RespRefMonth(nMonthCount,nffcount) *         &
			        (RespMonth(nMonthCount,nffcount)-RespRefMonth(nMonthCount,nffcount))/(VV(1,ip)-VVRef(1,ip)) 
				end if  
			end do
			write(28,"(a10,', ',i4,', ',20(e13.6,','))") VVName(ip,1),nMonthCount,TTime(1:ntt-1),VVRef(1,ip),RespMonth(nMonthCount,:)    
		    write(29,"(a10,', ',i4,', ',20(e13.6,','))") VVName(ip,1),nMonthCount,TTime(1:ntt-1),SensMonth(nMonthCount,:,nvvcount)


		  endif
		end do
	
		CALL CPU_TIME (cpuTimeNow)
		write(*,"('Done parameter ',a16,'      CPUtime (sec) = ',f9.2)")    &
		  trim(VVName(ip,1)), cpuTimeNow - cpuTimeBegin
	
		
	  end if   
	end do 

elseif (SensTestFlag.eq.5) then   ! include soil moisture range in sensitivity analysis
	! * Count cases
	nuusens = 0
	do iuu = 1, nuu
	  if (UUSensFlag(iuu) .eq. 1) nuusens = nuusens + 1
	end do  
	nffsens = 0
	do iff = 1, nff
	  if (FFSensFlag(iff) .eq. 1) nffsens = nffsens + 1
	end do   
	write(*,"(/,'START PERTURBATION RUNS FOR SENSITIVITY TEST 5')")
	write(*,"('Parameters varied         =',i10)") nuusens
	write(*,"('Flux sensitivites         =',i10)") nffsens+2
	write(*,"('Model runs (inc ref case) =',i10)") 2.*nuusens+1

	! * allocate arrays for sensitivity tests
	allocate(RespUmin(nffsens+2))
	allocate(RespUmax(nffsens+2))
	allocate(names(nffsens+2))
	allocate(RespRef(nffsens+2))  
	allocate(Sens(nffsens+2,nuusens))

	! * Save reference-case parameters and outputs, as UURef and RespRef(nffsens)
	UURef = UU              ! save reference UU
	nffcount = 0
	do iff = 1,nff  
	  if (FFSensFlag(iff) .eq. 1) then 
		nffcount = nffcount + 1
		RespRef(nffcount) = FFSpTAv(1,iff,3) 
		names(nffcount) = FFName(iff,1)
	  end if  
	end do
	  names(nffcount+1) = 'range_sm0008'
	  names(nffcount+2) = 'range_sm0090'
	  RespRef(nffcount+1) = maxsm0008(1)-minsm0008(1)
	  RespRef(nffcount+2) = maxsm0090(1)-minsm0090(1)
	write(*,"(/,'      Flux      RefValue')")
	write(*,"(a10,2x,e12.4)") (trim(names(iff)),RespRef(iff),iff=1,nffsens)
  
	! Open output files (in OutputSens directory)
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'RawFluxes.csv'
	OPEN (unit=28, file=trim(GenFileName))
	write(28,"('Domain:,,',  a)") trim(DomName)
	write(28,"('Run ID:,,',  a)") trim(RunID)
	write(28,"('Comment:,,', a)") trim(Comment)//' '//'SENSITIVITY TEST 5'
	write(28,"(20(a10,','))") 'Fluxes,  ',names(:)
	write(28,"(a15,' ',20(e13.6,','))") 'Default,   , ',RespRef
	GenFileName = trim(DataOutRootDir) // trim(DomName) // '/Run' //    &
				  trim(RunID) // '/OutputSens/' // 'Sensitivity.csv'
	OPEN (unit=29, file=trim(GenFileName))
	write(29,"('Domain:,,',  a)") trim(DomName)
	write(29,"('Run ID:,,',  a)") trim(RunID)
	write(29,"('Comment:,,', a)") trim(Comment)
	write(29,"(20(a10,','))") 'Fluxes ',names(:)

	! Loop over parameters varying one at a time
	do ip = 1, nuu
	  nuucount = 0
	  if (UUSensFlag(ip) .eq. 1) then
		write(*,"('Vary parameter ',a16)") trim(UUName(ip,1))
		nuucount = nuucount + 1
		UU = UUinitCtl          ! Initialise quantities with global scope (XX0, UU)
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1
		UU(ip) = UUminctl(ip)   ! set UU to minimum value
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  if (it.eq.1) then
			minsm0008 = Zsm0008
			maxsm0008 = Zsm0008
			minsm0090 = Zsm0090
			maxsm0090 = Zsm0090
		  else
			minsm0008 = min(Zsm0008,minsm0008)
			maxsm0008 = max(Zsm0008,maxsm0008)
			minsm0090 = min(Zsm0090,minsm0090)
			maxsm0090 = max(Zsm0090,maxsm0090)
		  endif
		  XX0 = XX1
		end do
		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			RespUmin(nffcount) = FFSpTAv(1,iff,3)        
		  end if  
		end do
			RespUmin(nffcount+1) = maxsm0008(1) - minsm0008(1)
			RespUmin(nffcount+2) = maxsm0090(1) - minsm0090(1)

		UU = UUinitCtl          ! Initialise quantities with global scope (XX0, UU)
		CALL InitStores (XX0)
		InitFlag = 0            ! Initialise WaterDyn on first step by calling InitModel
		EndFlag = 1
		UU(ip) = UUmaxctl(ip)   ! set UU to minimum value
		do it=0,nSteps
		  Time0 = real(it) * 1.0
		  if (it.eq.nSteps) EndFlag = 0
		  CALL WaterDynStep (Time0,Time0, UU,VV, XX0, XX1, ZZP1,ZZPh1, ZZC1, AAP1,AAPh1, AAC1, TTime)
		  if (it.eq.1) then
			minsm0008 = Zsm0008
			maxsm0008 = Zsm0008
			minsm0090 = Zsm0090
			maxsm0090 = Zsm0090
		  else
			minsm0008 = min(Zsm0008,minsm0008)
			maxsm0008 = max(Zsm0008,maxsm0008)
			minsm0090 = min(Zsm0090,minsm0090)
			maxsm0090 = max(Zsm0090,maxsm0090)
		  endif
		  XX0 = XX1
		end do
		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1
			RespUmax(nffcount) = FFSpTAv(1,iff,3)        
		  end if  
		end do
			RespUmax(nffcount+1) = maxsm0008(1) - minsm0008(1)
			RespUmax(nffcount+2) = maxsm0090(1) - minsm0090(1)

		nffcount = 0
		do iff = 1,nff  
		  if (FFSensFlag(iff) .eq. 1) then 
			nffcount = nffcount + 1  
			Sens(nffcount,nuucount) = (RespUmax(nffcount)-RespUmin(nffcount))/RespRef(nffcount)
		  end if  
		end do
		Sens(nffcount+1,nuucount) = (RespUmax(nffcount+1)-RespUmin(nffcount+1))/RespRef(nffcount+1)
		Sens(nffcount+2,nuucount) = (RespUmax(nffcount+2)-RespUmin(nffcount+2))/RespRef(nffcount+2)

		CALL CPU_TIME (cpuTimeNow)
		write(*,"('Done parameter ',a16,'      CPUtime (sec) = ',f9.2)")    &
		  trim(UUName(ip,1)), cpuTimeNow - cpuTimeBegin
	!    write(6,*) ' Resp =',Resp(:)
	!    write(6,*) ' Sens =',Sens(:,nuucount)
	!    write(6,*) ' '
		write(28,"(a10,', ',20(e13.6,','))") UUName(ip,1),UURef(ip),Resp(:)    
		write(29,"(a10,', ',20(e13.6,','))") UUName(ip,1),Sens(:,nuucount)
	  end if   
	end do 

endif

close(28)
close(29)
!close(7)
close(99)
close(66)



END PROGRAM WaterDynM


