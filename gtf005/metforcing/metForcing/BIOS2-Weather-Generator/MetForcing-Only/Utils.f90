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
