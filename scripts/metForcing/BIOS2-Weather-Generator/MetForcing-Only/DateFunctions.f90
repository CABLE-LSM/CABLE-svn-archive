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
! Define interfaces for the +, -, .eq., .ne., .lt., .gt., .le., .ge.
! operators to:
! * add, subtract a number of days (integer) from a dmydate, with
! AddDay, SubDay; 
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
! Extends the '.eq.' operator to enable comparison of two dates for
! equality.
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
! Extends the '.eq.' operator to enable comparison of two dates for
! equality.
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
! Extends the '.ne.' operator to enable comparison of two dates for
! inequality.
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
! Extends the '.lt.' operator to enable comparison of two dates for
! Date1 < Date2.
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
! Extends the '.gt.' operator to enable comparison of two dates for
! Date1 > Date2.
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
! Extends the '.le.' operator to enable comparison of two dates Date1 <=
! Date2.
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
! Extends the '.ge.' operator to enable comparison of two dates Date1 >=
! Date2.
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
! Extends the '.eq.' operator to enable comparison of two time-series
! months for equality.
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
! Extends the '.eq.' operator to enable comparison of two months of year
! for equality.
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
  JulianDay = D + int((153.0_8 * M - 457.0_8)/5.0_8) + (365.0_8*Y) + &
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
