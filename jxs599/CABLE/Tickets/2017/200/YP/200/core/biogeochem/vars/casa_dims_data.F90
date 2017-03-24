MODULE casadimension

   USE cable_def_types_mod, ONLY : mp, r_2, mvtype, ms

   IMPLICIT NONE



  INTEGER, PARAMETER :: mdyear=365         ! days per year
  INTEGER, PARAMETER :: mdmonth=30         ! days per month
  INTEGER, PARAMETER :: mdweek=7           ! days per week
  INTEGER, PARAMETER :: mmyear=12          ! month per year
  INTEGER, PARAMETER :: mt=36500           ! integration time step
  INTEGER, PARAMETER :: mpftmax=2          ! max. PFT/cell
  INTEGER, PARAMETER :: mplant = 3         ! plant pools
  INTEGER, PARAMETER :: mlitter= 3         ! litter pools
  INTEGER, PARAMETER :: msoil  = 3         ! soil pools
  INTEGER, PARAMETER :: mso    = 12        ! soil order number
  INTEGER, PARAMETER :: mhwp  = 1         ! harvested wood pools
  INTEGER, PARAMETER :: mclear  = 1         ! forest clearing pools
! BP put icycle into namelist file
  INTEGER            :: icycle
!  INTEGER, PARAMETER :: icycle=3           ! =1 for C, =2 for C+N; =3 for C+N+P
  INTEGER, PARAMETER :: mstart=1           ! starting time step
  INTEGER, PARAMETER :: mphase=4           ! phen. phases
  REAL(r_2),    PARAMETER :: deltcasa=1.0/365.0 ! year
  REAL(r_2),    PARAMETER :: deltpool=1.0       ! pool delt(1day)

END MODULE casadimension


