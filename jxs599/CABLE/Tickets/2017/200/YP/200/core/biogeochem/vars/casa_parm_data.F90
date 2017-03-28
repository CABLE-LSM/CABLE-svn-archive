MODULE casaparm
  USE casadimension

  IMPLICIT NONE
  INTEGER, PARAMETER :: initcasa= 1   ! =0 spin; 1 restart file
  INTEGER, PARAMETER :: iceland  = 17 !=13 for casa vegtype =15 for IGBP vegtype
  INTEGER, PARAMETER :: cropland = 9  ! 12 and 14 for IGBP vegtype
  INTEGER, PARAMETER :: croplnd2 =10  ! ditto
  INTEGER, PARAMETER :: forest  = 3
  INTEGER, PARAMETER :: shrub   = 2
  INTEGER, PARAMETER :: grass   = 1
  INTEGER, PARAMETER :: icewater= 0
  INTEGER, PARAMETER :: LEAF    = 1
  INTEGER, PARAMETER :: WOOD    = 2
  INTEGER, PARAMETER :: FROOT   = 3
!  INTEGER, PARAMETER :: LABILE  = 4
  INTEGER, PARAMETER :: METB    = 1
  INTEGER, PARAMETER :: STR     = 2
  INTEGER, PARAMETER :: CWD     = 3
  INTEGER, PARAMETER :: MIC     = 1
  INTEGER, PARAMETER :: SLOW    = 2
  INTEGER, PARAMETER :: PASS    = 3
  INTEGER, PARAMETER :: PLAB    = 1
  INTEGER, PARAMETER :: PSORB   = 2
  INTEGER, PARAMETER :: POCC    = 3
 !! vh_js !! LALLOC moved to bgcdriver to allow for value to be switchable
  !INTEGER, PARAMETER :: LALLOC  = 0      !=0 constant; 1 variable
  REAL(r_2), PARAMETER :: z30=0.3
  REAL(r_2), PARAMETER :: R0=0.3
  REAL(r_2), PARAMETER :: S0=0.3
  REAL(r_2), PARAMETER :: fixed_stem=1.0/3.0
  REAL(r_2), PARAMETER :: Q10alloc=2.0
  REAL(r_2), PARAMETER :: ratioNCstrfix = 1.0/150.0
  REAL(r_2), PARAMETER :: ratioNPstrfix = 25.0
  REAL(r_2), PARAMETER :: fracCbiomass = 0.50
  REAL(r_2), PARAMETER :: tsoilrefc=25.0
  REAL(r_2), PARAMETER :: tkzeroc=273.15
  REAL(r_2), PARAMETER :: frootparma = 0.3192
  REAL(r_2), PARAMETER :: frootparmb =-0.0485
  REAL(r_2), PARAMETER :: frootparmc = 0.1755
  REAL(r_2), PARAMETER :: xweightalloc = 0.2
END MODULE casaparm


