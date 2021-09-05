MODULE landuse_constant
  USE cable_def_types_mod,    ONLY: mvtype,mstype,mland,ms,msn,ncp,ncs
  USE casadimension,     ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  IMPLICIT NONE
!  integer, parameter                     :: sp =selected_real_kind(8)
!  integer, parameter                     :: dp =selected_real_kind(16)
  !  
  INTEGER,   PARAMETER                   :: mstate   = 12  ! number of land use states
  INTEGER,   PARAMETER                   :: mvmax    = mvtype 
  INTEGER,   PARAMETER                   :: mharvw   = 5           

  real*8,    PARAMETER                   :: thresh_frac  = 1.0e-6
  REAL*8,    PARAMETER, DIMENSION(mwood) :: fwoodprod    =(/0.3,0.4,0.4/)
  REAL*8,    PARAMETER, DIMENSION(mwood) :: ratewoodprod =(/1.0,0.1,0.01/)
  REAL*8,    PARAMETER, DIMENSION(mwood) :: fracwoodseed =(/0.4,0.3,0.3/)
  REAL*8,    PARAMETER, DIMENSION(mwood) :: fracgrassseed=(/0.6,0.0,0.4/)
  REAL*8,    PARAMETER                   :: fseedling = 0.001
END MODULE landuse_constant
