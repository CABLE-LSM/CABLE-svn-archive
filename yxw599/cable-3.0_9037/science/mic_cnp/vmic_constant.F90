module mic_constant
  USE cable_def_types_mod, ONLY : r_2
  IMPLICIT NONE
!  integer,  parameter  :: r_2 = SELECTED_REAL_KIND(12, 60)
  integer,  parameter   :: kinetics=1
  integer,  parameter   :: nyeqpool=500
  integer,  parameter   :: diag=0       ! =1 for printout 0 no prinout
  integer,  parameter   :: outp=1       ! output site
!  integer,  parameter  :: msite=182    ! number of sites
!  integer,  parameter  :: mp=1000      ! number of site the model runs for
!  integer,  parameter  :: ms=15        ! soil layers
  integer,  parameter  :: mcpool=7     ! number of C pools
!  integer,  parameter  :: nfvar=36     ! number of data input variables
  real,     parameter  :: delt= 1.0    ! one hour
!  integer,  parameter  :: mpft=17      ! number of PFTs
  real(r_2),PARAMETER  :: diffsoc  =(1.0/24.0)* 2.74e-3  !cm2/hour   
                                        ! m2/hour  ! see Table 1,Camino-Serrano et al. (2018)
  
end module mic_constant