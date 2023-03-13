MODULE vmic_constant_mod
  USE cable_def_types_mod, ONLY : r_2
  IMPLICIT NONE
  integer,  parameter   :: kinetics=1
  integer,  parameter   :: nyeqpool=500
  integer,  parameter   :: diag=0           ! =1 for printout 0 no prinout
  integer,  parameter   :: outp=14          ! output site
  integer,  parameter   :: mcpool=7         ! number of C pools
  real(r_2),parameter   :: deltvmic= 1.0    ! one hour

  real(r_2),PARAMETER   :: diffsoc  =(1.0/24.0)* 2.74e-3  !cm2/hour   
                                        ! m2/hour  ! see Table 1,Camino-Serrano et al. (2018)
  integer               :: vmicrobe         ! =0 for not used, 1 (C), 2(CN) and 3(CNP)
  
END MODULE vmic_constant_mod
