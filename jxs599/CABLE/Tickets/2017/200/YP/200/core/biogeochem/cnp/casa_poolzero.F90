MODULE casa__mod

USE cable_def_types_mod
USE casadimension
USE casaparm
USE casavariable
USE phenvariable
USE cable_common_module, only: cable_user ! Custom soil respiration: Ticket #42

IMPLICIT NONE
  REAL(r_2), PARAMETER :: zero = 0.0_r_2
  REAL(r_2), PARAMETER :: one  = 1.0_r_2

CONTAINS

SUBROUTINE casa_poolzero(n,ipool,casapool)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ipool
  TYPE (casa_pool), INTENT(INOUT) :: casapool

  WRITE(57,*) ' WARNING: negative pools are reset to ZERO!!'
  SELECT CASE(ipool)
  CASE(1)
     WRITE(57,*) 'plant carbon pool size negative!!'
     WRITE(57,*) 'plant C pools: ', n,casapool%cplant(n,:)
  CASE(2)
     WRITE(57,*) 'plant nitrogen pool size negative!!'
     WRITE(57,*) 'plant C pools: ',n,casapool%cplant(n,:)
     WRITE(57,*) 'plant N pools: ',n,casapool%nplant(n,:)
  CASE(3)
     WRITE(57,*) 'litter carbon pool size negative!!'
     WRITE(57,*) 'litter C pools: ',n,casapool%clitter(n,:)
  CASE(4)
     WRITE(57,*) 'litter nitrogen pool size negative!!'
     WRITE(57,*) 'carbon pool: ',n,casapool%clitter(n,:)
     WRITE(57,*) 'nitrogen pools: ',n,casapool%nlitter(n,:)
  CASE(5)
     WRITE(57,*) 'soil carbon pool size negative!!'
     WRITE(57,*) 'soil C pools: ',n,casapool%csoil(n,:)
  CASE(6)
     WRITE(57,*) 'soil nitrogen pool size negative!!'
     WRITE(57,*) 'soil C pools: ', n,casapool%csoil(n,:)
     WRITE(57,*) 'soil N pools: ', n,casapool%nsoil(n,:)
  END SELECT

END SUBROUTINE casa_poolzero



END MODULE casa__mod
