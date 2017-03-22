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

SUBROUTINE casa_pdummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Pplant(:,:) = casapool%Nplant(:,:) / casapool%ratioNPplant(:,:)

END SUBROUTINE casa_pdummy

END MODULE casa__mod
