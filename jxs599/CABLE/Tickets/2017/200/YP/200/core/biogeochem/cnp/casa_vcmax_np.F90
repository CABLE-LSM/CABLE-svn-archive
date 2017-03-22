MODULE casa_cnp_module
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

REAL FUNCTION vcmax_np(nleaf, pleaf)
IMPLICIT NONE
REAL, INTENT(IN) :: nleaf ! leaf N in g N m-2 leaf
REAL, INTENT(IN) :: pleaf ! leaf P in g P m-2 leaf

!Walker, A. P. et al.: The relationship of leaf photosynthetic traits – Vcmax and Jmax – 
!to leaf nitrogen, leaf phosphorus, and specific leaf area: 
!a meta-analysis and modeling study, Ecology and Evolution, 4, 3218-3235, 2014.
    vcmax_np = exp(3.946 + 0.921*log(nleaf) + 0.121*log(pleaf) + &
         0.282*log(pleaf)*log(nleaf)) * 1.0e-6 ! units of mol m-2 (leaf)


END FUNCTION vcmax_np

END MODULE casa_cnp_module
