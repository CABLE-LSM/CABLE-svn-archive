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

SUBROUTINE avgsoil(veg,soil,casamet)
! Get avg soil moisture, avg soil temperature
! need to estimate the land cell mean soil temperature and moisture weighted by the area fraction
! of each tile within the land cell

  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  INTEGER                     :: ns,nland

  casamet%tsoilavg   = 0.0
  casamet%moistavg   = 0.0
  casamet%btran      = 0.0

  DO ns = 1, ms
  DO nland=1,mp
    casamet%tsoilavg(nland)  = casamet%tsoilavg(nland)+veg%froot(nland,ns)  &
                             * casamet%tsoil(nland,ns)
    casamet%moistavg(nland)  = casamet%moistavg(nland)+ veg%froot(nland,ns) &
                           * min(soil%sfc(nland),casamet%moist(nland,ns))
    casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
            * (min(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland)) &
            /(soil%sfc(nland)-soil%swilt(nland))

 ! Ticket#121

    casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
            * (max(min(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland),0.0)) &
            /(soil%sfc(nland)-soil%swilt(nland))

  ENDDO
  ENDDO

END SUBROUTINE avgsoil


END MODULE casa__mod
