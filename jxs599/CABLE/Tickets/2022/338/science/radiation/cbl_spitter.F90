MODULE cbl_spitter_module

IMPLICIT NONE

PUBLIC :: spitter
PRIVATE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calc. beam fraction: Spitters et al. 1986, agric. for meteorol., 38:217-229
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION spitter(mp, cpi, doy, coszen, fsd) RESULT(fbeam)
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp         ! CABLE vector length of active tiles
REAL :: fbeam(mp)                 ! beam fraction of total SW (f^n result)
REAL, INTENT(IN) :: cpi           ! PI
REAL, INTENT(IN) :: coszen(mp)    ! cos(zenith angle of sun)
REAL, INTENT(IN) :: fsd(mp)       ! short wave down (positive) w/m^2
INTEGER, INTENT(IN) :: doy(mp)    ! day of year !typecast from integer
!local
REAL :: tmpr(mp)
REAL :: tmpk(mp)
REAL :: tmprat(mp)
REAL, PARAMETER :: solcon = 1370.0

fbeam = 0.0
tmpr = 0.847 + coszen * (1.04 * coszen - 1.61)
tmpk = (1.47 - tmpr) / 1.66

WHERE (coszen > 1.0e-10 .AND. fsd > 10.0)
  tmprat = fsd / ( solcon * ( 1.0 + 0.033 * COS( 2.0 * cpi* ( doy-10.0 )       &
           / 365.0 ) ) * coszen )
ELSE WHERE
  tmprat = 0.0
END WHERE

WHERE ( tmprat > 0.22 )
  fbeam = 6.4 * ( tmprat - 0.22 )**2
END WHERE
WHERE ( tmprat > 0.35 )
  fbeam = MIN( 1.66 * tmprat - 0.4728, 1.0 )
END WHERE
WHERE ( tmprat > tmpk )
  fbeam = MAX( 1.0 - tmpr, 0.0 )
END WHERE

RETURN
END FUNCTION spitter

END MODULE cbl_spitter_module
