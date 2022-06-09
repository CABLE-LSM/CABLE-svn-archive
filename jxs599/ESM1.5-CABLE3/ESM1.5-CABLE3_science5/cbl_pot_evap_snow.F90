MODULE cbl_pot_evap_snow_module

IMPLICIT NONE

PUBLIC Penman_Monteith
PUBLIC Humidity_deficit_method
PRIVATE

CONTAINS

FUNCTION Penman_Monteith( mp, Ctfrz, CRMH2o, Crmair, CTETENA, CTETENB,         &
                          CTETENC, air_dsatdk, air_psyc, air_rho, air_rlam,             & 
                          met_tvair, met_pmb, met_qvair,                       &
                          ground_H_flux, canopy_fns, ssnow_rtsoil, ssnow_isflag ) RESULT(ssnowpotev)
USE cbl_qsat_module, ONLY : qsatfjh 
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
!REAL, INTENT(OUT)   :: ssnowpotev(mp)     ! returned result of function
REAL   :: ssnowpotev(mp)     ! returned result of function

REAL, INTENT(IN) :: Ctfrz 
REAl :: CRMH2o
REAl :: Crmair
REAl :: CTETENA, CTETENB, CTETENC
REAL, INTENT(IN) :: air_dsatdk(mp)
REAL, INTENT(IN) :: air_psyc(mp)
REAL, INTENT(IN) :: air_rho(mp)
REAL, INTENT(IN) :: air_rlam(mp)
REAL, INTENT(IN) :: met_tvair(mp) 
REAL, INTENT(IN) :: met_qvair(mp) 
REAL, INTENT(IN) :: met_pmb(mp)
REAL, INTENT(IN) :: ground_H_flux(mp)
REAL, INTENT(IN) :: canopy_fns(mp)
REAL, INTENT(IN) :: ssnow_rtsoil(mp)
INTEGER, INTENT(IN) :: ssnow_isflag(mp)

!local vars
REAL :: sss(mp)          ! var for Penman-Monteith soil evap  
REAL :: cc1(mp)          ! var for Penman-Monteith soil evap                                                    &
REAL :: cc2(mp)          ! var for Penman-Monteith soil evap
REAL :: qsatfvar(mp)
INTEGER :: j

! Penman-Monteith formula
sss=air_dsatdk
cc1=sss/(sss+air_psyc )
cc2=air_psyc /(sss+air_psyc )

CALL qsatfjh( mp, qsatfvar, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC,         &
              met_tvair-CTfrz,met_pmb)
!ESM15!  CALL qsatfjh(qsatfvar,met%tvair-Ctfrz,met%pmb)

!INH 10-1-2017 - this P-M implementation is incorrect over snow.
!variable ssnowpotev is actually the latent heat flux associated with
!potential evaporation.
!Needs to be addressed/simplified at a later date - involves changes
!to HDM method and latent_heat_flux() and elsewhere

ssnowpotev = cc1 * (canopy_fns - ground_H_flux) + &
             cc2 * air_rho * air_rlam*(qsatfvar  - met_qvair)/ssnow_rtsoil

END FUNCTION Penman_Monteith


! ------------------------------------------------------------------------------
! method alternative to P-M formula above
FUNCTION Humidity_deficit_method( mp, Ctfrz,  &
                                 air_rho,air_rlam,           & 
                                 dq,& 
                                 ssnow_rtsoil,    &
                                 ssnow_snowd, ssnow_tgg &
                                 ) RESULT(ssnowpotev)

integer :: mp
REAL ::  ssnowpotev(mp)

REAL:: Ctfrz 
REAL :: air_rho(mp)    !
REAL :: air_rlam(mp)    !
REAL :: dq(mp)       ! sat spec hum diff.
REAL :: ssnow_snowd(mp)    ! 
REAL :: ssnow_tgg(mp)    ! 
REAL :: ssnow_rtsoil(mp)    !

!local vars
INTEGER :: j

DO j=1,mp

  IF( ssnow_snowd(j)>1.0 .OR. ssnow_tgg(j) .EQ. Ctfrz ) THEN
    dq(j) = MAX( -0.1e-3, dq(j))
  END IF

ENDDO

ssnowpotev = air_rho * air_rlam * dq / ssnow_rtsoil
   
END FUNCTION Humidity_deficit_method

END MODULE cbl_pot_evap_snow_module
