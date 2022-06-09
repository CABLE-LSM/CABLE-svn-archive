MODULE cable_latent_heat_module

  IMPLICIT NONE

  PUBLIC latent_heat_flux 
  PRIVATE

CONTAINS

SUBROUTINE Latent_heat_flux( mp, CTFRZ, dels, soil_zse, soil_swilt,           &
                             cable_user_l_new_reduce_soilevp, pwet, air_rlam,  &
                             ssnow_snowd, ssnow_wb, ssnow_wbice,             &
                             ssnow_pudsto, ssnow_pudsmx, ssnow_potev,          &
                             ssnow_wetfac, ssnow_evapfbl, ssnow_cls,          & 
                             ssnow_tss, canopy_fes, canopy_fess, canopy_fesp  )

USE cable_def_types_mod, ONLY : r_2
IMPLICIT NONE

INTEGER :: mp
REAL(KIND=r_2), INTENT(OUT) :: canopy_fes(mp)
REAL(KIND=r_2), INTENT(OUT) :: canopy_fess(mp)
REAL(KIND=r_2), INTENT(OUT) :: canopy_fesp(mp)
REAL, INTENT(OUT) :: ssnow_cls(mp)
REAL, INTENT(IN OUT) :: pwet(mp)
REAL, INTENT(IN OUT) :: ssnow_wetfac(mp)


REAL, INTENT(IN) :: CTFRZ
REAL, INTENT(IN) :: dels
REAL, INTENT(IN) :: soil_zse
REAL, INTENT(IN) :: soil_swilt(mp)
LOGICAL , INTENT(IN) :: cable_user_l_new_reduce_soilevp

REAL, INTENT(IN) :: air_rlam(mp)
REAL, INTENT(IN) :: ssnow_snowd(mp)
REAL, INTENT(IN) :: ssnow_pudsto(mp)
REAL, INTENT(IN) :: ssnow_pudsmx(mp)
REAL, INTENT(IN) :: ssnow_potev(mp)
REAL, INTENT(IN) :: ssnow_evapfbl(mp)
REAL(KIND=r_2), INTENT(IN) :: ssnow_wb(mp)
REAL(KIND=r_2), INTENT(IN) :: ssnow_wbice(mp)
REAL, INTENT(IN) :: ssnow_tss(mp)

REAL, DIMENSION(mp) ::                                                      &
  frescale,  flower_limit, fupper_limit

INTEGER :: j

   ! Soil latent heat:
   canopy_fess= ssnow_wetfac * ssnow_potev
   WHERE (ssnow_potev < 0. ) canopy_fess = ssnow_potev
   
   ! Reduce soil evap due to presence of puddle
   pwet = max(0.,min(0.2,ssnow_pudsto/max(1.,ssnow_pudsmx)))
   canopy_fess = canopy_fess * (1.-pwet)

   frescale = soil_zse * 1000. * air_rlam / dels         

   DO j=1,mp
      
      IF(ssnow_snowd(j) < 0.1 .AND. canopy_fess(j) .GT. 0. ) THEN

        IF (.not.cable_user_l_new_reduce_soilevp) THEN
         flower_limit(j) = REAL(ssnow_wb(j))-soil_swilt(j)/2.0
        ELSE
         ! E.Kowalczyk 2014 - reduces the soil evaporation
         flower_limit(j) = REAL(ssnow_wb(j))-soil_swilt(j)
        ENDIF
         fupper_limit(j) = MAX( 0._r_2,                                        &
                           flower_limit(j) * frescale(j)                       &
                           - ssnow_evapfbl(j)*air_rlam(j)/dels)

         canopy_fess(j) = MIN(canopy_fess(j), fupper_limit(j))
         
         fupper_limit(j) = REAL(ssnow_wb(j)-ssnow_wbice(j)) * frescale(j)

         canopy_fess(j) = min(canopy_fess(j), fupper_limit(j))

      END IF

      ssnow_cls(j)=1.
      
      IF (ssnow_snowd(j) >= 0.1 .and. ssnow_potev(j) > 0.) THEN

         ssnow_cls(j) = 1.1335
         canopy_fess(j) = MIN( (ssnow_wetfac(j)*ssnow_potev(j))*ssnow_cls(j), &
                          ssnow_snowd(j)/dels*air_rlam(j)*ssnow_cls(j))
      
      ENDIF

   ENDDO 
   
   ! Evaporation form soil puddle
   canopy_fesp = min(ssnow_pudsto/dels*air_rlam,max(pwet*ssnow_potev,0.))
   canopy_fes = canopy_fess + canopy_fesp

RETURN
END SUBROUTINE latent_heat_flux

END MODULE cable_latent_heat_module
