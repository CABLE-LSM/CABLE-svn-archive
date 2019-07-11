#if !defined(UM_JULES)
!-----------------------------------------------------------------------
! This subroutine describes a simple uptake
! by the oceans of atmospheric CO2
! The work is based upon the paper by Joos.
!
! The parameters chosen replicate those of the 3-D model, see
! Table 2 of Joos et al.
!-----------------------------------------------------------------------
!
! Huntingford and Cox (March 1999)

  SUBROUTINE OCEAN_CO2(                                           &
    IYEAR,YEAR_CO2,CO2_ATMOS_PPMV,CO2_ATMOS_INIT_PPMV,DT_OCEAN,   &
    FA_OCEAN,OCEAN_AREA,GRAD_CO2_ATMOS_PPMV,YEAR_RUN,             &
    T_OCEAN_INIT,NFARRAY,D_OCEAN_ATMOS                            &
  )

    USE jules_print_mgr, ONLY :                                     &
      jules_message,                                               &
      jules_print
    IMPLICIT NONE

    INTEGER ::                                                    &
      I,                                                          &
             !WORK Looping counter
      J,                                                          &
             !WORK Looping counter
      ISTEP,                                                      &
             !IN Looping counter from main subroutine that shows
             !   calls to Joos model
      IYEAR,                                                      &
             !IN Years since start of run
      YEAR_CO2,                                                   &
             !IN Years between updated atmospheric CO2
             !   concentration (yr)
      NCALLYR,                                                    &
             !IN Number of calls per year the ocean routine
             !   (with CO2 conc. fixed)
      YEAR_RUN,                                                   &
             !IN Number of years in simulation (yr)
      NFARRAY
             !IN Array size for FA_OCEAN

    PARAMETER(NCALLYR=20)

    REAL ::                                                       &
      DT_OCEAN,                                                   &
             !IN Mixed-layer temperature anomaly (K)
      CO2_ATMOS_PPMV,                                             &
             !WORK Atmospheric CO2 concentation during period
             !     YEAR_CO2 (ppm)
      CO2_ATMOS_SHORT,                                            &
             !WORK Atmospheric CO2 concentation on the shortime
             !     (1/NCALLYR) timescale (ppm)
      CO2_ATMOS_INIT_PPMV,                                        &
             !WORK Initial atmospheric concentration
      DCO2_ATMOS,                                                 &
             !WORK Perturbation in atmospheric concentration (ppm)
      TIMESTEP_CO2,                                               &
             !WORK Timestep (yr)
      RS(YEAR_RUN*NCALLYR),                                       &
             !WORK Response function values
      GRAD_CO2_ATMOS_PPMV
             !IN Gradient of atmospheric CO2 (ppm/y) during
             !   YEAR_CO2 period

    REAL ::                                                       &
      FA_OCEAN(NFARRAY),                                          &
             !IN/OUT CO2 fluxes from the atmosphere to the ocean
             !       (ie positive downwards) (ppm/m2/yr)
      DCO2_OCEAN,                                                 &
             !IN/OUT Carbon dioxide concentration perturbation
             !       in mixed-layer (ppm)
      CO2_OCEAN,                                                  &
             !OUT Carbon dioxide concentration in mixed-layer (ppm)
      CO2_OCEAN_INIT,                                             &
             !OUT Initial ocean carbon dioxide concentration (ppm)
      DCO2_OCEAN_MOL,                                             &
             !WORK Carbon dioxide concentration perturbation
             !     of mixed-layer (umol/kg)
      T_MIXED_INIT,                                               &
             !WORK Initial global mean ocean surface temperature (C)
      T_OCEAN_INIT,                                               &
             !IN Initial global surface temperature of the ocean (K)
      H,                                                          &
             !WORK Mixed-layer depth (m)
      K_G,                                                        &
             !WORK Gas exchange coefficient (/m2/yr)
      C,                                                          &
             !WORK Units conversion parameter (umol m3/ppm/
      OCEAN_AREA
             !WORK Ocean area (m2)

    REAL ::                                                       &
      D_OCEAN_ATMOS
             ! OUT Change in atmospheric CO2 concentration a
             !     result of ocean-atmosphere feedbacks between
             !     calls to SCENARIO (ppm/"YEAR_CO2")

    DATA K_G/0.1306/
    DATA H/40.0/
    DATA C/1.722E17/



    IF(YEAR_RUN*NCALLYR > NFARRAY) THEN
      WRITE(jules_message,*) 'Array size too small for FA_OCEAN'
      CALL jules_print('ocean_co2',jules_message)
    ENDIF

    CALL RESPONSE(NCALLYR,YEAR_RUN,RS)

    D_OCEAN_ATMOS = 0.0

!-----------------------------------------------------------------------
! Assume initial mixed-layer temperature identical to diagnosed global
! surface temperature.
!-----------------------------------------------------------------------
    T_MIXED_INIT= T_OCEAN_INIT - 273.15

    TIMESTEP_CO2 = 1.0 / FLOAT(NCALLYR)

! Loops internally within the model to cover iterations within
! period YEAR_CO2
    DO J = 1,NCALLYR*YEAR_CO2

      ISTEP = ((IYEAR-YEAR_CO2)*NCALLYR)+J
!                             !ISTEP is the integer that indexes calcula
!                             !using the Joos model. Recall this subrout
!                             !is called at end of period YEAR_CO2

      CO2_OCEAN_INIT = CO2_ATMOS_INIT_PPMV           !Assume starting
              !from equilibrium.

!-----------------------------------------------------------------------
! Introduce linear correction in atmospheric CO2 concentration
! down to timescale 1/NCALLYR
!-----------------------------------------------------------------------
      CO2_ATMOS_SHORT = CO2_ATMOS_PPMV + GRAD_CO2_ATMOS_PPMV *    &
                   ((FLOAT(J)/FLOAT(NCALLYR))-0.5*FLOAT(YEAR_CO2))

      DCO2_ATMOS = CO2_ATMOS_SHORT-CO2_ATMOS_INIT_PPMV

!-----------------------------------------------------------------------
! Calculate perturbation in dissolved inorganic carbon (Eqn (3) of Joos)
! for timestep indexed by ISTEP.
!-----------------------------------------------------------------------
      DCO2_OCEAN_MOL = 0.0  ! Initialised for loop

      IF(ISTEP >= 2) THEN
        DO I = 1,ISTEP-1
          DCO2_OCEAN_MOL = DCO2_OCEAN_MOL +                       &
                        (C/H)*FA_OCEAN(I)*RS(ISTEP-I)*TIMESTEP_CO2
        ENDDO
      ENDIF

!----------------------------------------------------------------------
! Relation between DCO2_OCEAN_MOL and DCO2_OCEAN:   Joos et al.
! ie Convert from umol/kg to ppm
!----------------------------------------------------------------------
      DCO2_OCEAN =                                                &
         (1.5568-(1.3993*1.0E-2*T_MIXED_INIT))*DCO2_OCEAN_MOL     &
      + (7.4706-(0.20207*T_MIXED_INIT))*1.0E-3*(DCO2_OCEAN_MOL**2)&
      - (1.2748-(0.12015*T_MIXED_INIT))*1.0E-5*(DCO2_OCEAN_MOL**3)&
      + (2.4491-(0.12639*T_MIXED_INIT))*1.0E-7*(DCO2_OCEAN_MOL**4)&
      - (1.5468-(0.15326*T_MIXED_INIT))*1.0E-10*(DCO2_OCEAN_MOL**5)

!----------------------------------------------------------------------
! Now incorporate correction suggested by Joos
!----------------------------------------------------------------------
      CO2_OCEAN = CO2_OCEAN_INIT + DCO2_OCEAN
      CO2_OCEAN = CO2_OCEAN * EXP(0.0423*DT_OCEAN)

      FA_OCEAN(ISTEP) = (K_G/OCEAN_AREA)                          &
                      * (CO2_ATMOS_SHORT-CO2_OCEAN)

!----------------------------------------------------------------------
! Now calculate D_OCEAN_ATMOS (D_OCEAN_ATMOS is positive when
! flux is upwards)
!----------------------------------------------------------------------
      D_OCEAN_ATMOS = D_OCEAN_ATMOS                               &
                    - FA_OCEAN(ISTEP)*OCEAN_AREA*TIMESTEP_CO2
    ENDDO

    RETURN

  END SUBROUTINE OCEAN_CO2
#endif
