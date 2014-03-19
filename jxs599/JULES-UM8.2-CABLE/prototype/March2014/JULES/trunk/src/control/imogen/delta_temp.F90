#if !defined(UM_JULES)
!-----------------------------------------------------------------------
! This subroutine calls a simple two-box thermal model in order to
! estimate the global mean land surface temperature.
! The code uses an implicit solver.
! C. Huntingford, 11/09/98.
!-----------------------------------------------------------------------
  SUBROUTINE DELTA_TEMP(                                          &
    N_OLEVS,F_OCEAN,KAPPA,LAMBDA_L,LAMBDA_O,MU,Q,DTEMP_L,DTEMP_O  &
  )

    USE jules_print_mgr, ONLY :                                     &
      jules_message,                                               &
      jules_print
    IMPLICIT NONE

    INTEGER ::                                                    &
      N_OLEVS  !IN Number of ocean thermal layers

    REAL ::                                                       &
      Q,                                                          &
               !IN Increase in global radiative forcing (W/m2)
      LAMBDA_L,                                                   &
               !IN Inverse climate sensitivity over land (W/K/m2)
      LAMBDA_O,                                                   &
               !IN Inverse climate sensitivity over ocean (W/K/m2)
      MU,                                                         &
               !IN Ratio of land-to-ocean temperature anomalies (.)
      F_OCEAN,                                                    &
               !IN Fractional coverage of planet with ocean
      KAPPA    !IN Ocean eddy diffusivity (W/m/K)

    INTEGER ::                                                    &
      ITER_PER_YEAR,                                              &
               !WORK Iterations per year
      TIMESTEPS,                                                  &
               !WORK Number of iterations per call to DELTA_TEMP
      I,J      !WORK Looping parameter
    PARAMETER(ITER_PER_YEAR = 20)

    REAL ::                                                       &
      DTEMP_L,                                                    &
               !OUT Land mean temperature anomaly
      DTEMP_O(N_OLEVS)
               !IN/OUT Land mean temperature anomaly

    REAL ::                                                       &
      RHOCP,                                                      &
               !WORK Rho times cp for the ocean (J/K/m3)
      FLUX_TOP,                                                   &
               !WORK Flux into the top of the ocean (W/m2)
      FLUX_BOTTOM
               !WORK Flux into the top of the ocean (W/m2)

    REAL ::                                                       &
      DTIME,                                                      &
               !WORK Timestep (s)
      DZ(1:N_OLEVS),                                              &
               !WORK Distance step (m)
      SEC_YEAR !WORK Seconds in a year (s)

    REAL ::                                                       &
      R1(1:N_OLEVS),                                              &
               !WORK Mesh ratio
      R2(1:N_OLEVS),                                              &
               !WORK Mesh ratio
      LAMBDA_OLD,                                                 &
               !WORK Inhomogenous component in implicit scheme (/s)
      LAMBDA_NEW,                                                 &
               !WORK Inhomogenous component in implicit scheme (/s)
      P        !WORK ``Dirchlet'' part of mixed boundary
               !     condition at ocean surface (/m)

    REAL ::                                                       &
      U_OLD(1:N_OLEVS),                                           &
               !WORK Old ocean temperature anomalies
      U_NEW(1:N_OLEVS)
               !WORK New ocean temperature anomalies

    REAL ::                                                       &
      FACTOR_DEP,                                                 &
               !WORK Factor by which layers are changed
      DEPTH,                                                      &
               !WORK Cumulated depth of layers (m)
      DZ_TOP   !WORK Depth of the top layer (m)

    PARAMETER(RHOCP=4.04E6)
    PARAMETER(SEC_YEAR = 3.1536E7)


!-----------------------------------------------------------------------
! Set up parameters for the model
!-----------------------------------------------------------------------

    DTIME = SEC_YEAR/FLOAT(ITER_PER_YEAR)
    TIMESTEPS = ITER_PER_YEAR


! Here the variable depths are calculated, along with the "R" values

! The mesh is as follows
!
!      -------------            U(1)
!                  (This is the top layer, called TEM(0) at points)
!          R(1)
!
!      ------------             U(2)
!
!          R(2)
!
!
!           |
!           |
!          \ /
!
!      ------------             U(N_OLEVS) = 0.0

    DZ_TOP = 2.0
    DEPTH = 0.0

    DO I = 1,N_OLEVS
      FACTOR_DEP = 2.71828**(DEPTH/500.0)
      DZ(I) = DZ_TOP*FACTOR_DEP
      DEPTH = DEPTH + DZ(I)
    ENDDO

! Now arrange that the lowest depth (here DEPTH = U(N_OLEVS+1)) is at 50
    DZ_TOP = DZ_TOP * (5000.0/DEPTH)
    DO I = 1,N_OLEVS
      DZ(I) = DZ(I) * (5000.0/DEPTH)
    ENDDO

    R1(1) = (KAPPA/RHOCP)*(DTIME/(DZ_TOP*DZ_TOP))
    R2(1) = 0.0
    DO I = 2,N_OLEVS
      R1(I) = (KAPPA/RHOCP)*(DTIME/(DZ(I-1)*(DZ(I)+DZ(I-1))))
      R2(I) = (KAPPA/RHOCP)*(DTIME/(DZ(I)*(DZ(I)+DZ(I-1))))
    ENDDO

! Reset the new values of U_OLD for the start of this run.
    DO I = 1,N_OLEVS
      U_OLD(I) = DTEMP_O(I)
    ENDDO

    DO I = 1,TIMESTEPS

      LAMBDA_OLD = -Q/(KAPPA*F_OCEAN)
      LAMBDA_NEW = -Q/(KAPPA*F_OCEAN)

      P = ((1.0-F_OCEAN)*LAMBDA_L*MU)/(F_OCEAN*KAPPA)             &
        + (LAMBDA_O/KAPPA)

      CALL INVERT(                                                &
        U_OLD,U_NEW,P,LAMBDA_OLD,LAMBDA_NEW,R1,R2,DZ_TOP,N_OLEVS  &
      )

! Now check that the flux out of the bottom is not too large.
      FLUX_TOP =  - 0.5*(KAPPA*LAMBDA_OLD*DTIME)                  &
                  - 0.5*(KAPPA*LAMBDA_NEW*DTIME)                  &
                  - 0.5*(P*KAPPA*U_OLD(1)*DTIME)                  &
                  - 0.5*(P*KAPPA*U_NEW(1)*DTIME)
      FLUX_BOTTOM = (DTIME / (2.0 * DZ(N_OLEVS))) * KAPPA         &
                  * (U_OLD(N_OLEVS) + U_NEW(N_OLEVS))

      IF(ABS(FLUX_BOTTOM/(FLUX_TOP+0.0000001)) > 0.01) THEN
        WRITE(jules_message,*) 'Flux at bottom of the ocean is greater ' //   &
                   'than 0.01 of top'
        CALL jules_print('delta_temp',jules_message)
      ENDIF


! Set calculated values at now ``old'' values
      DO J = 1,N_OLEVS
        U_OLD(J) = U_NEW(J)
      ENDDO

    ENDDO

! At end of this model run, reset values of LAMBDA_L and LAMBDA_O

    DO J = 1,N_OLEVS
      DTEMP_O(J) = U_NEW(J)
    ENDDO
    DTEMP_L = MU*DTEMP_O(1)

    RETURN

  END SUBROUTINE DELTA_TEMP
#endif
