#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_FIT_FSAT-------------------------------------------
!
!   Purpose: To speed up the large scale hydrology code (LTOP=TRUE)
!            dramatically. This is done by fitting exponential
!            functions to the incomplete gamma function for each grid
!            box and the complete range of possible "water table"
!            (top_crit) cases - see documentation.
!            Estimates the fitted parameters for Fsat=function(ZW)
!            and  Fwet=function(ZW) for each land grid point.
!            (Calculating the incomplete gamma function for each grid
!            box at each time step was very time consuming).
!                                                             !
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  6.4    10/01/07   New Deck         Nic Gedney
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.

SUBROUTINE calc_fit_fsat()

  !Use in relevant subroutines
USE calc_fsat_mod, ONLY: calc_fsat

!Use in relevant varaibles
USE ancil_info,           ONLY:                                               &
  land_pts, soil_pts, soil_index, nsoilt

USE top_pdm,              ONLY:                                               &
  a_fsat_soilt, a_fwet_soilt, c_fsat_soilt, c_fwet_soilt,  fexp_soilt,        &
  gamtot_soilt, ti_mean_soilt, ti_sig_soilt

USE jules_hydrology_mod,  ONLY:                                               &
  l_top, zw_max, nfita

USE jules_soil_mod,       ONLY:                                               &
  dzsoil, sm_levels

USE jules_print_mgr,      ONLY:                                               &
  jules_message, jules_print

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Subroutine arguments - NONE.
!-----------------------------------------------------------------------------

! Work variables
REAL(KIND=real_jlslsm) :: zdepth           ! Standard Soil model DEPTH.

! Local scalars:
INTEGER, PARAMETER :: nzw = 100  ! Number of ZW values used in fit.
                                 ! Maximum value for a significant
                                 ! improvement in the fit.

INTEGER :: i,j,iz,n,ifita,m  ! Loop counters.

REAL(KIND=real_jlslsm) ::                                                     &
  dzw,                                                                        &
                 ! WORK ZW increment ; defined by zw_max and nzw.
  rms,                                                                        &
                 ! WORK rms errors for given fsat fit values.
  rmsw,                                                                       &
                 ! WORK rms errors for given fwet fit values.
  rmsold,                                                                     &
                 ! WORK rms errors for given fsat fit values.
                 !      for best fit so far.
  rmswold,                                                                    &
                 ! WORK rms errors for given fwet fit values
                 !      for best fit so far.
  cfit,                                                                       &
                 ! WORK CFit value for given loop.
  cfitmax
                 ! Maximum possible value for Cfit.

REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  thr_err = 5.0e-3,                                                           &
                 ! Error threshold value

  cfitmin = 0.0
                 ! Minimum possible value for Cfit.

! Local arrays:
REAL(KIND=real_jlslsm) ::                                                     &
  fsat_calc_soilt(land_pts,nsoilt,nzw),                                       &
                            ! WORK Surface saturation fraction.
  fsat_fit(nzw),                                                              &
                            ! WORK Fitted surface saturation fraction.
  fwet_calc_soilt(land_pts,nsoilt,nzw),                                       &
                            ! WORK Wetland fraction.
  fwet_fit(nzw),                                                              &
                            ! WORK Fitted wetland fraction.
  dumzw(nzw),                                                                 &
                            ! WORK Dummy water table depth (m).
  dumfsat_soilt(land_pts,nsoilt),                                             &
                            ! WORK Dummy surface saturation fraction.
  dumfwetl_soilt(land_pts,nsoilt)
                            ! WORK Dummy wetland fraction.
! DBC We could use this local version of gamtot in this version,
!     but using the module version to match other versions of code.
!     &,gamtot_soilt(land_pts)       ! WORK Integrated complete Gamma function

REAL(KIND=real_jlslsm) ::                                                     &
  top_crit_soilt(land_pts,nsoilt,nzw),                                        &
                            ! WORK LOG(QBASE_MAX/QBASE) -see document.
  top_crit1z_soilt(land_pts,nsoilt),                                          &
                            ! WORK As above but for an individual zw.
  top_min_soilt(land_pts,nsoilt),                                             &
                            ! WORK value for when zw=zw_max.
  wutot_soilt(land_pts,nsoilt)           ! WORK Dummy (set to 1.0).

INTEGER :: errorstatus
CHARACTER(LEN=80) :: cmessage
CHARACTER(LEN=13) :: routinename


!-----------------------------------------------------------------------------
routinename='CALC_FIT_FSAT'

IF ( l_top ) THEN

  cfitmax = 0.15 * nfita

  ! Define the water table depths to be used in the fitting process:
  dzw = 1.0 / REAL(nzw) * zw_max
  DO iz = 1,nzw
    dumzw(iz) = REAL(iz-1) * dzw
  END DO

  DO i = 1,land_pts        ! initialise to zero
    wutot_soilt(i,:)    = 1.0
    dumfsat_soilt(i,:)  = 0.0
    dumfwetl_soilt(i,:) = 0.0
  END DO
  zdepth = 0.0

  ! Calculate total soil depth
  DO n = 1,sm_levels
    zdepth = zdepth + dzsoil(n)
  END DO

  ! Calculate Gamtot
  DO m = 1,nsoilt
    gamtot_soilt(:,m) = 0.0
    CALL calc_fsat( .TRUE., soil_pts, soil_index, land_pts,                   &
                  ti_mean_soilt(:,m),                                         &
                  ti_sig_soilt(:,m), wutot_soilt(:,m),                        &
                  top_crit1z_soilt(:,m), gamtot_soilt(:,m),                   &
                  dumfsat_soilt(:,m), dumfwetl_soilt(:,m))
  END DO

  ! Calculate top_crit for the water table depths:
  DO iz = 1,nzw
    DO m = 1, nsoilt
      DO j = 1,soil_pts
        i = soil_index(j)

        fsat_calc_soilt(i,m,iz) = 0.0
        fwet_calc_soilt(i,m,iz) = 0.0
        top_crit_soilt(i,m,iz)  = 0.0

        IF ( ti_mean_soilt(i,m) > 0.0 .AND. ti_sig_soilt(i,m) > 0.0 ) THEN
          top_min_soilt(i,m) = 1.0 / fexp_soilt(i,m)                          &
                               * EXP(-fexp_soilt(i,m) * (zw_max - zdepth))

          IF ( dumzw(iz) <= zdepth ) THEN
            top_crit1z_soilt(i,m) = LOG(zdepth+1.0                            &
                                      / fexp_soilt(i,m) - top_min_soilt(i,m)) &
                                  - LOG(zdepth - dumzw(iz) + 1.0              &
                                      / fexp_soilt(i,m) - top_min_soilt(i,m))
          END IF

          IF ( dumzw(iz) >  zdepth ) THEN
            top_crit1z_soilt(i,m) = LOG(zdepth+1.0                            &
                                    / fexp_soilt(i,m) - top_min_soilt(i,m))   &
                                    - LOG(1 / fexp_soilt(i,m)                 &
                                    * EXP(-fexp_soilt(i,m)                    &
                                    * (dumzw(iz) - zdepth))                   &
                                    - top_min_soilt(i,m))
          END IF
        END IF
      END DO !soil_pts

      ! Calculate FSAT and FWET for one ZW at all soil land_pts:
      CALL calc_fsat( .FALSE., soil_pts, soil_index, land_pts,                &
                     ti_mean_soilt(:,m),                                      &
                     ti_sig_soilt(:,m), wutot_soilt(:,m),                     &
                     top_crit1z_soilt(:,m), gamtot_soilt(:,m),                &
                     dumfsat_soilt(:,m), dumfwetl_soilt(:,m))

      DO j = 1,soil_pts
        i = soil_index(j)

        IF ( ti_mean_soilt(i,m) > 0.0 .AND. ti_sig_soilt(i,m) > 0.0 ) THEN
          fsat_calc_soilt(i,m,iz) = dumfsat_soilt(i,m)
          fwet_calc_soilt(i,m,iz) = dumfwetl_soilt(i,m)
          top_crit_soilt(i,m,iz)  = top_crit1z_soilt(i,m)

          IF ( dumzw(iz) <  dzw ) THEN ! Values at zw=0m
            a_fsat_soilt(i,m) = fsat_calc_soilt(i,m,iz)
            a_fwet_soilt(i,m) = fwet_calc_soilt(i,m,iz)
          END IF
        END IF
      END DO
    END DO !nsoilt
  END DO                  !ZW calc_fsat loop

  ! Now carry out fit for FSAT, where FSAT=function(ZW). (Likewise FWET)
  DO m = 1, nsoilt
    DO j = 1,soil_pts
      i = soil_index(j)

      IF ( ti_mean_soilt(i,m) > 0.0 .AND. ti_sig_soilt(i,m) > 0.0 ) THEN
        rmsold  = 1.0e10
        rmswold = 1.0e10

        DO ifita = 1,nfita
          cfit = cfitmax * (ifita) / REAL(nfita)

          rms  = 0.0
          rmsw = 0.0
          ! top_crit=TI_MAX when zw=zw_max
          DO iz = 1,nzw
            fsat_fit(iz) = a_fsat_soilt(i,m)                                  &
                            * EXP(-cfit * top_crit_soilt(i,m,iz))
            fwet_fit(iz) = a_fwet_soilt(i,m)                                  &
                            * EXP(-cfit * top_crit_soilt(i,m,iz))
            rms  = rms  + (fsat_calc_soilt(i,m,iz) - fsat_fit(iz))**2
            rmsw = rmsw + (fwet_calc_soilt(i,m,iz) - fwet_fit(iz))**2
          END DO         !ZW

          rms  = SQRT(rms) / REAL(nzw)
          rmsw = SQRT(rmsw) / REAL(nzw)

          IF ( rms < rmsold ) THEN
            rmsold = rms
            c_fsat_soilt(i,m) = cfit
          END IF
          IF ( rmsw < rmswold ) THEN
            rmswold = rmsw
            c_fwet_soilt(i,m) = cfit
          END IF
        END DO

        DO iz = 1,nzw
          fsat_fit(iz) = a_fsat_soilt(i,m) * EXP(-c_fsat_soilt(i,m) *         &
                         top_crit_soilt(i,m,iz))
          fwet_fit(iz) = a_fwet_soilt(i,m) * EXP(-c_fwet_soilt(i,m) *         &
                         top_crit_soilt(i,m,iz))
        END DO            !ZW

        IF ( rmsold >= thr_err ) THEN
          IF ( c_fsat_soilt(i,m) <= cfitmin .OR.                              &
              c_fsat_soilt(i,m) >= cfitmax ) THEN

            WRITE(jules_message,*) 'ERROR cfit FSAT', i,                      &
                c_fsat_soilt(i,m), cfitmin, cfitmax
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*)                                            &
              'If c_fsat=cfitmax try increasing nfita'
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fsat_calc=',                              &
                fsat_calc_soilt(i,m,1), fsat_calc_soilt(i,m,3),               &
                fsat_calc_soilt(i,m,5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fsat_fit=',                               &
                fsat_fit(1), fsat_fit(3), fsat_fit(5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'rms=', rmsold
            CALL jules_print('calc_fit_fsat',jules_message)
            ErrorStatus = 35
            WRITE(CMessage, *) 'Error in cfit FSAT in LSH model setup'
          END IF
        END IF

        IF ( rmswold >= thr_err ) THEN
          IF ( c_fwet_soilt(i,m) <= cfitmin .OR.                              &
              c_fwet_soilt(i,m) >= cfitmax ) THEN

            WRITE(jules_message,*) 'ERROR cfit FWET', i,                      &
                c_fwet_soilt(i,m), cfitmin, cfitmax
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fwet_calc=',                              &
                fwet_calc_soilt(i,m,1), fwet_calc_soilt(i,m,3),               &
                fwet_calc_soilt(i,m,5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fwet_fit=',                               &
                fwet_fit(1), fwet_fit(3), fwet_fit(5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'rmsw=', rmswold
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) '(fsat_calc=)',                            &
                fsat_calc_soilt(i,m,1), fsat_calc_soilt(i,m,3),               &
                fsat_calc_soilt(i,m,5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) '(fsat_fit=)',                             &
                fsat_fit(1), fsat_fit(3), fsat_fit(5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) '(rms=)', rmsold
            CALL jules_print('calc_fit_fsat',jules_message)
            ErrorStatus = 35
            ErrorStatus = 40
            WRITE(CMessage,*) 'Error in cfit FWET in LSH model setup'
          END IF
        END IF

        IF ( rmsold > thr_err ) THEN
          WRITE(jules_message,*) 'Warning LSH rms Error in fit:',             &
              rmsold, rmswold
          CALL jules_print('calc_fit_fsat',jules_message)
        END IF
      END IF

    END DO                  ! land_pts
  END DO !nsoilt
ELSE
  DO m = 1,nsoilt
    DO j = 1,soil_pts
      i = soil_index(j)

      a_fsat_soilt(i,m) = 0.0
      c_fsat_soilt(i,m) = 0.0
      a_fwet_soilt(i,m) = 0.0
      c_fwet_soilt(i,m) = 0.0
      gamtot_soilt(i,m) = 0.0
    END DO
  END DO
END IF

RETURN

END SUBROUTINE calc_fit_fsat
#endif
