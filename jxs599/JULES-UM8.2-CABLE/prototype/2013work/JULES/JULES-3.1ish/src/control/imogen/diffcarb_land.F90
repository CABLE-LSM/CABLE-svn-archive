#if !defined(UM_JULES)
!**********************************************************************
! Calculates the global change in atmospheric CO2 due to 
! grid box changes in total carbon content.
!
! This code corresponds tothe output of HADCM3, with a longitudinal
! spacing of 3.75 degrees (96 points around a given latitudinal "band"
! and a latitude spacing of 2.5 degrees with a top and bottom box
! of 1.25 degrees.
!
! Written by C. Huntingford (November 1999)
!**********************************************************************
  SUBROUTINE DIFFCARB_LAND(LAND_PTS,D_LAND_ATMOS,LAT,DCTOT,CONV)

    USE c_pi, ONLY : pi
    
    USE earth_constants_mod, ONLY : RAD_EARTH => earth_radius

    IMPLICIT NONE

    INTEGER ::                                                    &
      LAND_PTS,                                                   &
                   ! IN Number of land points
      L            !    Looping parameter

    REAL ::                                                       &
      D_LAND_ATMOS,                                               &
                   ! OUT Change in atmospheric CO2 concentration
                   !     as a result of land-atmosphere feedbacks
                   !     between calls to SCENARIO (ppm/"YEAR_CO2")
      LAT(LAND_PTS),                                              &
                   ! IN Latitute (degrees)
      DCTOT(LAND_PTS),                                            &
                   ! IN Change in total surface gridbox carbon
                   !    content during a period "YEAR_CO2" (kg C/m2)
      CONV,                                                       &
                   ! IN Converts global emission of C (Gt/yr)
                   !    into atmospheric CO2 (ppm/GtC)
      AREA_EARTH,                                                 &
                   ! WORK Surface area of the Earth (m2)
      DAREA,                                                      &
                   ! WORK Gridbox area (m2)
      LAT_MAX,                                                    &
                   ! Latitute of "top" of GCM box (degrees)
      LAT_MIN,                                                    &
                   ! Latitute of "bottom" of GCM box (degrees)
      LAT_MAX_RAD,                                                &
                   ! Latitute of "top" of GCM box (rad)
      LAT_MIN_RAD,                                                &
                   ! Latitute of "bottom" of GCM box (rad)
      LAND_GAIN    ! Total gain in carbon by land (kg C)

    
    
    LAND_GAIN = 0.0
    AREA_EARTH = 4.0*PI*(RAD_EARTH**2)

    DO L = 1,LAND_PTS
      IF(LAT(L) >= 88.75) THEN
        LAT_MAX = 90.0
        LAT_MIN = 88.75
      ELSE IF(LAT(L) <= -88.75) THEN
        LAT_MAX = -88.75
        LAT_MIN = -90.0
      ELSE
        LAT_MAX = LAT(L) + 1.25
        LAT_MIN = LAT(L) - 1.25
      ENDIF
      
      LAT_MAX_RAD = (LAT_MAX*PI)/180.0
      LAT_MIN_RAD = (LAT_MIN*PI)/180.0

      DAREA =                                                     &
        (2*PI*(RAD_EARTH**2)*(SIN(LAT_MAX_RAD)-SIN(LAT_MIN_RAD)))/96.0

      LAND_GAIN = LAND_GAIN + DAREA*DCTOT(L)

    ENDDO

    D_LAND_ATMOS = -(LAND_GAIN/1.0E12)*CONV

    RETURN

  END SUBROUTINE DIFFCARB_LAND
#endif
