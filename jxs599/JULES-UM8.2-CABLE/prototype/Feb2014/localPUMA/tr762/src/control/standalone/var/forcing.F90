#if !defined(UM_JULES)
! Module containing all of the driving (atmospheric forcing) variables.

MODULE forcing

  IMPLICIT NONE

!-------------------------------------------------------------------------------

! The forcing variables.
  REAL, ALLOCATABLE ::                                                        &
    qw_1(:,:),                                                                &
                        !  Total water content (Kg/Kg)
    tl_1(:,:),                                                                &
                        ! Ice/liquid water temperature (k)
    u_0(:,:),                                                                 &
                        ! W'ly component of surface current (m/s)
    v_0(:,:),                                                                 &
                        ! S'ly component of surface current (m/s)
    u_1(:,:),                                                                 &
                        ! W'ly wind component (m/s)
    v_1(:,:),                                                                 &
                        ! S'ly wind component (m/s)
    pstar(:,:),                                                               &
                        ! Surface pressure (Pascals)
    ls_rain(:,:),                                                             &
                        ! Large-scale rain (kg/m2/s)
    con_rain(:,:),                                                            &
                        ! Convective rain (kg/m2/s)
    ls_snow(:,:),                                                             &
                        ! Large-scale snowfall (kg/m2/s)
    con_snow(:,:),                                                            &
                        ! Convective snowfall (kg/m2/s)
    sw_down(:,:),                                                             &
                        ! Surface downward SW radiation (W/m2)
    lw_down(:,:)        ! Surface downward LW radiation (W/m2)

! Variables that aid in the calculation of the actual forcing variables
  REAL, ALLOCATABLE ::                                                        &
    diff_rad(:,:)       ! Input diffuse radiation (W/m2)

!-------------------------------------------------------------------------------

END MODULE forcing
#endif
