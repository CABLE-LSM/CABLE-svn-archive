! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine to calculate albedos of land-surface tiles and gridbox-mean
! albedo for MOSES II.


SUBROUTINE tile_albedo (                                          &
 p_field,land_field,land_index,ntiles,tile_pts,                   &
 tile_index,l_aggregate,l_spec_albedo,albsoil,                    &
 cosz,frac,lai_in,rgrain,snow_tile,soot,tstar_tile,               &
 z0_tile,alb_tile,land_albedo,can_rad_mod                         &
 )

USE nstypes, ONLY : npft,ntype, urban_canyon, urban_roof
USE c_0_dg_c
USE nvegparm
USE pftparm
USE rad_param,  ONLY : kland,maskd,tcland
USE snow_param, ONLY : rho_snow_const                             &
                      ,cansnowtile
USE switches,   ONLY : l_point_data                               &
                      ,can_model                                  &
                      ,l_snowdep_surf, L_cable
USE switches_urban, ONLY : l_urban2t, l_moruses_albedo
USE urban_param, ONLY : albwl, albrd, hwr, ztm, albsnc_c,         &
                        albsnc_rf, albsnf_c, albsnf_rf

#if defined(UM_JULES)
USE jules_mod,   ONLY : smvcst_levs                               &
                       ,snowdep_surf
#else
USE p_s_parms, ONLY : smvcst_levs => smvcst
USE jules_mod, ONLY : snowdep_surf
#endif
USE prognostics, ONLY : snowdepth

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
!CABLE:
USE model_time_mod, ONLY : timestep
USE cable_vars_mod, ONLY : cosz_cable

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER                                                           &
 p_field                                                          &
                             ! Total number of grid points.
,land_field                                                       &
                             ! No. of land points.
,ntiles                      ! Number of surface tiles.

LOGICAL                                                           &
 l_aggregate                                                      &
                             ! IN Logical to set aggregate
!                                  !    surface scheme
,l_spec_albedo               ! .TRUE. for spectral albedos
!                                  ! and prognostic snow albedo.

!   Array arguments with intent(in):
INTEGER                                                           &
 land_index(land_field)                                           &
                             ! Index of land points.
,tile_pts(ntype)                                                  &
                             ! Number of land points which
!                                  ! include the nth surface type.
,tile_index(land_field,ntype)! Indices of land points which
!                                  ! include the nth surface type.

REAL                                                              &
 albsoil(land_field)                                              &
                             ! Soil albedo.
,cosz(p_field)                                                    &
                             ! Cosine of the zenith angle.
,frac(land_field,ntype)                                           &
                             ! Fractional cover of each
!                                  ! surface type.
,lai_in(land_field,npft)                                          &
                             ! Leaf area index.
,rgrain(land_field,ntiles)                                        &
                             ! Snow grain size on tiles
!                                  ! (microns).
,snow_tile(land_field,ntiles)                                     &
                             ! Canopy snow on tiles (kg/m2)
,soot(p_field)                                                    &
                             ! Snow soot content (kg/kg).
,tstar_tile(land_field,ntiles)                                    &
                             ! Tile surface temperatures (K).
,z0_tile(land_field,ntiles)  ! Surface roughness on tiles (m).

!   Array arguments with intent(out):
REAL                                                              &
 alb_tile(land_field,ntiles,4)                                    &
                              !Albedos for surface tiles.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR
,land_albedo(p_field,4)      ! GBM albedos.

! Local arrays:
REAL                                                              &
 albsnc(land_field,ntype)                                         &
                             ! Snow-covered albedo of surf types.
,albsnf(land_field,ntype)                                         &
                             ! Snow-free albedo of surf types.
,alb_type(land_field,ntype,4)                                     &
                             ! Albedos of surface types.
,alb_snow(land_field,ntype,4)                                     &
                             ! Snow albedos.
,fsnow(land_field)                                                &
                             ! Weighting factor for albedo.
,lai(land_field,npft)                                             &
                             ! Adjusted leaf area index.
,snowd(land_field)                                                &
                             ! Work variable (snow depth).
,tstar(land_field)                                                &
                             ! Copy of TSTAR_TILE.
,z0(land_field)              ! Copy of Z0_TILE.

INTEGER, PARAMETER ::       ilayers_dummy=1


! This variable is not used in this routine, but it is required in
! the subroutine argument list for the UM implementation of JULES
INTEGER                                                           &
 can_rad_mod                       ! Which canopy radiation model
                                   ! we're using



REAL                                                              &
 fapar_dir_dummy(land_field,npft,ilayers_dummy)                   &
!                                 ! Profile of absorbed PAR -
!                                 ! Direct beam - DUMMY
,fapar_dif_dummy(land_field,npft,ilayers_dummy)                   &
!                                 ! Profile of absorbed PAR -
!                                 ! Diffuse beam -DUMMY
,fapar_dir2dif_dummy(land_field,npft,ilayers_dummy)               &
!                                 ! DUMMY
,fapar_dif2dif_dummy(land_field,npft,ilayers_dummy)               &
!                                 ! DUMMY
,fapar_dir2dir_dummy(land_field,npft,ilayers_dummy)               &
!                                 ! DUMMY
,fsun_dummy(land_field,npft,ilayers_dummy)
!                                 ! DUMMY

! Local scalars:
REAL                                                              &
 dsa                                                              &
                             ! Deep-snow albedo.
,flit
                             ! Weighting factor for albedo.
INTEGER                                                           &
 band,i,j,l,n                ! Loop counters

LOGICAL                                                           &
 pointflag                   ! Switch for treatment of snow

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! declare parameters used in ccontrol
!#include "chsunits.h"
!#include "csubmmax.h"
! We need this file for CAN_MODEL
!#include "ccontrol.h"


IF (lhook) CALL dr_hook('TILE_ALBEDO',zhook_in,zhook_handle)
DO n=1,ntiles
  DO band=1,4
    DO l=1,land_field
      alb_tile(l,n,band) = 0.
    END DO
  END DO
END DO
DO n=1,ntype
  DO band=1,4
    DO l=1,land_field
      alb_type(l,n,band) = 0.
      alb_snow(l,n,band) = 0.
    END DO
  END DO
END DO

! Impose minimum LAI for bare vegetation
DO n=1,npft
  DO j=1,tile_pts(n)
    l = tile_index(j,n)
    lai(l,n) = MAX( lai_in(l,n), 0.5 )
  END DO
END DO

! Equivalent snowdepth for surface calculations.
snowdep_surf(:,:) = snowdepth(:,:)
DO n=1,ntiles
  IF ( (can_model==4).AND.cansnowtile(n).AND.l_snowdep_surf ) THEN
    DO l=1,land_field
      snowdep_surf(l,n)=snow_tile(l,n)/rho_snow_const
    END DO
  END IF
END DO

IF (l_spec_albedo) THEN
!----------------------------------------------------------------------
! Spectral albedo scheme with prognostic snow albedo
!----------------------------------------------------------------------

! Set albedos of vegetated surface types
! The logical argument (getProfile in albpft) is FALSE to indicate
! that profiles through the canopy should not be calculated.
! DEPENDS ON: albpft
  CALL albpft       (p_field,land_field,                          &
                     land_index,tile_index,tile_pts,              &
                     ilayers_dummy,.FALSE.,                       &
                     albsoil,cosz,lai,alb_type,                   &
                     fapar_dir_dummy,fapar_dif_dummy,             &
                     fapar_dir2dif_dummy,fapar_dif2dif_dummy,     &
                     fapar_dir2dir_dummy,fsun_dummy )

! Set albedos of non-vegetated surface types
  DO band=1,4
    DO n=npft+1,ntype
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        alb_type(l,n,band) = albsnf_nvg(n-npft)
        IF ( albsnf_nvg(n-npft) <  0. )                           &
                                               ! Soil tile
          alb_type(l,n,band) = albsoil(l)
      END DO
    END DO
  END DO

! Set snow-free albedos for urban_2T. Cannot have a canyon without a roof so
! these are done at the same time to avoid overwriting ice tiles. These are
! over-written if l_moruses_albedo.
  IF ( l_urban2t ) THEN
    n = urban_canyon
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      alb_type(l,n,:)          = albsnf_c
      alb_type(l,urban_roof,:) = albsnf_rf
    END DO
  END IF

! Calculate snow albedos
! DEPENDS ON: albsnow
  CALL albsnow(p_field,land_field,land_index,                     &
               ntiles,tile_index,tile_pts,l_aggregate,            &
               cosz,rgrain,snowdep_surf,soot,alb_snow)

! Adjust surface type albedos for snow cover
  DO l=1,land_field
    snowd(l) = snowdep_surf(l,1)
    z0(l) = z0_tile(l,1)
  END DO
  DO n=1,ntype
    IF (.NOT. l_aggregate) THEN
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        snowd(l) = snowdep_surf(l,n)
        z0(l) = z0_tile(l,n)

      END DO
    END IF
! Calculate snow albedo weighting factor.
    fsnow(:) = 0.0
    IF ( l_point_data ) THEN
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        IF ( snowd(l) .gt. 0.) fsnow(l) = 1.0 - EXP( -50.0*snowd(l) )
      END DO
    ELSE
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        IF ( snowd(l) .gt. 0.) fsnow(l) = snowd(l) / ( snowd(l) + 10.*z0(l) )
      END DO
    END IF
! Calculate weighted tile albedo.
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      DO band=1,4
        alb_type(l,n,band) = fsnow(l)*alb_snow(l,n,band)        &
                            + (1. - fsnow(l))*alb_type(l,n,band)
      END DO

! MORUSES: overwrite alb_type with values from canyonalb or albsnf_rf. Use
! test on smvcst_levs to check for land-ice points. Loop around band could be
! taken out as no snow involved at the moment.
      IF ( l_moruses_albedo ) THEN
        IF ( n == urban_canyon ) THEN
          i = land_index(l)
          DO band = 1,4
! DEPENDS ON: canyonalb
            CALL canyonalb(cosz(i),hwr(l),albwl(l),albrd(l),      &
               alb_type(l,n,band))
          END DO
        ELSE IF ( n == urban_roof                                 &
                  .AND. smvcst_levs(l,1) > 0.0 ) THEN
! MORUSES here removes the effects of snow; snow scheme for MORUSES needs
! to be added
          alb_type(l,n,:) = albsnf_rf
        END IF
      END IF

    END DO
  END DO   !  ntype

ELSE
!----------------------------------------------------------------------
! Non-spectral albedo scheme with diagnosed snow albedo
!----------------------------------------------------------------------

! Set albedos of vegetated surface types
  DO n=1,npft
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      flit = 1.0 - EXP(-kext(n)*lai(l,n))
      albsnc(l,n) = albsnc_min(n)*(1 - flit) + albsnc_max(n)*flit
      albsnf(l,n) = albsoil(l)*(1 - flit) + albsnf_max(n)*flit
    END DO
  END DO

! Set albedos of non-vegetated surface types
  DO n=npft+1,ntype
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      albsnc(l,n) = albsnc_nvg(n-npft)
      albsnf(l,n) = albsnf_nvg(n-npft)
      IF ( albsnf_nvg(n-npft) <  0. ) albsnf(l,n) = albsoil(l)
    END DO
  END DO

! Set canyon & roof albedos for two-tile urban at the same time to avoid
! over-writing land-ice points. Note alb_type is over-written if
! l_moruses_albedo.
  IF ( l_urban2t ) THEN
    n = urban_canyon
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      albsnc(l,n) = albsnc_c
      albsnf(l,n) = albsnf_c
      albsnc(l,urban_roof) = albsnc_rf
      albsnf(l,urban_roof) = albsnf_rf
    END DO
  END IF

! Adjust surface type albedos for snow cover
  DO l=1,land_field
    tstar(l) = tstar_tile(l,1)
    snowd(l) = snowdep_surf(l,1)
  END DO
  DO n=1,ntype
    IF (.NOT. l_aggregate) THEN
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        tstar(l) = tstar_tile(l,n)
        snowd(l) = snowdep_surf(l,n)
      END DO
    END IF
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      IF ( tstar(l)  <   tcland ) THEN
        dsa = albsnc(l,n)
      ELSE IF ( tstar(l)  <   tm ) THEN
        dsa = albsnc(l,n) + kland*(albsnf(l,n) - albsnc(l,n))     &
                                 *(tstar(l) - tcland)
      ELSE
        dsa = albsnc(l,n) + kland*(albsnf(l,n) - albsnc(l,n))     &
                                 *(tm - tcland)
      END IF
      alb_type(l,n,1) = albsnf(l,n) + (dsa - albsnf(l,n)) *       &
                                    ( 1. - EXP(-maskd*snowd(l)) )
    END DO
  END DO

  IF ( l_moruses_albedo ) THEN
    n = urban_canyon
    DO j = 1,tile_pts(n)
      l = tile_index(j,n)
      i = land_index(l)
      ! DEPENDS ON: canyonalb
      CALL canyonalb( cosz(i), hwr(l), albwl(l), albrd(l),        &
         alb_type(l,n,1) )
! MORUSES here removes the effects of snow; snow scheme for MORUSES needs
! to be added
      alb_type(l,urban_roof,1) = albsnf_rf
    END DO
  END IF

! Copy albedo to all bands
  DO band=2,4
    DO n=1,ntype
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        alb_type(l,n,band) = alb_type(l,n,1)
      END DO
    END DO
  END DO

END IF       ! Spectral or non-spectral albedo schemes

!----------------------------------------------------------------------
! Calculate GBM surface albedo
!----------------------------------------------------------------------

DO band=1,4
  DO i=1,p_field
    land_albedo(i,band) = 0.
  END DO
  DO n=1,ntype
    DO j=1,tile_pts(n)
      l = tile_index(j,n)
      i = land_index(l)
      land_albedo(i,band) = land_albedo(i,band) +                 &
                            frac(l,n)*alb_type(l,n,band)
    END DO
  END DO
END DO

!----------------------------------------------------------------------
! Copy albedos as required for aggregate or distinct tiles
!----------------------------------------------------------------------

IF (l_aggregate) THEN
  DO band=1,4
    DO l=1,land_field
      i = land_index(l)
      alb_tile(l,1,band) = land_albedo(i,band)
    END DO
  END DO
ELSE
  DO band=1,4
    DO n=1,ntype
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        alb_tile(l,n,band) = alb_type(l,n,band)
      END DO
    END DO
  END DO
END IF

IF ( l_cable ) then

! land_alb n cable_RADNum is only set on points where there is incoming
! SW radiation at the previous timestep, therefore it will be zero over some
! land points
          if(timestep .gt. 1) then
! DEPENDS ON: cable_rad_driver.o
                  CALL cable_rad_driver(                &
! IN atmospheric forcing
!     &           surf_down_sw, cos_zenith_angle, &
     &           cosz_cable,                       &
! IN soil/vegetation/land surface data :
     !&           SNOW_TILE,SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,            &
     !&           ISNOW_FLG3L,                                           &
! IN  Time stepping infomation:
     &           albsoil,                            &
     !&           LAND_ALBEDO,ALB_TILE,LAND_ALB                          &
     &           ALB_TILE&
     &  )
         endif  ! if if(istep_cur .gt. 1)
endif  ! if l_cable



IF (lhook) CALL dr_hook('TILE_ALBEDO',zhook_out,zhook_handle)
RETURN
END SUBROUTINE tile_albedo
