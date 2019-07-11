#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine JULES_INIT ----------------------------------
!
! Description: Initialisation of JULES.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land


SUBROUTINE jules_init(  &
#include "arglndm.h"
                        land_field,ntiles,sm_levels                   &
                       ,nice,nice_use,frac_land                       &
                       ,l_snow_albedo                                 &
                       ,rgrain,snow_tile,deep_soil_temp,snow_grnd     &
                       ,nsnow,rgrainl,rho_snow_grnd,sice,sliq         &
                       ,snowdepth,ds,tsnow )

USE nstypes, ONLY : ntype, lake

! JULES variables to be initialised
USE jules_mod, ONLY :  clapp_levs                                     &
                     , sathh_levs                                     &
                     ,  hcap_levs                                     &
                     ,  hcon_levs                                     &
                     ,satcon_levs                                     &
                     ,smvccl_levs                                     &
                     ,smvcwt_levs                                     &
                     ,smvcst_levs

! UM variables from which to initialise the JULES variables
USE atm_fields_mod, ONLY : clapp_horn                                 &
                          ,sat_soilw_suction                          &
                          ,sat_soil_cond                              &
                          ,therm_cap                                  &
                          ,therm_cond                                 &
                          ,vol_smc_crit                               &
                          ,vol_smc_wilt                               &
                          ,vol_smc_sat

! JULES module snow depth
USE prognostics, ONLY :                                               &
        snowdepth_jules => snowdepth

! JULES switches
      USE switches, ONLY : l_um_jules                                 &
                          ,l_aggregate                                &
                          ,l_flake_model

! JULES max no. of snow layers
  USE ancil_info, ONLY :  &
     nsmax

  USE dyn_coriolis_mod,  ONLY : f3_at_u
  USE theta_field_sizes, ONLY : t_i_length

  USE lake_mod, ONLY :    &
     coriolis_param       &
    ,nusselt              &
    ,nusselt_0

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine arguments

!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) ::                                                &
  land_field,                                                         &
                 ! Number of land points
  ntiles,                                                             &
                 ! Number of surface tiles
  sm_levels,                                                          &
                 ! Number of soil layers
  nice,                                                               &
                 ! Number of sea ice categories
  nice_use 
                 ! Number of sea ice cats used in radiation and
                 !  explicit part of surface exchange      

  REAL ::  frac_land( land_field,ntype)       &
!
          ,rgrain(    land_field,ntiles)      &
          ,snow_tile( land_field,ntiles)      &
          ,snow_grnd( land_field,ntiles)      &
!
          ,deep_soil_temp(land_field,sm_levels)  &
!
          ,snowdepth(    land_field,ntiles)          &
          ,nsnow(        land_field,ntiles)          &
          ,rho_snow_grnd(land_field,ntiles)          &
!
          ,tsnow(        land_field,ntiles,nsmax)    &
          ,rgrainl(      land_field,ntiles,nsmax)    &
          ,sice(         land_field,ntiles,nsmax)    &
          ,sliq(         land_field,ntiles,nsmax)    &
          ,ds(           land_field,ntiles,nsmax)    &
          ,rho_snow(     land_field,ntiles,nsmax)

! land_index
#include "typlndm.h"

 LOGICAL :: l_snow_albedo

! WORK variables:
INTEGER :: i,j,k,l,n

  INTEGER ::                                &
    nsnow_integer(land_field,ntiles)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!
! START OF EXECUTABLE CODE
!

! Dimension the JULES fields.
IF (lhook) CALL dr_hook('JULES_INIT',zhook_in,zhook_handle)
ALLOCATE( clapp_levs(land_field,  sm_levels))
ALLOCATE( sathh_levs(land_field,  sm_levels))
ALLOCATE(  hcap_levs(land_field,  sm_levels))
ALLOCATE(  hcon_levs(land_field,0:sm_levels))
ALLOCATE(satcon_levs(land_field,0:sm_levels))
ALLOCATE(smvccl_levs(land_field,  sm_levels))
ALLOCATE(smvcwt_levs(land_field,  sm_levels))
ALLOCATE(smvcst_levs(land_field,  sm_levels))

#if ! defined(SCMA)
! Initialise JULES arrays.
! Soil properties on soil moisture levels.
DO i = 1, land_field
  DO j = 1, sm_levels
      clapp_levs(i,j) = clapp_horn(       i)
      sathh_levs(i,j) = sat_soilw_suction(i)
       hcap_levs(i,j) = therm_cap(        i)
     smvccl_levs(i,j) = vol_smc_crit(     i)
     smvcwt_levs(i,j) = vol_smc_wilt(     i)
     smvcst_levs(i,j) = vol_smc_sat(      i)
   END DO
   DO j = 0, sm_levels
       hcon_levs(i,j) = therm_cond(   i)
     satcon_levs(i,j) = sat_soil_cond(i)
  END DO
END DO
#endif

#if defined(UM_JULES)
IF (l_um_jules) THEN

 PRINT*,'JULES_INIT: nsmax = ',nsmax

  IF (nsmax > 0) THEN
! initialise multi-layer snow variables
!---------------------------------------
! NOTE that this initialisation is for the pointers
! to the D1 array, NOT the variables in the JULES
! prognostics(etc) modules. The JULES variables are
! mapped to the UM prognostics in AP2. However an
! initial value of SNOWDEPTH is needed before then,
! see below.
!---------------------------------------
!
  DO i = 1, land_field
    DO n = 1, ntiles
      nsnow_integer(i,n) = NINT(nsnow(i,n))
    END DO
  END DO

! DEPENDS ON: total_snow_um
        CALL total_snow_um(  land_field     &
                            ,ntype          &
                            ,ntiles         &
                            ,frac_land      &
                            ,sm_levels      &
                            ,l_snow_albedo  &
                            ,rgrain         &
                            ,snow_tile      &
                            ,deep_soil_temp &
                            ,snow_grnd      &
                            ,nsnow_integer  &
                            ,rgrainl        &
                            ,rho_snow_grnd  &
                            ,sice           &
                            ,sliq           &
                            ,snowdepth      &
                            ,ds             &
                            ,tsnow  )
  DO i = 1, land_field
    DO n = 1, ntiles
      nsnow(i,n) = REAL(nsnow_integer(i,n))
    END DO
  END DO

  END IF

! snowdepth needed in AP1 for JULES radiation
!---------------------------------------------
  snowdepth_jules = snowdepth

END IF
#endif

! FLake model
!--------------
IF (     l_flake_model                       &
    .AND.(.NOT.l_aggregate)) THEN

! initialise the Nusselt number
    nusselt(:) = nusselt_0

    DO l=1,land_field

        j=(land_index(l)-1)/t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length

! set the Coriolis parameter : ABSOLUTE VALUE
!
! To get the value at theta points,
! average the adjacent values at u points.
!
        coriolis_param(l) = ABS( (f3_at_u(i,j)+f3_at_u(i-1,j))/2.0 )

    END DO

END IF

IF (lhook) CALL dr_hook('JULES_INIT',zhook_out,zhook_handle)
RETURN

END SUBROUTINE jules_init
#endif
