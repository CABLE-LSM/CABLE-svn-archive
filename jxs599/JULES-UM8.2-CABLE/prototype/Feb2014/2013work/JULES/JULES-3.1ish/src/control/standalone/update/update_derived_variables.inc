#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/update/update_derived_variables.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE update_derived_variables()

  USE datetime_mod, ONLY : SECS_IN_DAY

  USE model_time_mod, ONLY : current_time, timestep_len

  USE switches, ONLY : l_imogen, l_point_data, l_aggregate, can_model

  USE ancil_info, ONLY : land_pts, land_index, ntiles, tile_pts, tile_index,  &
                         frac

  USE theta_field_sizes, ONLY : t_i_length

  USE forcing

  USE u_v_grid

  USE imogen_drive_vars

  USE prognostics, ONLY : canht_ft, lai

  USE p_s_parms, ONLY : satcon, catch_snow, catch, infil_tile, z0_tile

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Updates variables that are derived from those given in time-varying files
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: insd  ! Timestep in day - used to index IMOGEN arrays

  INTEGER :: i,j,l  ! Index variables


!------------------------------------------------------------------------------


  IF( l_imogen ) THEN
!-------------------------------------------------------------------------------
! If IMOGEN is enabled, copy the correct timestep of the current climatology
! into the driving variables
!-------------------------------------------------------------------------------
! Get the timestep in the day that we are on
    insd = (current_time%time / timestep_len) + 1

    DO l = 1,land_pts
      j = (land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j - 1) * t_i_length

      pstar(i,j)    = PSTAR_OUT(l,current_time%month,current_time%day,insd)
      u_1(i,j)      = WIND_OUT(l,current_time%month,current_time%day,insd)
      v_1(i,j)      = 0.0
      u_0(i,j)      = 0.0
      v_0(i,j)      = 0.0
      con_rain(i,j) = CONV_RAIN_OUT(l,current_time%month,current_time%day,insd) &
                                                             / REAL(SECS_IN_DAY)
      con_snow(i,j) = CONV_SNOW_OUT(l,current_time%month,current_time%day,insd) &
                                                             / REAL(SECS_IN_DAY)
      ls_rain(i,j)  = LS_RAIN_OUT(l,current_time%month,current_time%day,insd)   &
                                                             / REAL(SECS_IN_DAY)
      ls_snow(i,j)  = LS_SNOW_OUT(l,current_time%month,current_time%day,insd)   &
                                                             / REAL(SECS_IN_DAY)
      sw_down(i,j)  = SW_OUT(l,current_time%month,current_time%day,insd)
      lw_down(i,j)  = LW_OUT(l,current_time%month,current_time%day,insd)
      qw_1(i,j)     = QHUM_OUT(l,current_time%month,current_time%day,insd)
      tl_1(i,j)     = T_OUT(l,current_time%month,current_time%day,insd)
    ENDDO

  ELSE
!-----------------------------------------------------------------------------
! Otherwise, update the main driving variables based on what was given in the
! input file(s) - the update of those variables has already taken place
!-----------------------------------------------------------------------------
! Initialise wind based on io_wind_speed
! This is just a case of initialising variables that are not set via file
! u_1 is always set from file, v_1 is only set from file if using both
! wind components
    IF ( io_wind_speed ) v_1(:,:) = 0.0
    u_0(:,:) = 0.0
    v_0(:,:) = 0.0


!-----------------------------------------------------------------------------
! Update precipitation variables
! What we need to do depends on io_precip_type
!-----------------------------------------------------------------------------
    SELECT CASE ( io_precip_type )
      CASE ( 1 )
! Total precipitation is given in file and stored in ls_rain

! First partition into rain and snow - all snow is assumed to be large-scale
        ls_snow(:,:) = 0.0
        WHERE ( tl_1(:,:) <= t_for_snow )
          ls_snow = ls_rain
          ls_rain = 0.0
        END WHERE

! Now ls_rain contains total rainfall
! If using point data assume all rain is large-scale
        con_rain(:,:) = 0.0
        IF ( .NOT. l_point_data ) THEN
! Otherwise partition into convective and large scale based on t_for_con_rain
          WHERE ( tl_1(:,:) >= t_for_con_rain )
            con_rain = ls_rain
            ls_rain  = 0.0
          END WHERE
        END IF

! Convective snow is assumed to be 0 unless explicitly provided
        con_snow(:,:) = 0.0

      CASE ( 2 )
! Total rainfall given in file and currently stored in ls_rain

! If using point data assume all rain is large-scale
        con_rain(:,:) = 0.0
        IF ( .NOT. l_point_data ) THEN
! Otherwise partition into convective and large scale based on t_for_con_rain
          WHERE ( tl_1(:,:) >= t_for_con_rain )
            con_rain = ls_rain
            ls_rain  = 0.0
          END WHERE
        END IF

! Total snow was given in file, and all assumed to be large scale

! Convective snow is assumed to be 0 unless explicitly provided
        con_snow(:,:) = 0.0

      CASE ( 3 )
! Convective and large-scale rain given in file

! Total snow was given in file, and all assumed to be large-scale

! Convective snow is assumed to be 0 unless explicitly provided
        con_snow(:,:) = 0.0

      CASE ( 4 )
! There is nothing to be done for case 4, since all components are input from
! file
    END SELECT

!-----------------------------------------------------------------------------
! Radiation variables are updated in CONTROL, since the new albedos need to be
! known before the update takes place
!-----------------------------------------------------------------------------
  END IF

!-------------------------------------------------------------------------------
! Copy information to U, V and T grids (assume that att grids are the same)
!-------------------------------------------------------------------------------
  u_0_p(:,:) = u_0(:,:)
  v_0_p(:,:) = v_0(:,:)
  u_1_p(:,:) = u_1(:,:)
  v_1_p(:,:) = v_1(:,:)


!-----------------------------------------------------------------------------
! Update any variables dependent on variables being prescribed that are not
! driving variables
!-----------------------------------------------------------------------------
  IF ( have_prescribed_veg )                                                  &
    CALL SPARM(LAND_PTS, NTILES, CAN_MODEL, L_AGGREGATE, TILE_PTS,            &
               TILE_INDEX, FRAC, CANHT_FT, LAI, SATCON, CATCH_SNOW, CATCH,    &
               INFIL_TILE, Z0_TILE)

  RETURN

END SUBROUTINE update_derived_variables
#endif
