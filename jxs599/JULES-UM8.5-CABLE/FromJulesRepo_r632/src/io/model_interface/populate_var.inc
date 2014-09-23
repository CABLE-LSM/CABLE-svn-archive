#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE populate_var(var_id, cube, const_val)

  USE string_utils_mod, ONLY : to_string

  USE data_cube_mod, ONLY : data_cube, cube_get_data, cube_free

  USE ancil_info, ONLY : land_pts, sm_levels, frac, ntiles

  USE model_grid_mod, ONLY : latitude, longitude

  USE coastal, ONLY : flandg

  USE c_elevate, ONLY : surf_hgt

  USE p_s_parms, ONLY : albsoil, albobs_sw, albobs_vis, albobs_nir, b, sathh, &
                        satcon, smvcst, smvccl, smvcwt, hcap, hcon, sthu

  USE top_pdm, ONLY : fexp, ti_mean, ti_sig, sthzw, zw

  USE trifctl, ONLY : frac_agr, cv

  USE urban_param, ONLY : wrr, hwr, hgt, ztm, disp, albwl, albrd, emisw, emisr

  USE prognostics, ONLY : canopy, cs, gs, snow_tile, t_soil, tstar_tile, lai, &
                          canht_ft, rgrain, rho_snow_grnd, snowdepth,         &
                          snow_grnd, nsnow, sice, sliq, tsnow, rgrainL

  USE snow_param, ONLY : ds

  USE forcing, ONLY : pstar, qw_1, tl_1, lw_down, sw_down, diff_rad, ls_rain, &
                      ls_snow, con_rain, con_snow, u_1, v_1

  USE ozone_vars, ONLY : o3

  IMPLICIT NONE

! Argument types
  INTEGER, INTENT(IN) :: var_id  ! Identifies the variable to fill
  TYPE(data_cube), INTENT(IN), OPTIONAL :: cube
                       ! The data to put in to the variable as a cube
  REAL, INTENT(IN), OPTIONAL :: const_val
                       ! A constant value to fill all elements of the variable
                       ! with


! Work variables
  TYPE(data_cube) :: cube_land  ! Workspace cube for land data
                                ! This is required so that it can be deallocated
                                ! to avoid memory leaks

  REAL :: nsnow_real(land_pts, ntiles)  ! Real version of nsnow


!-----------------------------------------------------------------------------


  IF ( .NOT. PRESENT(cube) .AND. .NOT. PRESENT(const_val) )                   &
    CALL log_fatal("populate_var",                                            &
                   "Neither data cube or const_val have been provided for " //&
                   "variable '" // TRIM(get_string_identifier(var_id)) // "'")

  IF ( PRESENT(cube) .AND. PRESENT(const_val) )                               &
    CALL log_warn("populate_var",                                             &
                  "data cube and const_val both provided for variable '" //   &
                  TRIM(get_string_identifier(var_id)) //                      &
                  "' - using data in preference")


  SELECT CASE ( var_id )
    CASE ( var_id_latitude )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, latitude)
      ELSE
        latitude(:,:) = const_val
      END IF

    CASE ( var_id_longitude )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, longitude)
      ELSE
        longitude(:,:) = const_val
      END IF

    CASE ( var_id_land_fraction )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, flandg)
      ELSE
        flandg(:,:) = const_val
      END IF

    CASE ( var_id_surf_hgt )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, surf_hgt)
      ELSE
        surf_hgt(:,:) = const_val
      END IF

    CASE ( var_id_frac )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, frac)
      ELSE
        frac(:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Soil properties have a corresponding 0-level version for input of constant
! z values. albsoil and albobs_sw etc have no levels anyway
!-----------------------------------------------------------------------------
    CASE ( var_id_albsoil )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, albsoil)
      ELSE
        albsoil(:) = const_val
      END IF

    CASE ( var_id_b )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, b(:,1:sm_levels))
      ELSE
        b(:,:) = const_val
      END IF

    CASE ( var_id_b_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, b(:,1))
        b(:,1:sm_levels) = SPREAD(b(:,1), 2, sm_levels)
      ELSE
        b(:,:) = const_val
      END IF

    CASE ( var_id_sathh )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, sathh(:,1:sm_levels))
      ELSE
        sathh(:,:) = const_val
      END IF

    CASE ( var_id_sathh_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, sathh(:,1))
        sathh(:,1:sm_levels) = SPREAD(sathh(:,1), 2, sm_levels)
      ELSE
        sathh(:,:) = const_val
      END IF

    CASE ( var_id_satcon )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, satcon(:,1:sm_levels))
      ELSE
        satcon(:,:) = const_val
      END IF

    CASE ( var_id_satcon_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, satcon(:,1))
        satcon(:,1:sm_levels) = SPREAD(satcon(:,1), 2, sm_levels)
      ELSE
        satcon(:,:) = const_val
      END IF

    CASE ( var_id_sm_sat )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, smvcst(:,1:sm_levels))
      ELSE
        smvcst(:,:) = const_val
      END IF

    CASE ( var_id_sm_sat_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, smvcst(:,1))
        smvcst(:,1:sm_levels) = SPREAD(smvcst(:,1), 2, sm_levels)
      ELSE
        smvcst(:,:) = const_val
      END IF

    CASE ( var_id_sm_crit )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, smvccl(:,1:sm_levels))
      ELSE
        smvccl(:,:) = const_val
      END IF

    CASE ( var_id_sm_crit_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, smvccl(:,1))
        smvccl(:,1:sm_levels) = SPREAD(smvccl(:,1), 2, sm_levels)
      ELSE
        smvccl(:,:) = const_val
      END IF

    CASE ( var_id_sm_wilt )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, smvcwt(:,1:sm_levels))
      ELSE
        smvcwt(:,:) = const_val
      END IF

    CASE ( var_id_sm_wilt_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, smvcwt(:,1))
        smvcwt(:,1:sm_levels) = SPREAD(smvcwt(:,1), 2, sm_levels)
      ELSE
        smvcwt(:,:) = const_val
      END IF

    CASE ( var_id_hcap )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, hcap(:,1:sm_levels))
      ELSE
        hcap(:,:) = const_val
      END IF

    CASE ( var_id_hcap_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, hcap(:,1))
        hcap(:,1:sm_levels) = SPREAD(hcap(:,1), 2, sm_levels)
      ELSE
        hcap(:,:) = const_val
      END IF

    CASE ( var_id_hcon )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, hcon(:,1:sm_levels))
      ELSE
        hcon(:,:) = const_val
      END IF

    CASE ( var_id_hcon_const_z )
      IF ( PRESENT(cube) ) THEN
! Copy the values from the first vertical level to all others
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, hcon(:,1))
        hcon(:,1:sm_levels) = SPREAD(hcon(:,1), 2, sm_levels)
      ELSE
        hcon(:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! TOPMODEL variables are just on land points
!-----------------------------------------------------------------------------
    CASE ( var_id_fexp )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, fexp)
      ELSE
        fexp(:) = const_val
      END IF

    CASE ( var_id_ti_mean )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, ti_mean)
      ELSE
        ti_mean(:) = const_val
      END IF

    CASE ( var_id_ti_sig )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, ti_sig)
      ELSE
        ti_sig(:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Agricultural fraction is on land points only
!-----------------------------------------------------------------------------
    CASE ( var_id_frac_agr )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, frac_agr)
      ELSE
        frac_agr(:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Urban variables are on land points only
!-----------------------------------------------------------------------------
    CASE ( var_id_wrr )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, wrr)
      ELSE
        wrr(:) = const_val
      END IF

    CASE ( var_id_hwr )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, hwr)
      ELSE
        hwr(:) = const_val
      END IF

    CASE ( var_id_hgt )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, hgt)
      ELSE
        hgt(:) = const_val
      END IF

    CASE ( var_id_ztm )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, ztm)
      ELSE
        ztm(:) = const_val
      END IF

    CASE ( var_id_disp )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, disp)
      ELSE
        disp(:) = const_val
      END IF

    CASE ( var_id_albwl )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, albwl)
      ELSE
        albwl(:) = const_val
      END IF

    CASE ( var_id_albrd )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, albrd)
      ELSE
        albrd(:) = const_val
      END IF

    CASE ( var_id_emisw )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, emisw)
      ELSE
        emisw(:) = const_val
      END IF

    CASE ( var_id_emisr )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, emisr)
      ELSE
        emisr(:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Variables that are set as initial conditions
! These are all land points only, but with varying numbers of levels
!-----------------------------------------------------------------------------
    CASE ( var_id_canopy )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, canopy)
      ELSE
        canopy(:,:) = const_val
      END IF

    CASE ( var_id_cs )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, cs)
      ELSE
        cs(:,:) = const_val
      END IF

    CASE ( var_id_gs )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, gs)
      ELSE
        gs(:) = const_val
      END IF

    CASE ( var_id_snow_tile )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, snow_tile)
      ELSE
        snow_tile(:,:) = const_val
      END IF

    CASE ( var_id_sthuf )
! sthu is used as storage for total wetness
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, sthu)
      ELSE
        sthu(:,:) = const_val
      END IF

    CASE ( var_id_t_soil )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, t_soil)
      ELSE
        t_soil(:,:) = const_val
      END IF

    CASE ( var_id_tstar_tile )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, tstar_tile)
      ELSE
        tstar_tile(:,:) = const_val
      END IF

    CASE ( var_id_lai )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, lai)
      ELSE
        lai(:,:) = const_val
      END IF

    CASE ( var_id_canht )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, canht_ft)
      ELSE
        canht_ft(:,:) = const_val
      END IF

    CASE ( var_id_sthzw )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, sthzw)
      ELSE
        sthzw(:) = const_val
      END IF

    CASE ( var_id_zw )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, zw)
      ELSE
        zw(:) = const_val
      END IF

    CASE ( var_id_rgrain )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, rgrain)
      ELSE
        rgrain(:,:) = const_val
      END IF

    CASE ( var_id_cv )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, cv)
      ELSE
        cv(:) = const_val
      END IF

    CASE ( var_id_rho_snow )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, rho_snow_grnd)
      ELSE
        rho_snow_grnd(:,:) = const_val
      END IF

    CASE ( var_id_snow_depth )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, snowDepth)
      ELSE
        snowDepth(:,:) = const_val
      END IF

    CASE ( var_id_snow_grnd )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, snow_grnd)
      ELSE
        snow_grnd(:,:) = const_val
      END IF

    CASE ( var_id_nsnow )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, nsnow_real)
        nsnow(:,:) = NINT(nsnow_real)
      ELSE
        nsnow(:,:) = NINT(const_val)
      END IF

    CASE ( var_id_snow_ds )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, ds)
      ELSE
        ds(:,:,:) = const_val
      END IF

    CASE ( var_id_snow_ice )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, sice)
      ELSE
        sice(:,:,:) = const_val
      END IF

    CASE ( var_id_snow_liq )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, sliq)
      ELSE
        sliq(:,:,:) = const_val
      END IF

    CASE ( var_id_tsnow )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, tsnow)
      ELSE
        tsnow(:,:,:) = const_val
      END IF

    CASE ( var_id_rgrainl )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, rgrainL)
      ELSE
        rgrainL(:,:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Forcing variables
!-----------------------------------------------------------------------------
    CASE ( var_id_pstar )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, pstar)
      ELSE
        pstar(:,:) = const_val
      END IF

    CASE ( var_id_q )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, qw_1)
      ELSE
        qw_1(:,:) = const_val
      END IF
      qw_1(:,:) = MAX(qw_1, 0.0)

    CASE ( var_id_t )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, tl_1)
      ELSE
        tl_1(:,:) = const_val
      END IF

    CASE ( var_id_rad_net )
! Net downward radiation is stored in lw_down until it is processed
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, lw_down)
      ELSE
        lw_down(:,:) = const_val
      END IF

    CASE ( var_id_lw_net )
! Net LW downward radiation is stored in lw_down until it is processed
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, lw_down)
      ELSE
        lw_down(:,:) = const_val
      END IF

    CASE ( var_id_sw_net )
! Net SW downward radiation is stored in sw_down until it is processed
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, sw_down)
      ELSE
        sw_down(:,:) = const_val
      END IF

    CASE ( var_id_lw_down )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, lw_down)
      ELSE
        lw_down(:,:) = const_val
      END IF

    CASE ( var_id_sw_down )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, sw_down)
      ELSE
        sw_down(:,:) = const_val
      END IF

    CASE ( var_id_diff_rad )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, diff_rad)
      ELSE
        diff_rad(:,:) = const_val
      END IF

    CASE ( var_id_precip )
! Store total precip as large-scale rainfall until it is partitioned
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, ls_rain)
      ELSE
        ls_rain(:,:) = const_val
      END IF

    CASE ( var_id_tot_rain )
! Store total rainfall as large-scale until it is partitioned
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, ls_rain)
      ELSE
        ls_rain(:,:) = const_val
      END IF

    CASE ( var_id_tot_snow )
! If given total snow, we assume it is all large-scale
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, ls_snow)
      ELSE
        ls_snow(:,:) = const_val
      END IF

    CASE ( var_id_con_rain )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, con_rain)
      ELSE
        con_rain(:,:) = const_val
      END IF

    CASE ( var_id_ls_rain )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, ls_rain)
      ELSE
        ls_rain(:,:) = const_val
      END IF

    CASE ( var_id_con_snow )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, con_snow)
      ELSE
        con_snow(:,:) = const_val
      END IF

    CASE ( var_id_ls_snow )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, ls_snow)
      ELSE
        ls_snow(:,:) = const_val
      END IF

    CASE ( var_id_wind )
! Wind speed just goes directly into u component
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, u_1)
      ELSE
        u_1(:,:) = const_val
      END IF

    CASE ( var_id_u )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, u_1)
      ELSE
        u_1(:,:) = const_val
      END IF

    CASE ( var_id_v )
      IF ( PRESENT(cube) ) THEN
        CALL cube_get_data(cube, v_1)
      ELSE
        v_1(:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Other variables that might be prescribed
!-----------------------------------------------------------------------------
    CASE ( var_id_ozone )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, o3)
      ELSE
        o3(:) = const_val
      END IF

    CASE ( var_id_albobs_sw )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, albobs_sw)
      ELSE
        albobs_sw(:) = const_val
      END IF

    CASE ( var_id_albobs_vis )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, albobs_vis)
      ELSE
        albobs_vis(:) = const_val
      END IF

    CASE ( var_id_albobs_nir )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, albobs_nir)
      ELSE
        albobs_nir(:) = const_val
      END IF

    CASE DEFAULT
      CALL log_fatal("populate_var",                                          &
                     "Unrecognised variable for input - '" //                 &
                     TRIM(get_string_identifier(var_id)) // "'. " //          &
                     "See docs for available variables")
  END SELECT


! Free the land cube
  CALL cube_free(cube_land)


  RETURN

END SUBROUTINE populate_var
#endif
