#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/populate_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE populate_var(var_id, data, const_val)

  USE string_utils_mod, ONLY : to_string

  USE ancil_info, ONLY : land_pts, sm_levels, frac

  USE latlon_mod, ONLY : latitude, longitude

  USE coastal, ONLY : flandg

  USE c_elevate, ONLY : surf_hgt

  USE p_s_parms, ONLY : albsoil, b, sathh, satcon, smvcst, smvccl, smvcwt,    &
                        hcap, hcon, sthu

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
  REAL, INTENT(IN), OPTIONAL :: data(:,:,:)
                       ! The data to put in to the variable
                       ! The dimensions are x, y and levels
  REAL, INTENT(IN), OPTIONAL :: const_val
                       ! A constant value to fill all elements of the variable
                       ! with


! Work variables
  REAL, ALLOCATABLE :: data_land(:,:)
                                   ! The data mapped onto land points


!-----------------------------------------------------------------------------

  IF ( .NOT. PRESENT(data) .AND. .NOT. PRESENT(const_val) )                   &
    CALL log_fatal("populate_var",                                            &
                   "Neither data or const_val have been provided for " //     &
                   "variable '" // TRIM(get_string_identifier(var_id)) // "'")

  IF ( PRESENT(data) .AND. PRESENT(const_val) )                               &
    CALL log_warn("populate_var",                                             &
                  "data and const_val both provided for variable '" //        &
                  TRIM(get_string_identifier(var_id)) //                      &
                  "' - using data in preference")


! Create a version of the data on land points only if possible
! I.e. if the data variable is present and there are some land points defined
! to map to
  IF ( PRESENT(data) .AND. land_pts > 0 ) THEN
    ALLOCATE(data_land(land_pts,SIZE(data,3)))
    CALL map_to_land(data, data_land)
  END IF


  SELECT CASE ( var_id )
    CASE ( var_id_latitude )
      IF ( PRESENT(data) ) THEN
        latitude(:,:) = data(:,:,1)
      ELSE
        latitude(:,:) = const_val
      END IF

    CASE ( var_id_longitude )
      IF ( PRESENT(data) ) THEN
        longitude(:,:) = data(:,:,1)
      ELSE
        longitude(:,:) = const_val
      END IF

    CASE ( var_id_land_fraction )
      IF ( PRESENT(data) ) THEN
        flandg(:,:) = data(:,:,1)
      ELSE
        flandg(:,:) = const_val
      END IF

    CASE ( var_id_surf_hgt )
      IF ( PRESENT(data) ) THEN
        surf_hgt(:,:) = data_land(:,:)
      ELSE
        surf_hgt(:,:) = const_val
      END IF

    CASE ( var_id_frac )
      IF ( PRESENT(data) ) THEN
        frac(:,:) = data_land(:,:)
      ELSE
        frac(:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Soil properties have a corresponding 0-level version for input of constant
! z values. albsoil has no levels anyway
!-----------------------------------------------------------------------------
    CASE ( var_id_albsoil )
      IF ( PRESENT(data) ) THEN
        albsoil(:) = data_land(:,1)
      ELSE
        albsoil(:) = const_val
      END IF

    CASE ( var_id_b )
      IF ( PRESENT(data) ) THEN
        b(:,1:sm_levels) = data_land(:,:)
      ELSE
        b(:,:) = const_val
      END IF

    CASE ( var_id_b_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        b(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        b(:,:) = const_val
      END IF

    CASE ( var_id_sathh )
      IF ( PRESENT(data) ) THEN
        sathh(:,1:sm_levels) = data_land(:,:)
      ELSE
        sathh(:,:) = const_val
      END IF

    CASE ( var_id_sathh_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        sathh(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        sathh(:,:) = const_val
      END IF

    CASE ( var_id_satcon )
      IF ( PRESENT(data) ) THEN
        satcon(:,1:sm_levels) = data_land(:,:)
      ELSE
        satcon(:,:) = const_val
      END IF

    CASE ( var_id_satcon_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        satcon(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        satcon(:,:) = const_val
      END IF

    CASE ( var_id_sm_sat )
      IF ( PRESENT(data) ) THEN
        smvcst(:,1:sm_levels) = data_land(:,:)
      ELSE
        smvcst(:,:) = const_val
      END IF

    CASE ( var_id_sm_sat_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        smvcst(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        smvcst(:,:) = const_val
      END IF

    CASE ( var_id_sm_crit )
      IF ( PRESENT(data) ) THEN
        smvccl(:,1:sm_levels) = data_land(:,:)
      ELSE
        smvccl(:,:) = const_val
      END IF

    CASE ( var_id_sm_crit_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        smvccl(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        smvccl(:,:) = const_val
      END IF

    CASE ( var_id_sm_wilt )
      IF ( PRESENT(data) ) THEN
        smvcwt(:,1:sm_levels) = data_land(:,:)
      ELSE
        smvcwt(:,:) = const_val
      END IF

    CASE ( var_id_sm_wilt_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        smvcwt(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        smvcwt(:,:) = const_val
      END IF

    CASE ( var_id_hcap )
      IF ( PRESENT(data) ) THEN
        hcap(:,1:sm_levels) = data_land(:,:)
      ELSE
        hcap(:,:) = const_val
      END IF

    CASE ( var_id_hcap_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        hcap(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        hcap(:,:) = const_val
      END IF

    CASE ( var_id_hcon )
      IF ( PRESENT(data) ) THEN
        hcon(:,1:sm_levels) = data_land(:,:)
      ELSE
        hcon(:,:) = const_val
      END IF

    CASE ( var_id_hcon_const_z )
      IF ( PRESENT(data) ) THEN
! Copy the values from the first vertical level to all others
        hcon(:,1:sm_levels) = SPREAD(data_land(:,1), 2, sm_levels)
      ELSE
        hcon(:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! TOPMODEL variables are just on land points
!-----------------------------------------------------------------------------
    CASE ( var_id_fexp )
      IF ( PRESENT(data) ) THEN
        fexp(:) = data_land(:,1)
      ELSE
        fexp(:) = const_val
      END IF

    CASE ( var_id_ti_mean )
      IF ( PRESENT(data) ) THEN
        ti_mean(:) = data_land(:,1)
      ELSE
        ti_mean(:) = const_val
      END IF

    CASE ( var_id_ti_sig )
      IF ( PRESENT(data) ) THEN
        ti_sig(:) = data_land(:,1)
      ELSE
        ti_sig(:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Agricultural fraction is on land points only
!-----------------------------------------------------------------------------
    CASE ( var_id_frac_agr )
      IF ( PRESENT(data) ) THEN
        frac_agr(:) = data_land(:,1)
      ELSE
        frac_agr(:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Urban variables are on land points only
!-----------------------------------------------------------------------------
    CASE ( var_id_wrr )
      IF ( PRESENT(data) ) THEN
        wrr(:) = data_land(:,1)
      ELSE
        wrr(:) = const_val
      END IF

    CASE ( var_id_hwr )
      IF ( PRESENT(data) ) THEN
        hwr(:) = data_land(:,1)
      ELSE
        hwr(:) = const_val
      END IF

    CASE ( var_id_hgt )
      IF ( PRESENT(data) ) THEN
        hgt(:) = data_land(:,1)
      ELSE
        hgt(:) = const_val
      END IF

    CASE ( var_id_ztm )
      IF ( PRESENT(data) ) THEN
        ztm(:) = data_land(:,1)
      ELSE
        ztm(:) = const_val
      END IF

    CASE ( var_id_disp )
      IF ( PRESENT(data) ) THEN
        disp(:) = data_land(:,1)
      ELSE
        disp(:) = const_val
      END IF

    CASE ( var_id_albwl )
      IF ( PRESENT(data) ) THEN
        albwl(:) = data_land(:,1)
      ELSE
        albwl(:) = const_val
      END IF

    CASE ( var_id_albrd )
      IF ( PRESENT(data) ) THEN
        albrd(:) = data_land(:,1)
      ELSE
        albrd(:) = const_val
      END IF

    CASE ( var_id_emisw )
      IF ( PRESENT(data) ) THEN
        emisw(:) = data_land(:,1)
      ELSE
        emisw(:) = const_val
      END IF

    CASE ( var_id_emisr )
      IF ( PRESENT(data) ) THEN
        emisr(:) = data_land(:,1)
      ELSE
        emisr(:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Variables that are set as initial conditions
! These are all land points only, but with varying numbers of levels
!-----------------------------------------------------------------------------
    CASE ( var_id_canopy )
      IF ( PRESENT(data) ) THEN
        canopy(:,:) = data_land(:,:)
      ELSE
        canopy(:,:) = const_val
      END IF

    CASE ( var_id_cs )
      IF ( PRESENT(data) ) THEN
        cs(:,:) = data_land(:,:)
      ELSE
        cs(:,:) = const_val
      END IF

    CASE ( var_id_gs )
      IF ( PRESENT(data) ) THEN
        gs(:) = data_land(:,1)
      ELSE
        gs(:) = const_val
      END IF

    CASE ( var_id_snow_tile )
      IF ( PRESENT(data) ) THEN
        snow_tile(:,:) = data_land(:,:)
      ELSE
        snow_tile(:,:) = const_val
      END IF

    CASE ( var_id_sthuf )
! sthu is used as storage for total wetness
      IF ( PRESENT(data) ) THEN
        sthu(:,:) = data_land(:,:)
      ELSE
        sthu(:,:) = const_val
      END IF

    CASE ( var_id_t_soil )
      IF ( PRESENT(data) ) THEN
        t_soil(:,:) = data_land(:,:)
      ELSE
        t_soil(:,:) = const_val
      END IF

    CASE ( var_id_tstar_tile )
      IF ( PRESENT(data) ) THEN
        tstar_tile(:,:) = data_land(:,:)
      ELSE
        tstar_tile(:,:) = const_val
      END IF

    CASE ( var_id_lai )
      IF ( PRESENT(data) ) THEN
        lai(:,:) = data_land(:,:)
      ELSE
        lai(:,:) = const_val
      END IF

    CASE ( var_id_canht )
      IF ( PRESENT(data) ) THEN
        canht_ft(:,:) = data_land(:,:)
      ELSE
        canht_ft(:,:) = const_val
      END IF

    CASE ( var_id_sthzw )
      IF ( PRESENT(data) ) THEN
        sthzw(:) = data_land(:,1)
      ELSE
        sthzw(:) = const_val
      END IF

    CASE ( var_id_zw )
      IF ( PRESENT(data) ) THEN
        zw(:) = data_land(:,1)
      ELSE
        zw(:) = const_val
      END IF

    CASE ( var_id_rgrain )
      IF ( PRESENT(data) ) THEN
        rgrain(:,:) = data_land(:,:)
      ELSE
        rgrain(:,:) = const_val
      END IF

    CASE ( var_id_cv )
      IF ( PRESENT(data) ) THEN
        cv(:) = data_land(:,1)
      ELSE
        cv(:) = const_val
      END IF

    CASE ( var_id_rho_snow )
      IF ( PRESENT(data) ) THEN
        rho_snow_grnd(:,:) = data_land(:,:)
      ELSE
        rho_snow_grnd(:,:) = const_val
      END IF

    CASE ( var_id_snow_depth )
      IF ( PRESENT(data) ) THEN
        snowDepth(:,:) = data_land(:,:)
      ELSE
        snowDepth(:,:) = const_val
      END IF

    CASE ( var_id_snow_grnd )
      IF ( PRESENT(data) ) THEN
        snow_grnd(:,:) = data_land(:,:)
      ELSE
        snow_grnd(:,:) = const_val
      END IF

    CASE ( var_id_nsnow )
      IF ( PRESENT(data) ) THEN
        nsnow(:,:) = NINT(data_land(:,:))
      ELSE
        nsnow(:,:) = NINT(const_val)
      END IF

! CABLE 
    CASE ( var_id_isnow_flg3l )
      IF ( PRESENT(data) ) THEN
        isnow_flg3l(:,:) = (data_land(:,:))
      ELSE
        isnow_flg3l(:,:) = (const_val)
      END IF

    CASE ( var_id_snow_rho1l )
      IF ( PRESENT(data) ) THEN
        snow_rho1l(:,:) = (data_land(:,:))
      ELSE
        snow_rho1l(:,:) = const_val 
      END IF

    CASE ( var_id_snage_tile )
      IF ( PRESENT(data) ) THEN
        snage_tile(:,:) = (data_land(:,:))
      ELSE
        snage_tile(:,:) = const_val 
      END IF

    CASE ( var_id_snow_rho1 )
      IF ( PRESENT(data) ) THEN
        snow_rho1(:,:) = (data_land(:,:))
      ELSE
        snow_rho1(:,:) = const_val 
      END IF

    CASE ( var_id_snow_rho2 )
      IF ( PRESENT(data) ) THEN
        snow_rho2(:,:) = (data_land(:,:))
      ELSE
        snow_rho2(:,:) = const_val 
      END IF

    CASE ( var_id_snow_rho3 )
      IF ( PRESENT(data) ) THEN
        snow_rho3(:,:) = (data_land(:,:))
      ELSE
        snow_rho3(:,:) = const_val 
      END IF

    CASE ( var_id_snow_mass1 )
      IF ( PRESENT(data) ) THEN
        snow_mass1(:,:) = (data_land(:,:))
      ELSE
        snow_mass1(:,:) = const_val 
      END IF

    CASE ( var_id_snow_mass2 )
      IF ( PRESENT(data) ) THEN
        snow_mass2(:,:) = (data_land(:,:))
      ELSE
        snow_mass2(:,:) = const_val 
      END IF

    CASE ( var_id_snow_mass3 )
      IF ( PRESENT(data) ) THEN
        snow_mass3(:,:) = (data_land(:,:))
      ELSE
        snow_mass3(:,:) = const_val 
      END IF

    CASE ( var_id_snow_depth1 )
      IF ( PRESENT(data) ) THEN
        snow_depth1(:,:) = (data_land(:,:))
      ELSE
        snow_depth1(:,:) = const_val 
      END IF

    CASE ( var_id_snow_depth2 )
      IF ( PRESENT(data) ) THEN
        snow_depth2(:,:) = (data_land(:,:))
      ELSE
        snow_depth2(:,:) = const_val 
      END IF

    CASE ( var_id_snow_depth3 )
      IF ( PRESENT(data) ) THEN
        snow_depth3(:,:) = (data_land(:,:))
      ELSE
        snow_depth3(:,:) = const_val 
      END IF

    CASE ( var_id_snow_tmp1 )
      IF ( PRESENT(data) ) THEN
        snow_tmp1(:,:) = (data_land(:,:))
      ELSE
        snow_tmp1(:,:) = const_val 
      END IF

    CASE ( var_id_snow_tmp2 )
      IF ( PRESENT(data) ) THEN
        snow_tmp2(:,:) = (data_land(:,:))
      ELSE
        snow_tmp2(:,:) = const_val 
      END IF

    CASE ( var_id_snow_tmp3 )
      IF ( PRESENT(data) ) THEN
        snow_tmp3(:,:) = (data_land(:,:))
      ELSE
        snow_tmp3(:,:) = const_val 
      END IF

    CASE ( var_id_sthu_tile1 )
      IF ( PRESENT(data) ) THEN
        sthu_tile1(:,:) = (data_land(:,:))
      ELSE
        sthu_tile1(:,:) = const_val 
      END IF

    CASE ( var_id_sthu_tile2 )
      IF ( PRESENT(data) ) THEN
        sthu_tile2(:,:) = (data_land(:,:))
      ELSE
        sthu_tile2(:,:) = const_val 
      END IF

    CASE ( var_id_sthu_tile3 )
      IF ( PRESENT(data) ) THEN
        sthu_tile3(:,:) = (data_land(:,:))
      ELSE
        sthu_tile3(:,:) = const_val 
      END IF

    CASE ( var_id_sthu_tile4 )
      IF ( PRESENT(data) ) THEN
        sthu_tile4(:,:) = (data_land(:,:))
      ELSE
        sthu_tile4(:,:) = const_val 
      END IF

    CASE ( var_id_sthu_tile5 )
      IF ( PRESENT(data) ) THEN
        sthu_tile5(:,:) = (data_land(:,:))
      ELSE
        sthu_tile5(:,:) = const_val 
      END IF

    CASE ( var_id_sthu_tile6 )
      IF ( PRESENT(data) ) THEN
        sthu_tile6(:,:) = (data_land(:,:))
      ELSE
        sthu_tile6(:,:) = const_val 
      END IF

    CASE ( var_id_sthf_tile1 )
      IF ( PRESENT(data) ) THEN
        sthf_tile1(:,:) = (data_land(:,:))
      ELSE
        sthf_tile1(:,:) = const_val 
      END IF

    CASE ( var_id_sthf_tile2 )
      IF ( PRESENT(data) ) THEN
        sthf_tile2(:,:) = (data_land(:,:))
      ELSE
        sthf_tile2(:,:) = const_val 
      END IF

    CASE ( var_id_sthf_tile3 )
      IF ( PRESENT(data) ) THEN
        sthf_tile3(:,:) = (data_land(:,:))
      ELSE
        sthf_tile3(:,:) = const_val 
      END IF

    CASE ( var_id_sthf_tile4 )
      IF ( PRESENT(data) ) THEN
        sthf_tile4(:,:) = (data_land(:,:))
      ELSE
        sthf_tile4(:,:) = const_val 
      END IF

    CASE ( var_id_sthf_tile5 )
      IF ( PRESENT(data) ) THEN
        sthf_tile5(:,:) = (data_land(:,:))
      ELSE
        sthf_tile5(:,:) = const_val 
      END IF

    CASE ( var_id_sthf_tile6 )
      IF ( PRESENT(data) ) THEN
        sthf_tile6(:,:) = (data_land(:,:))
      ELSE
        sthf_tile6(:,:) = const_val 
      END IF

    CASE ( var_id_smcl_tile1 )
      IF ( PRESENT(data) ) THEN
        smcl_tile1(:,:) = (data_land(:,:))
      ELSE
        smcl_tile1(:,:) = const_val 
      END IF

    CASE ( var_id_smcl_tile2 )
      IF ( PRESENT(data) ) THEN
        smcl_tile2(:,:) = (data_land(:,:))
      ELSE
        smcl_tile2(:,:) = const_val 
      END IF

    CASE ( var_id_smcl_tile3 )
      IF ( PRESENT(data) ) THEN
        smcl_tile3(:,:) = (data_land(:,:))
      ELSE
        smcl_tile3(:,:) = const_val 
      END IF

    CASE ( var_id_smcl_tile4 )
      IF ( PRESENT(data) ) THEN
        smcl_tile4(:,:) = (data_land(:,:))
      ELSE
        smcl_tile4(:,:) = const_val 
      END IF

    CASE ( var_id_smcl_tile5 )
      IF ( PRESENT(data) ) THEN
        smcl_tile5(:,:) = (data_land(:,:))
      ELSE
        smcl_tile5(:,:) = const_val 
      END IF

    CASE ( var_id_smcl_tile6 )
      IF ( PRESENT(data) ) THEN
        smcl_tile6(:,:) = (data_land(:,:))
      ELSE
        smcl_tile6(:,:) = const_val 
      END IF

    CASE ( var_id_tsoil_tile1 )
      IF ( PRESENT(data) ) THEN
        tsoil_tile1(:,:) = (data_land(:,:))
      ELSE
        tsoil_tile1(:,:) = const_val 
      END IF

    CASE ( var_id_tsoil_tile2 )
      IF ( PRESENT(data) ) THEN
        tsoil_tile2(:,:) = (data_land(:,:))
      ELSE
        tsoil_tile2(:,:) = const_val 
      END IF

    CASE ( var_id_tsoil_tile3 )
      IF ( PRESENT(data) ) THEN
        tsoil_tile3(:,:) = (data_land(:,:))
      ELSE
        tsoil_tile3(:,:) = const_val 
      END IF

    CASE ( var_id_tsoil_tile4 )
      IF ( PRESENT(data) ) THEN
        tsoil_tile4(:,:) = (data_land(:,:))
      ELSE
        tsoil_tile4(:,:) = const_val 
      END IF

    CASE ( var_id_tsoil_tile5 )
      IF ( PRESENT(data) ) THEN
        tsoil_tile5(:,:) = (data_land(:,:))
      ELSE
        tsoil_tile5(:,:) = const_val 
      END IF

    CASE ( var_id_tsoil_tile6 )
      IF ( PRESENT(data) ) THEN
        tsoil_tile6(:,:) = (data_land(:,:))
      ELSE
        tsoil_tile6(:,:) = const_val 
      END IF


! Snow layer variables are only available to be set as constants
    CASE ( var_id_snow_ds )
      IF ( PRESENT(const_val) ) THEN
        ds(:,:,:) = const_val
      ELSE
        CALL log_fatal("populate_var",                                        &
                       "'snow_ds' can only be populated using a " //          &
                       "constant value or from a dump file")
      END IF

    CASE ( var_id_snow_ice )
      IF ( PRESENT(const_val) ) THEN
        sice(:,:,:) = const_val
      ELSE
        CALL log_fatal("populate_var",                                        &
                       "'snow_ice' can only be populated using a " //         &
                       "constant value or from a dump file")
      END IF

    CASE ( var_id_snow_liq )
      IF ( PRESENT(const_val) ) THEN
        sliq(:,:,:) = const_val
      ELSE
        CALL log_fatal("populate_var",                                        &
                       "'snow_liq' can only be populated using a " //         &
                       "constant value or from a dump file")
      END IF

    CASE ( var_id_tsnow )
      IF ( PRESENT(const_val) ) THEN
        tsnow(:,:,:) = const_val
      ELSE
        CALL log_fatal("populate_var",                                        &
                       "'tsnow' can only be populated using a " //            &
                       "constant value or from a dump file")
      END IF

    CASE ( var_id_rgrainl )
      IF ( PRESENT(const_val) ) THEN
        rgrainL(:,:,:) = const_val
      ELSE
        CALL log_fatal("populate_var",                                        &
                       "'rgrainl' can only be populated using a " //          &
                       "constant value or from a dump file")
      END IF

!-----------------------------------------------------------------------------
! Forcing variables
!-----------------------------------------------------------------------------
    CASE ( var_id_pstar )
      IF ( PRESENT(data) ) THEN
        pstar(:,:) = data(:,:,1)
      ELSE
        pstar(:,:) = const_val
      END IF

    CASE ( var_id_q )
      IF ( PRESENT(data) ) THEN
        qw_1(:,:) = data(:,:,1)
      ELSE
        qw_1(:,:) = const_val
      END IF
      qw_1(:,:) = MAX(qw_1, 0.0)

    CASE ( var_id_t )
      IF ( PRESENT(data) ) THEN
        tl_1(:,:) = data(:,:,1)
      ELSE
        tl_1(:,:) = const_val
      END IF

    CASE ( var_id_rad_net )
! Net downward radiation is stored in lw_down until it is processed
      IF ( PRESENT(data) ) THEN
        lw_down(:,:) = data(:,:,1)
      ELSE
        lw_down(:,:) = const_val
      END IF

    CASE ( var_id_lw_net )
! Net LW downward radiation is stored in lw_down until it is processed
      IF ( PRESENT(data) ) THEN
        lw_down(:,:) = data(:,:,1)
      ELSE
        lw_down(:,:) = const_val
      END IF

    CASE ( var_id_sw_net )
! Net SW downward radiation is stored in sw_down until it is processed
      IF ( PRESENT(data) ) THEN
        sw_down(:,:) = data(:,:,1)
      ELSE
        sw_down(:,:) = const_val
      END IF

    CASE ( var_id_lw_down )
      IF ( PRESENT(data) ) THEN
        lw_down(:,:) = data(:,:,1)
      ELSE
        lw_down(:,:) = const_val
      END IF

    CASE ( var_id_sw_down )
      IF ( PRESENT(data) ) THEN
        sw_down(:,:) = data(:,:,1)
      ELSE
        sw_down(:,:) = const_val
      END IF

    CASE ( var_id_diff_rad )
      IF ( PRESENT(data) ) THEN
        diff_rad(:,:) = data(:,:,1)
      ELSE
        diff_rad(:,:) = const_val
      END IF

    CASE ( var_id_precip )
! Store total precip as large-scale rainfall until it is partitioned
      IF ( PRESENT(data) ) THEN
        ls_rain(:,:) = data(:,:,1)
      ELSE
        ls_rain(:,:) = const_val
      END IF

    CASE ( var_id_tot_rain )
! Store total rainfall as large-scale until it is partitioned
      IF ( PRESENT(data) ) THEN
        ls_rain(:,:) = data(:,:,1)
      ELSE
        ls_rain(:,:) = const_val
      END IF

    CASE ( var_id_tot_snow )
! If given total snow, we assume it is all large-scale
      IF ( PRESENT(data) ) THEN
        ls_snow(:,:) = data(:,:,1)
      ELSE
        ls_snow(:,:) = const_val
      END IF

    CASE ( var_id_con_rain )
      IF ( PRESENT(data) ) THEN
        con_rain(:,:) = data(:,:,1)
      ELSE
        con_rain(:,:) = const_val
      END IF

    CASE ( var_id_ls_rain )
      IF ( PRESENT(data) ) THEN
        ls_rain(:,:) = data(:,:,1)
      ELSE
        ls_rain(:,:) = const_val
      END IF

    CASE ( var_id_con_snow )
      IF ( PRESENT(data) ) THEN
        con_snow(:,:) = data(:,:,1)
      ELSE
        con_snow(:,:) = const_val
      END IF

    CASE ( var_id_ls_snow )
      IF ( PRESENT(data) ) THEN
        ls_snow(:,:) = data(:,:,1)
      ELSE
        ls_snow(:,:) = const_val
      END IF

    CASE ( var_id_wind )
! Wind speed just goes directly into u component
      IF ( PRESENT(data) ) THEN
        u_1(:,:) = data(:,:,1)
      ELSE
        u_1(:,:) = const_val
      END IF

    CASE ( var_id_u )
      IF ( PRESENT(data) ) THEN
        u_1(:,:) = data(:,:,1)
      ELSE
        u_1(:,:) = const_val
      END IF

    CASE ( var_id_v )
      IF ( PRESENT(data) ) THEN
        v_1(:,:) = data(:,:,1)
      ELSE
        v_1(:,:) = const_val
      END IF

!-----------------------------------------------------------------------------
! Other variables that might be prescribed
!-----------------------------------------------------------------------------
    CASE ( var_id_ozone )
      IF ( PRESENT(data) ) THEN
        o3(:) = data_land(:,1)
      ELSE
        o3(:) = const_val
      END IF

    CASE DEFAULT
      CALL log_fatal("populate_var",                                          &
                     "Unrecognised variable for input - '" //                 &
                     TRIM(get_string_identifier(var_id)) // "'. " //          &
                     "See docs for available variables")
  END SELECT

  IF ( ALLOCATED(data_land) ) DEALLOCATE(data_land)

  RETURN

END SUBROUTINE populate_var
#endif
