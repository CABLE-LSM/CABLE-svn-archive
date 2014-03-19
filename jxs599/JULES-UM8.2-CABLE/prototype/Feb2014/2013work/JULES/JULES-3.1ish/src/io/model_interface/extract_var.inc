#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/extract_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-09-27 09:33:34 +0100 (Thu, 27 Sep 2012) $
!!   $LastChangedRevision: 609 $
!!
!!****************************************************************************

SUBROUTINE extract_var(var_id, data)

  USE io_constants, ONLY : MDI

  USE string_utils_mod, ONLY : to_string

  USE max_dimensions, ONLY : ntiles_max

  USE water_constants_mod, ONLY : RHO_WATER, TM

  USE csigma, ONLY : sbcon

  USE nstypes, ONLY : ntype, npft, lake

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE switches, ONLY : l_aggregate, l_point_data

  USE latlon_mod, ONLY : latitude, longitude

  USE ancil_info, ONLY : frac, land_pts, land_index, lice_index, ntiles,      &
                         sm_levels, soil_pts, soil_index, tile_index, tile_pts

  USE p_s_parms, ONLY : b, sathh, satcon, smvcst, smvccl, smvcwt, hcap, hcon, &
                        z0_tile, sthu, sthf, catch, cosz

  USE fluxes, ONLY : alb_tile, fqw_tile, hf_snow_melt, sub_surf_roff,         &
                     surf_roff, snomlt_sub_htf, snow_melt, tot_tfall, ecan,   &
                     ei, esoil, fqw_1, ftl_1, land_albedo, latent_heat,       &
                     snomlt_surf_htf, surf_ht_flux, taux_1, tauy_1, tstar,    &
                     melt_tile, anthrop_heat, emis_tile, ext, fsmc, ftl_tile, &
                     le_tile, radnet_tile, surf_ht_store, surf_htf_tile,      &
                     ecan_tile, ei_tile, esoil_tile


  USE prognostics, ONLY : canopy_gb, cs, t_soil, gs, sice, sliq, snow_mass,   &
                          snow_tile, snow_grnd, rgrainL, tsnow, canht_ft,     &
                          canopy, gc, lai, rgrain, rho_snow_grnd, smcl,       &
                          snowdepth, tstar_tile, nsnow

  USE trifctl, ONLY : cv, gpp, lit_c_mn, npp, resp_p, resp_s, resp_s_dr_out,  &
                      c_veg, g_leaf, g_leaf_day, g_leaf_dr_out, g_leaf_phen,  &
                      gpp_ft, lit_c, npp_dr_out, npp_ft, resp_p_ft,           &
                      resp_w_dr_out, resp_w_ft, lai_phen

  USE soil_param, ONLY : dzsoil, zsmc

  USE top_pdm, ONLY : drain, fch4_wetl, fsat, fwetl, qbase, qbase_zw,         &
                      dun_roff, sthzw, zw

  USE jules_mod, ONLY : snowdep_surf

  USE forcing, ONLY : con_rain, con_snow, ls_rain, ls_snow, lw_down, pstar,   &
                      qw_1, sw_down, tl_1, u_1, v_1

  USE screen, ONLY : q1p5m, t1p5m, u10m, v10m, q1p5m_tile, t1p5m_tile

  USE snow_param, ONLY : canSnowTile, ds

  USE surf_param, ONLY : diff_frac

  USE ozone_vars, ONLY : flux_o3_ft, fo3_ft

  USE bvoc_vars, ONLY : isoprene, isoprene_ft, terpene, terpene_ft,           &
                        methanol, methanol_ft, acetone, acetone_ft

  IMPLICIT NONE

! Argument types
  INTEGER, INTENT(IN) :: var_id
                       ! Identifies the variable to extract data from
  REAL, INTENT(OUT) :: data(:,:,:)  ! The extracted data
                                    ! The dimensions are model grid x,
                                    !                    model grid y,
                                    !                    levels

! Work variables
  REAL :: workspace(land_pts, 1)  ! Used to map variables from land points
                                  ! to the full grid that don't have a
                                  ! vertical levels dimension

  REAL :: workspace_tile(land_pts, ntiles)
                                  ! Used in calculation of tile variables

  REAL :: ones(land_pts, sm_levels)  ! An array full of ones to aid with
                                     ! calculation of first frozen/unfrozen
                                     ! layer
  INTEGER :: layer(land_pts)  ! Used in calculation of first frozen/unfrozen
                              ! soil layer

  REAL :: sum_frac(land_pts)  ! Used in calculation of fsmc_gb - the sum
                              ! of frac over all pft tiles

  REAL :: snowd(land_pts), z0(land_pts), fsnow(land_pts,ntype)
                              ! Used in calculation of snow_frac_alb

  REAL :: dz  ! Used in calculation of available soil moisture
              ! Size of current soil layer
  REAL :: ztop  ! Used in calculation of available soil moisture
                ! Depth to top of current soil layer

  INTEGER :: i,j,k,l,n  ! Index variables


!-----------------------------------------------------------------------------

! Initialise data to missing data value
  data(:,:,:) = MDI
  workspace(:,:) = MDI
  workspace_tile(:,:) = MDI

  ones(:,:) = 1.0

  SELECT CASE ( var_id )
    CASE ( var_id_latitude )
      data(:,:,1) = latitude(:,:)

    CASE ( var_id_longitude )
      data(:,:,1) = longitude(:,:)

    CASE ( var_id_frac )
      CALL map_from_land(frac, data)

    CASE ( var_id_b )
      CALL map_from_land(b(:,1:sm_levels), data)

    CASE ( var_id_sathh )
      CALL map_from_land(sathh(:,1:sm_levels), data)

    CASE ( var_id_satcon )
      CALL map_from_land(satcon(:,1:sm_levels), data)

    CASE ( var_id_sm_sat )
      CALL map_from_land(smvcst(:,1:sm_levels), data)

    CASE ( var_id_sm_crit )
      CALL map_from_land(smvccl(:,1:sm_levels), data)

    CASE ( var_id_sm_wilt )
      CALL map_from_land(smvcwt(:,1:sm_levels), data)

    CASE ( var_id_hcap )
      CALL map_from_land(hcap(:,1:sm_levels), data)

    CASE ( var_id_hcon )
      CALL map_from_land(hcon(:,1:sm_levels), data)

    CASE ( var_id_albedo_land )
! Calculate the albedo as used in subroutine control when calculating the net
! shortwave on tiles
! Here we take the average of diffuse albedos in VIS and NIR
      workspace(:,1) = tiles_to_gbm(0.5 * (alb_tile(:,:,2) + alb_tile(:,:,4)))
      CALL map_from_land(workspace, data)

    CASE ( var_id_canopy_gb )
      workspace(:,1) = canopy_gb(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_cs_gb )
      workspace(:,1) = SUM(cs(:,:), 2)
      CALL map_from_land(workspace, data)

    CASE ( var_id_cv )
      workspace(:,1) = cv(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_depth_frozen )
! Get frozen depth from surface
! Start by assuming that there are no frozen layers anywhere, and so frozen
! depth is 0
      workspace(:,:) = 0.0
! Get the first unfrozen layer for every land point using MINLOC
! We give an array of ones as the array to take mins from, so that we get the
! minimum index at which the mask is true
      layer(:) = MINLOC(ones, 2, t_soil(:,:) >= tm)
! If the layer found above is 0 for any point, that means no unfrozen layers
! were found and hence the whole column is frozen
      WHERE ( layer(:) == 0 )
        workspace(:,1) = SUM(dzsoil(:))
      END WHERE
! At points where the layer found above is 1, that means no frozen layers so we
! can ignore those points
! Otherwise, interpolate to estimate depth of zero degC isotherm
! Do all points with the same unfrozen layer at once
      DO i = 2,sm_levels
        WHERE ( layer(:) == i )
          workspace(:,1) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                   &
                           (tm - t_soil(:,i-1)) / (t_soil(:,i) - t_soil(:,i-1))
        END WHERE
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_depth_unfrozen )
! Get unfrozen depth from surface
! See depth_frozen for method description, but swapping frozen for unfrozen
! and vica-versa
      workspace(:,:) = 0.0
      layer(:) = MINLOC(ones, 2, t_soil(:,:) < tm)
      WHERE ( layer(:) == 0 )
        workspace(:,1) = SUM(dzsoil(:))
      END WHERE
      DO i = 2,sm_levels
        WHERE ( layer(:) == i )
          workspace(:,1) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                   &
                           (tm - t_soil(:,i-1)) / (t_soil(:,i) - t_soil(:,i-1))
        END WHERE
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_drain )
      workspace(:,1) = drain(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_elake )
      workspace(:,:) = 0.0
      IF ( .NOT. l_aggregate .AND. lake > 0 )                                 &
        workspace(:,1) = fqw_tile(:,lake) * frac(:,lake)
      CALL map_from_land(workspace, data)

    CASE ( var_id_emis_gb )
      workspace(:,1) = tiles_to_gbm(emis_tile)
      CALL map_from_land(workspace, data)

    CASE ( var_id_fch4_wetl )
      workspace(:,1) = fch4_wetl(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_fsat )
      workspace(:,1) = fsat(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_fsmc_gb )
! Calculate gridbox mean over PFTs.
! Calculate the weighted sum over pfts
      workspace(:,1) = SUM(fsmc * frac(:,1:npft), 2)
      sum_frac(:) = SUM(frac(:,1:npft), 2)
! Normalise to the vegetation fraction
      WHERE ( sum_frac > EPSILON(1.0) )
        workspace(:,1) = workspace(:,1) / sum_frac(:)
      ELSEWHERE
! Where there is no veg, set an impossible value
        workspace(:,1) = MDI
      END WHERE
! Copy the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_fwetl )
      workspace(:,1) = fwetl(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_gpp_gb )
      workspace(:,1) = gpp(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_gs )
      workspace(:,1) = gs(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_hf_snow_melt )
      workspace(:,1) = hf_snow_melt(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_land_index )
      workspace(:,1) = REAL(land_index)
      CALL map_from_land(workspace, data)

    CASE ( var_id_lice_index )
      workspace(:,1) = REAL(lice_index)
      CALL map_from_land(workspace, data)

    CASE ( var_id_lit_c_mean )
      workspace(:,1) = lit_c_mn(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_lw_net )
! Calculate gridbox mean upwards longwave
      workspace(:,1) = sbcon * tiles_to_gbm(emis_tile * tstar_tile**4)
! Now get the net flux.
      DO l = 1,land_pts
        j = (land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j - 1) * t_i_length
        workspace(l,1) = lw_down(i,j) - workspace(l,1)
      ENDDO
! Map to full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_lw_up )
      workspace(:,1) = sbcon * tiles_to_gbm(emis_tile * tstar_tile**4)
      CALL map_from_land(workspace, data)

    CASE ( var_id_npp_gb )
      workspace(:,1) = npp(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_qbase )
      workspace(:,1) = qbase(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_qbase_zw )
      workspace(:,1) = qbase_zw(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_rad_net )
      workspace(:,1) = tiles_to_gbm(radnet_tile)
      CALL map_from_land(workspace, data)

    CASE ( var_id_resp_p_gb )
      workspace(:,1) = resp_p(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_resp_s_gb )
      workspace(:,1) = SUM(resp_s(:,:), 2)
      CALL map_from_land(workspace, data)

    CASE ( var_id_resp_s_dr_out )
! HACK: We only output the total respiration for now
      workspace(:,1) = resp_s_dr_out(:,5)
      CALL map_from_land(workspace, data)

    CASE ( var_id_runoff )
      workspace(:,1) = sub_surf_roff(:) + surf_roff(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_sat_excess_roff )
      workspace(:,1) = dun_roff(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_smc_avail_top )
      workspace(:,:) = 0.0
! We maintain the depth of the top of the current layer as we go down through
! the soil column
      ztop = 0.0
      DO k = 1,sm_levels
! If the top of this layer is below where we want to calculate to, we are done
        IF ( ztop >= zsmc ) EXIT

! Calculate the amount of this layer that we want to take into account
        dz = dzsoil(k)
! If the layer contains the calculation limit, update dz
        IF ( ztop + dz > zsmc ) dz = zsmc - ztop

! Add on the contribution for this layer
        DO j = 1,soil_pts
          i = soil_index(j)
          workspace(i,1) = workspace(i,1) + RHO_WATER * dz *                  &
                               MAX(0.0, sthu(i,k) * smvcst(i,k) - smvcwt(i,k))
        END DO
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_smc_avail_tot )
! This is the same as smc_avail_top, but for whole column
      workspace(:,:) = 0.0
      DO k = 1,sm_levels
        DO j = 1,soil_pts
          i = soil_index(j)
          workspace(i,1) = workspace(i,1) + RHO_WATER * dzsoil(k) *           &
                               MAX(0.0, sthu(i,k) * smvcst(i,k) - smvcwt(i,k))
        END DO
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_smc_tot )
      workspace(:,:) = 0.0
      DO k = 1,sm_levels
        DO j = 1,soil_pts
          i = soil_index(j)
          workspace(i,1) = workspace(i,1) + RHO_WATER * dzsoil(k) *           &
                               MAX(0.0, (sthu(i,k) + sthf(i,k)) * smvcst(i,k))
        END DO
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_snomlt_sub_htf )
      workspace(:,1) = snomlt_sub_htf(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_can_gb )
! Only include tiles where canopy snow model is used
      workspace(:,1) = tiles_to_gbm(snow_tile, canSnowTile)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_depth_gb )
      workspace(:,1) = tiles_to_gbm(snowDepth)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_frac )
! Sum frac over tiles with snow.
      workspace(:,:) = 0.0
      IF ( l_aggregate ) THEN
        WHERE ( snow_tile(:,1) + snow_grnd(:,1) > EPSILON(1.0) )
          workspace(:,1) = 1.0
        END WHERE
      ELSE
        workspace(:,1) = SUM(frac, 2, snow_tile + snow_grnd > EPSILON(1.0))
      END IF
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_frac_alb )
! This calculation is the same as the one in tile_albedo
      snowd(:) = snowdep_surf(:,1)
      z0(:) = z0_tile(:,1)

      fsnow(:,:) = 0.0
      DO n = 1,ntype
        IF ( .NOT. l_aggregate ) THEN
! We only need to update these for each type if l_aggregate=F
          DO j = 1,tile_pts(n)
            l = tile_index(j,n)
            snowd(l) = snowdep_surf(l,n)
            z0(l) = z0_tile(l,n)
          END DO
        END IF
! Calculate snow albedo weighting factor.
        IF ( l_point_data ) THEN
          DO j = 1,tile_pts(n)
            l = tile_index(j,n)
            IF ( snowd(l) > 0.0 ) fsnow(l,n) = 1.0 - EXP(-50.0 * snowd(l))
          END DO
        ELSE
          DO j = 1,tile_pts(n)
            l = tile_index(j,n)
            IF ( snowd(l) > 0.0 ) fsnow(l,n) = snowd(l) / (snowd(l) + 10.0 * z0(l))
          END DO
        END IF
      END DO

      workspace(:,1) = tiles_to_gbm(fsnow)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_grnd_gb )
! Only include tiles where canopy snow model is used
      workspace(:,1) = tiles_to_gbm(snow_grnd, canSnowTile)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_ice_gb )
! Calculate sum of sice along the snow layers dimension
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sice(i,n,1:nsnow(i,n)))
        END DO
      END DO
      workspace(:,1) = tiles_to_gbm(workspace_tile)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_liq_gb )
! Calculate sum of sliq along the snow layers dimension
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sliq(i,n,1:nsnow(i,n)))
        END DO
      END DO
      workspace(:,1) = tiles_to_gbm(workspace_tile)
      CALL map_from_land(workspace, data)

    CASE ( var_id_snow_melt_gb )
      workspace(:,1) = snow_melt(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_soil_index )
      DO l = 1,soil_pts
        workspace(soil_index(l),1) = REAL(soil_index(l))
      ENDDO
      CALL map_from_land(workspace, data)

    CASE ( var_id_sthzw )
      workspace(:,1) = sthzw(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_sub_surf_roff )
      workspace(:,1) = sub_surf_roff(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_surf_roff )
      workspace(:,1) = surf_roff(:)
      CALL map_from_land(workspace, data)

!    CASE ( var_id_surf_roff_inf )
!      workspace(:,1) = roffInfDiag(:)
!      CALL map_from_land(workspace, data)

    CASE ( var_id_swet_liq_tot )
! Divide column unfrozen moisture content by saturated moisture content
      workspace(:,:) = 0.0
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace(i,1) = SUM(dzsoil(:) * MAX(0.0, sthu(i,:) * smvcst(i,:)))   &
                       / SUM(dzsoil(:) * MAX(0.0, smvcst(i,:)))
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_swet_tot )
! Divide column total moisture content by saturated moisture content
      workspace(:,:) = 0.0
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace(i,1) = SUM(dzsoil(:) *                                      &
                             MAX(0.0, (sthu(i,:) + sthf(i,:)) * smvcst(i,:))) &
                       / SUM(dzsoil(:) * MAX(0.0, smvcst(i,:)))
      END DO
! Map the result onto the full grid
      CALL map_from_land(workspace, data)

    CASE ( var_id_sw_net )
! Calculate the albedo as used in subroutine control when calculating
! the net shortwave on tiles
! Here we take the average of diffuse albedos in VIS and NIR.
      DO l = 1,land_pts
        j = (land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j - 1) * t_i_length
        workspace(l,1) = (1.0 -                                               &
                          0.5 * (land_albedo(i,j,2) + land_albedo(i,j,4)))    &
                       * sw_down(i,j)
      ENDDO
      CALL map_from_land(workspace, data)

    CASE ( var_id_tfall )
      workspace(:,1) = tot_tfall(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_trad )
! Assuming emissivity=1.
      workspace(:,1) = ( tiles_to_gbm(tstar_tile**4) )**0.25
      CALL map_from_land(workspace, data)

!    CASE ( var_id_w_flux_sfc )
!      workspace(:,1) = wfluxSfcDiag(:)
!      CALL map_from_land(workspace, data)

    CASE ( var_id_zw )
      workspace(:,1) = zw(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_c_veg )
      CALL map_from_land(c_veg, data)

    CASE ( var_id_canht )
      CALL map_from_land(canht_ft, data)

!    CASE ( var_id_ci )
!      CALL map_from_land(ciDiag, data)

    CASE ( var_id_flux_o3_stom )
      CALL map_from_land(flux_o3_ft, data)

    CASE ( var_id_fsmc )
      CALL map_from_land(fsmc, data)

    CASE ( var_id_g_leaf )
      CALL map_from_land(g_leaf, data)

    CASE ( var_id_g_leaf_day )
      CALL map_from_land(g_leaf_day, data)

    CASE ( var_id_g_leaf_dr_out )
      CALL map_from_land(g_leaf_dr_out, data)

    CASE ( var_id_g_leaf_phen )
      CALL map_from_land(g_leaf_phen, data)

!    CASE ( var_id_g_stom )
!      CALL map_from_land(gstomDiag, data)

    CASE ( var_id_gpp )
      CALL map_from_land(gpp_ft, data)

    CASE ( var_id_lai )
      CALL map_from_land(lai, data)

    CASE ( var_id_lai_phen )
      CALL map_from_land(lai_phen, data)

    CASE ( var_id_lit_c )
      CALL map_from_land(lit_c, data)

    CASE ( var_id_npp_dr_out )
      CALL map_from_land(npp_dr_out, data)

    CASE ( var_id_npp )
      CALL map_from_land(npp_ft, data)

    CASE ( var_id_o3_exp_fac )
      CALL map_from_land(fo3_ft, data)

!    CASE ( var_id_rdc )
!      CALL map_from_land(rdcDiag, data)

    CASE ( var_id_resp_p )
      CALL map_from_land(resp_p_ft, data)

    CASE ( var_id_resp_w_dr_out )
      CALL map_from_land(resp_w_dr_out, data)

    CASE ( var_id_resp_w )
      CALL map_from_land(resp_w_ft, data)

    CASE ( var_id_cs )
      CALL map_from_land(cs, data)

    CASE ( var_id_resp_s )
      CALL map_from_land(resp_s, data)

    CASE ( var_id_con_rain )
      data(:,:,1) = con_rain(:,:)

    CASE ( var_id_con_snow )
      data(:,:,1) = con_snow(:,:)

    CASE ( var_id_cosz )
      data(:,:,1) = RESHAPE(cosz, (/ t_i_length, t_j_length /))

    CASE ( var_id_diff_frac )
      data(:,:,1) = RESHAPE(diff_frac, (/ t_i_length, t_j_length /))

    CASE ( var_id_ecan_gb )
      data(:,:,1) = ecan(:,:)

    CASE ( var_id_ei_gb )
      data(:,:,1) = ei(:,:)

    CASE ( var_id_esoil_gb )
      data(:,:,1) = esoil(:,:)

    CASE ( var_id_fqw_gb )
      data(:,:,1) = fqw_1(:,:)

    CASE ( var_id_ftl_gb )
      data(:,:,1) = ftl_1(:,:)

    CASE ( var_id_land_albedo_1 )
      data(:,:,1) = land_albedo(:,:,1)

    CASE ( var_id_land_albedo_2 )
      data(:,:,1) = land_albedo(:,:,2)

    CASE ( var_id_land_albedo_3 )
      data(:,:,1) = land_albedo(:,:,3)

    CASE ( var_id_land_albedo_4 )
      data(:,:,1) = land_albedo(:,:,4)

    CASE ( var_id_latent_heat )
      data(:,:,1) = latent_heat(:,:)

    CASE ( var_id_ls_rain )
      data(:,:,1) = ls_rain(:,:)

    CASE ( var_id_ls_snow )
      data(:,:,1) = ls_snow(:,:)

    CASE ( var_id_lw_down )
      data(:,:,1) = lw_down(:,:)

    CASE ( var_id_precip )
      data(:,:,1) = ls_rain(:,:) + con_rain(:,:) + ls_snow(:,:) + con_snow(:,:)

    CASE ( var_id_pstar )
      data(:,:,1) = pstar(:,:)

    CASE ( var_id_q1p5m_gb )
      data(:,:,1) = q1p5m(:,:)

    CASE ( var_id_qw1 )
      data(:,:,1) = qw_1(:,:)

    CASE ( var_id_rainfall )
      data(:,:,1) = ls_rain(:,:) + con_rain(:,:)

    CASE ( var_id_snomlt_surf_htf )
      data(:,:,1) = snomlt_surf_htf(:,:)

    CASE ( var_id_snowfall )
      data(:,:,1) = ls_snow(:,:) + con_snow(:,:)

    CASE ( var_id_snow_mass_gb )
      data(:,:,1) = snow_mass(:,:)

    CASE ( var_id_surf_ht_flux_gb )
      data(:,:,1) = surf_ht_flux(:,:)

    CASE ( var_id_sw_down )
      data(:,:,1) = sw_down(:,:)

    CASE ( var_id_t1p5m_gb )
      data(:,:,1) = t1p5m(:,:)

    CASE ( var_id_taux1 )
      data(:,:,1) = taux_1(:,:)

    CASE ( var_id_tauy1 )
      data(:,:,1) = tauy_1(:,:)

    CASE ( var_id_tl1 )
      data(:,:,1) = tl_1(:,:)

    CASE ( var_id_tstar_gb )
      data(:,:,1) = tstar(:,:)

    CASE ( var_id_u1 )
      data(:,:,1) = u_1(:,:)

    CASE ( var_id_u10m )
      data(:,:,1) = u10m(:,:)

    CASE ( var_id_v1 )
      data(:,:,1) = v_1(:,:)

    CASE ( var_id_v10m )
      data(:,:,1) = v10m(:,:)

    CASE ( var_id_wind )
      data(:,:,1) = SQRT( u_1**2 + v_1**2 )

    CASE ( var_id_ext )
      CALL map_from_land(ext, data)

    CASE ( var_id_smcl )
      CALL map_from_land(smcl, data)

    CASE ( var_id_soil_wet )
      CALL map_from_land(sthu + sthf, data)

    CASE ( var_id_sthf )
      CALL map_from_land(sthf, data)

    CASE ( var_id_sthu )
      CALL map_from_land(sthu, data)

    CASE ( var_id_t_soil )
      CALL map_from_land(t_soil, data)

!    CASE ( var_id_w_flux )
!      CALL map_from_land(wfluxDiag, data)

    CASE ( var_id_alb_tile_1 )
      CALL map_from_land(alb_tile(:,:,1), data)

    CASE ( var_id_alb_tile_2 )
      CALL map_from_land(alb_tile(:,:,2), data)

    CASE ( var_id_alb_tile_3 )
      CALL map_from_land(alb_tile(:,:,3), data)

    CASE ( var_id_alb_tile_4 )
      CALL map_from_land(alb_tile(:,:,4), data)

    CASE ( var_id_anthrop_heat )
      CALL map_from_land(anthrop_heat, data)

    CASE ( var_id_canopy )
      CALL map_from_land(canopy, data)

    CASE ( var_id_catch )
      CALL map_from_land(catch, data)

    CASE ( var_id_ecan )
      CALL map_from_land(ecan_tile, data)

    CASE ( var_id_ei )
      CALL map_from_land(ei_tile, data)

    CASE ( var_id_emis )
      CALL map_from_land(emis_tile, data)

    CASE ( var_id_esoil )
      CALL map_from_land(esoil_tile, data)

    CASE ( var_id_fqw )
! Note that fqw_tile does not (always) do this job!
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          l = tile_index(j,n)
          workspace_tile(l,n) = ecan_tile(l,n) + ei_tile(l,n) + esoil_tile(l,n)
! Add lake evaporation
          IF ( n == lake )                                                    &
            workspace_tile(l,n) = workspace_tile(l,n) + fqw_tile(l,n)
        ENDDO
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_ftl )
      CALL map_from_land(ftl_tile, data)

    CASE ( var_id_gc )
      CALL map_from_land(gc, data)

    CASE ( var_id_le )
      CALL map_from_land(le_tile, data)

    CASE ( var_id_nsnow )
      CALL map_from_land(REAL(nsnow), data)

    CASE ( var_id_q1p5m )
      CALL map_from_land(q1p5m_tile, data)

    CASE ( var_id_rad_net_tile )
      CALL map_from_land(radnet_tile, data)

    CASE ( var_id_rgrain )
      CALL map_from_land(rgrain, data)

    CASE ( var_id_snow_can_melt )
! Only include tiles where canopy snow model is used
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) workspace_tile(:,i) = melt_tile(:,i)
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_can )
! Only include tiles where canopy snow model is used
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) workspace_tile(:,i) = snow_tile(:,i)
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_depth )
      CALL map_from_land(snowDepth, data)

!    CASE ( var_id_snow_grnd_can_melt )
!      CALL map_from_land(snowGMeltDiag, data)

    CASE ( var_id_snow_grnd_rho )
      CALL map_from_land(rho_snow_grnd, data)

    CASE ( var_id_snow_grnd )
! Only include tiles where canopy snow model is used
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) workspace_tile(:,i) = snow_grnd(:,i)
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_ground )
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) THEN
          workspace_tile(:,i) = snow_grnd(:,i)
        ELSE
          workspace_tile(:,i) = snow_tile(:,i)
        ENDIF
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_ice_tile )
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sice(i,n,1:nsnow(i,n)))
        ENDDO
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_liq_tile )
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sliq(i,n,1:nsnow(i,n)))
        ENDDO
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_mass )
      workspace_tile(:,:) = snow_tile(:,:)
! Add snow below canopy
      DO n = 1,ntiles
        IF ( canSnowTile(n) )                                                 &
          workspace_tile(:,n) = workspace_tile(:,n) + snow_grnd(:,n)
      ENDDO
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_snow_melt )
      workspace_tile(:,:) = melt_tile(:,:)
! Add melting of snow below canopy.
!     IF ( can_model == 4 ) workspace(:,:) = workspace(:,:) + snowGMeltDiag(:,:)
      CALL map_from_land(workspace_tile, data)

    CASE ( var_id_surf_ht_flux )
      CALL map_from_land(surf_htf_tile, data)

    CASE ( var_id_surf_ht_store )
      CALL map_from_land(surf_ht_store, data)

    CASE ( var_id_t1p5m )
      CALL map_from_land(t1p5m_tile, data)

    CASE ( var_id_tstar )
      CALL map_from_land(tstar_tile, data)

    CASE ( var_id_z0 )
      CALL map_from_land(z0_tile, data)

    CASE ( var_id_tile_index )
      CALL map_from_land(REAL(tile_index), data)

    CASE ( var_id_isoprene_gb )
      workspace(:,1) = isoprene(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_isoprene )
      CALL map_from_land(isoprene_ft, data)

    CASE ( var_id_terpene_gb )
      workspace(:,1) = terpene(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_terpene )
      CALL map_from_land(terpene_ft, data)

    CASE ( var_id_methanol_gb )
      workspace(:,1) = methanol(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_methanol )
      CALL map_from_land(methanol_ft, data)

    CASE ( var_id_acetone_gb )
      workspace(:,1) = acetone(:)
      CALL map_from_land(workspace, data)

    CASE ( var_id_acetone )
      CALL map_from_land(acetone_ft, data)

!-----------------------------------------------------------------------------
! Snow layer variables can be identified using ranges of ids
!-----------------------------------------------------------------------------
    CASE ( var_id_rgrainl:(var_id_rgrainl + ntiles_max - 1) )
      CALL map_from_land(rgrainL(:,var_id - var_id_rgrainl + 1,:), data)

    CASE ( var_id_snow_ds:(var_id_snow_ds + ntiles_max - 1) )
      CALL map_from_land(ds(:,var_id - var_id_snow_ds + 1,:), data)

    CASE ( var_id_snow_ice:(var_id_snow_ice + ntiles_max - 1) )
      CALL map_from_land(sice(:,var_id - var_id_snow_ice + 1,:), data)

    CASE ( var_id_snow_liq:(var_id_snow_liq + ntiles_max - 1) )
      CALL map_from_land(sliq(:,var_id - var_id_snow_liq + 1,:), data)

    CASE ( var_id_tsnow:(var_id_tsnow + ntiles_max - 1) )
      CALL map_from_land(tsnow(:,var_id - var_id_tsnow + 1,:), data)

    CASE DEFAULT
      CALL log_fatal("extract_var",                                           &
                     "Unrecognised variable for output - '" //                &
                     TRIM(get_string_identifier(var_id)) // "'. " //          &
                     "See docs for available variables")
  END SELECT

  RETURN

END SUBROUTINE extract_var
#endif
