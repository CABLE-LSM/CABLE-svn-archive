#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION extract_var(var_id) RESULT(cube)

  USE io_constants, ONLY : MDI

  USE string_utils_mod, ONLY : to_string

  USE data_cube_mod, ONLY : data_cube, cube_from_array, cube_safe_copy,       &
                            cube_free, OPERATOR (*), OPERATOR (-)

  USE max_dimensions, ONLY : ntiles_max

  USE water_constants_mod, ONLY : RHO_WATER, TM

  USE csigma, ONLY : sbcon

  USE nstypes, ONLY : ntype, npft, lake

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE switches, ONLY : l_aggregate, l_point_data

  USE model_grid_mod, ONLY : latitude, longitude

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


  USE prognostics, ONLY : cs, t_soil, gs, sice, sliq, snow_mass,              &
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

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the data currently
!   associated with that variable as a cube
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: var_id
                       ! Identifies the variable to extract data from


! Return
  TYPE(data_cube) :: cube  ! The extracted data


! Work variables
  TYPE(data_cube) :: cube_land  ! Workspace cube for land data
                                ! This is required so that it can be deallocated
                                ! to avoid memory leaks

  TYPE(data_cube) :: upward_cube, emis_cube, downward_cube, lw_down_cube
                                ! Work cubes for use in calculation of lw_net
                                ! Required to avoid memory leaks

  REAL :: workspace_land(land_pts)  ! Used as a space for calculations
                                    ! before creating a cube

  REAL :: workspace_tile(land_pts, ntiles)
                                  ! Used in calculation of tile variables

  REAL :: ones(land_pts, sm_levels)  ! An array full of ones to aid with
                                     ! calculation of first frozen/unfrozen
                                     ! layer
  INTEGER :: layer(land_pts)  ! Used in calculation of first frozen/unfrozen
                              ! soil layer

  REAL :: sum_frac(land_pts)  ! Used in calculation of fsmc_gb - the sum
                              ! of frac over all pft tiles

  REAL :: dz  ! Used in calculation of available soil moisture
              ! Size of current soil layer
  REAL :: ztop  ! Used in calculation of available soil moisture
                ! Depth to top of current soil layer

  INTEGER :: i,j,k,l,n  ! Index variables


!-----------------------------------------------------------------------------

! Initialise data to missing data value
  workspace_land(:) = MDI
  workspace_tile(:,:) = MDI

  ones(:,:) = 1.0

  SELECT CASE ( var_id )
    CASE ( var_id_latitude )
      cube = cube_from_array(latitude)

    CASE ( var_id_longitude )
      cube = cube_from_array(longitude)

    CASE ( var_id_frac )
      cube_land = cube_from_array(frac)
      cube = map_from_land(cube_land)

    CASE ( var_id_b )
      cube_land = cube_from_array(b(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_sathh )
      cube_land = cube_from_array(sathh(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_satcon )
      cube_land = cube_from_array(satcon(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_sm_sat )
      cube_land = cube_from_array(smvcst(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_sm_crit )
      cube_land = cube_from_array(smvccl(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_sm_wilt )
      cube_land = cube_from_array(smvcwt(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_hcap )
      cube_land = cube_from_array(hcap(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_hcon )
      cube_land = cube_from_array(hcon(:,1:sm_levels))
      cube = map_from_land(cube_land)

    CASE ( var_id_albedo_land )
! Calculate the albedo as used in subroutine control when calculating the net
! shortwave on tiles
! Here we take the average of diffuse albedos in VIS and NIR
      cube_land = cube_from_array(                                            &
        tiles_to_gbm(0.5 * (alb_tile(:,:,2) + alb_tile(:,:,4)))               &
      )
      cube = map_from_land(cube_land)

    CASE ( var_id_canopy_gb )
!     Don't use the canopy_gb variable, as this is calculated part-way through
!     a timestep.
      cube_land = cube_from_array(tiles_to_gbm( canopy ))
      cube = map_from_land(cube_land)

    CASE ( var_id_cs_gb )
      cube_land = cube_from_array(SUM(cs(:,:), 2))
      cube = map_from_land(cube_land)

    CASE ( var_id_cv )
      cube_land = cube_from_array(cv)
      cube = map_from_land(cube_land)

    CASE ( var_id_depth_frozen )
! Get frozen depth from surface
! Start by assuming that there are no frozen layers anywhere, and so frozen
! depth is 0
      workspace_land(:) = 0.0
! Get the first unfrozen layer for every land point using MINLOC
! We give an array of ones as the array to take mins from, so that we get the
! minimum index at which the mask is true
      layer(:) = MINLOC(ones, 2, t_soil(:,:) >= tm)
! If the layer found above is 0 for any point, that means no unfrozen layers
! were found and hence the whole column is frozen
      WHERE ( layer(:) == 0 )
        workspace_land(:) = SUM(dzsoil(:))
      END WHERE
! At points where the layer found above is 1, that means no frozen layers so we
! can ignore those points
! Otherwise, interpolate to estimate depth of zero degC isotherm
! Do all points with the same unfrozen layer at once
      DO i = 2,sm_levels
        WHERE ( layer(:) == i )
          workspace_land(:) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                &
                           (tm - t_soil(:,i-1)) / (t_soil(:,i) - t_soil(:,i-1))
        END WHERE
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_depth_unfrozen )
! Get unfrozen depth from surface
! See depth_frozen for method description, but swapping frozen for unfrozen
! and vica-versa
      workspace_land(:) = 0.0
      layer(:) = MINLOC(ones, 2, t_soil(:,:) < tm)
      WHERE ( layer(:) == 0 )
        workspace_land(:) = SUM(dzsoil(:))
      END WHERE
      DO i = 2,sm_levels
        WHERE ( layer(:) == i )
          workspace_land(:) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                &
                           (tm - t_soil(:,i-1)) / (t_soil(:,i) - t_soil(:,i-1))
        END WHERE
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_drain )
      cube_land = cube_from_array(drain)
      cube = map_from_land(cube_land)

    CASE ( var_id_elake )
      workspace_land(:) = 0.0
      IF ( .NOT. l_aggregate .AND. lake > 0 )                                 &
        workspace_land(:) = fqw_tile(:,lake) * frac(:,lake)
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_emis_gb )
      cube_land = cube_from_array(tiles_to_gbm(emis_tile))
      cube = map_from_land(cube_land)

    CASE ( var_id_fch4_wetl )
      cube_land = cube_from_array(fch4_wetl)
      cube = map_from_land(cube_land)

    CASE ( var_id_fsat )
      cube_land = cube_from_array(fsat)
      cube = map_from_land(cube_land)

    CASE ( var_id_fsmc_gb )
! Calculate gridbox mean over PFTs.
! Calculate the weighted sum over pfts
      workspace_land(:) = SUM(fsmc * frac(:,1:npft), 2)
      sum_frac(:) = SUM(frac(:,1:npft), 2)
! Normalise to the vegetation fraction
      WHERE ( sum_frac > EPSILON(1.0) )
        workspace_land(:) = workspace_land(:) / sum_frac(:)
      ELSEWHERE
! Where there is no veg, set an impossible value
        workspace_land(:) = MDI
      END WHERE
! Copy the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_fwetl )
      cube_land = cube_from_array(fwetl)
      cube = map_from_land(cube_land)

    CASE ( var_id_gpp_gb )
      cube_land = cube_from_array(gpp)
      cube = map_from_land(cube_land)

    CASE ( var_id_gs )
      cube_land = cube_from_array(gs)
      cube = map_from_land(cube_land)

    CASE ( var_id_hf_snow_melt )
      cube_land = cube_from_array(hf_snow_melt)
      cube = map_from_land(cube_land)

    CASE ( var_id_land_index )
      cube_land = cube_from_array(REAL(land_index))
      cube = map_from_land(cube_land)

    CASE ( var_id_lice_index )
      cube_land = cube_from_array(REAL(lice_index))
      cube = map_from_land(cube_land)

    CASE ( var_id_lit_c_mean )
      cube_land = cube_from_array(lit_c_mn)
      cube = map_from_land(cube_land)

    CASE ( var_id_lw_net )
! Calculate gridbox mean upwards longwave
      cube_land = cube_from_array(                                            &
        sbcon * tiles_to_gbm(emis_tile * tstar_tile**4)                       &
      )
      upward_cube = map_from_land(cube_land)
      CALL cube_free(cube_land)

! Calculate gridbox mean downward longwave
      cube_land = cube_from_array(tiles_to_gbm(emis_tile))
      emis_cube = map_from_land(cube_land)
      lw_down_cube = cube_from_array(lw_down)
      downward_cube = emis_cube * lw_down_cube

! Now calculate the net flux
      cube = downward_cube - upward_cube

! Free work cubes
      CALL cube_free(upward_cube)
      CALL cube_free(emis_cube)
      CALL cube_free(lw_down_cube)
      CALL cube_free(downward_cube)

    CASE ( var_id_lw_up )
      cube_land = cube_from_array(                                            &
        sbcon * tiles_to_gbm(emis_tile * tstar_tile**4)                       &
      )
      cube = map_from_land(cube_land)

    CASE ( var_id_npp_gb )
      cube_land = cube_from_array(npp)
      cube = map_from_land(cube_land)

    CASE ( var_id_qbase )
      cube_land = cube_from_array(qbase)
      cube = map_from_land(cube_land)

    CASE ( var_id_qbase_zw )
      cube_land = cube_from_array(qbase_zw)
      cube = map_from_land(cube_land)

    CASE ( var_id_rad_net )
      cube_land = cube_from_array(tiles_to_gbm(radnet_tile))
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_p_gb )
      cube_land = cube_from_array(resp_p)
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_s_gb )
      cube_land = cube_from_array(SUM(resp_s(:,:), 2))
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_s_dr_out )
! HACK: We only output the total respiration for now
      cube_land = cube_from_array(resp_s_dr_out(:,5))
      cube = map_from_land(cube_land)

    CASE ( var_id_runoff )
      cube_land = cube_from_array(sub_surf_roff(:) + surf_roff(:))
      cube = map_from_land(cube_land)

    CASE ( var_id_sat_excess_roff )
      cube_land = cube_from_array(dun_roff)
      cube = map_from_land(cube_land)

    CASE ( var_id_smc_avail_top )
      workspace_land(:) = 0.0
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
          workspace_land(i) = workspace_land(i) + RHO_WATER * dz *            &
                              MAX(0.0, sthu(i,k) * smvcst(i,k) - smvcwt(i,k))
        END DO
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_smc_avail_tot )
! This is the same as smc_avail_top, but for whole column
      workspace_land(:) = 0.0
      DO k = 1,sm_levels
        DO j = 1,soil_pts
          i = soil_index(j)
          workspace_land(i) = workspace_land(i) + RHO_WATER * dzsoil(k) *     &
                              MAX(0.0, sthu(i,k) * smvcst(i,k) - smvcwt(i,k))
        END DO
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_smc_tot )
      workspace_land(:) = 0.0
      DO k = 1,sm_levels
        DO j = 1,soil_pts
          i = soil_index(j)
          workspace_land(i) = workspace_land(i) + RHO_WATER * dzsoil(k) *     &
                              MAX(0.0, (sthu(i,k) + sthf(i,k)) * smvcst(i,k))
        END DO
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_snomlt_sub_htf )
      cube_land = cube_from_array(snomlt_sub_htf)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_can_gb )
! Only include tiles where canopy snow model is used
      cube_land = cube_from_array(tiles_to_gbm(snow_tile, canSnowTile))
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_depth_gb )
      cube_land = cube_from_array(tiles_to_gbm(snowDepth))
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_frac )
! Sum frac over tiles with snow.
      workspace_land(:) = 0.0
      IF ( l_aggregate ) THEN
        WHERE ( snow_tile(:,1) + snow_grnd(:,1) > EPSILON(1.0) )
          workspace_land(:) = 1.0
        END WHERE
      ELSE
        workspace_land(:) = SUM(frac, 2, snow_tile + snow_grnd > EPSILON(1.0))
      END IF
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_grnd_gb )
! Only include tiles where canopy snow model is used
      cube_land = cube_from_array(tiles_to_gbm(snow_grnd, canSnowTile))
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_ice_gb )
! Calculate sum of sice along the snow layers dimension
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sice(i,n,1:nsnow(i,n)))
        END DO
      END DO
      cube_land = cube_from_array(tiles_to_gbm(workspace_tile))
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_liq_gb )
! Calculate sum of sliq along the snow layers dimension
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sliq(i,n,1:nsnow(i,n)))
        END DO
      END DO
      cube_land = cube_from_array(tiles_to_gbm(workspace_tile))
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_melt_gb )
      cube_land = cube_from_array(snow_melt)
      cube = map_from_land(cube_land)

    CASE ( var_id_soil_index )
      DO l = 1,soil_pts
        workspace_land(soil_index(l)) = REAL(soil_index(l))
      ENDDO
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_sthzw )
      cube_land = cube_from_array(sthzw)
      cube = map_from_land(cube_land)

    CASE ( var_id_sub_surf_roff )
      cube_land = cube_from_array(sub_surf_roff)
      cube = map_from_land(cube_land)

    CASE ( var_id_surf_roff )
      cube_land = cube_from_array(surf_roff)
      cube = map_from_land(cube_land)

    CASE ( var_id_swet_liq_tot )
! Divide column unfrozen moisture content by saturated moisture content
      workspace_land(:) = 0.0
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace_land(i) = SUM(dzsoil(:) * MAX(0.0, sthu(i,:) * smvcst(i,:)))&
                          / SUM(dzsoil(:) * MAX(0.0, smvcst(i,:)))
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_swet_tot )
! Divide column total moisture content by saturated moisture content
      workspace_land(:) = 0.0
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace_land(i) = SUM(dzsoil(:) *                                   &
                             MAX(0.0, (sthu(i,:) + sthf(i,:)) * smvcst(i,:))) &
                          / SUM(dzsoil(:) * MAX(0.0, smvcst(i,:)))
      END DO
! Map the result onto the full grid
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_sw_net )
! Calculate the albedo as used in subroutine control when calculating
! the net shortwave on tiles
! Here we take the average of diffuse albedos in VIS and NIR.
      DO l = 1,land_pts
        j = (land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j - 1) * t_i_length
        workspace_land(l) = (1.0 -                                            &
                             0.5 * (land_albedo(i,j,2) + land_albedo(i,j,4))) &
                          * sw_down(i,j)
      ENDDO
      cube_land = cube_from_array(workspace_land)
      cube = map_from_land(cube_land)

    CASE ( var_id_tfall )
      cube_land = cube_from_array(tot_tfall)
      cube = map_from_land(cube_land)

    CASE ( var_id_trad )
! Assuming emissivity=1.
      cube_land = cube_from_array(( tiles_to_gbm(tstar_tile**4) )**0.25)
      cube = map_from_land(cube_land)

    CASE ( var_id_zw )
      cube_land = cube_from_array(zw)
      cube = map_from_land(cube_land)

    CASE ( var_id_c_veg )
      cube_land = cube_from_array(c_veg)
      cube = map_from_land(cube_land)

    CASE ( var_id_canht )
      cube_land = cube_from_array(canht_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_flux_o3_stom )
      cube_land = cube_from_array(flux_o3_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_fsmc )
      cube_land = cube_from_array(fsmc)
      cube = map_from_land(cube_land)

    CASE ( var_id_g_leaf )
      cube_land = cube_from_array(g_leaf)
      cube = map_from_land(cube_land)

    CASE ( var_id_g_leaf_day )
      cube_land = cube_from_array(g_leaf_day)
      cube = map_from_land(cube_land)

    CASE ( var_id_g_leaf_dr_out )
      cube_land = cube_from_array(g_leaf_dr_out)
      cube = map_from_land(cube_land)

    CASE ( var_id_g_leaf_phen )
      cube_land = cube_from_array(g_leaf_phen)
      cube = map_from_land(cube_land)

    CASE ( var_id_gpp )
      cube_land = cube_from_array(gpp_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_lai )
      cube_land = cube_from_array(lai)
      cube = map_from_land(cube_land)

    CASE ( var_id_lai_phen )
      cube_land = cube_from_array(lai_phen)
      cube = map_from_land(cube_land)

    CASE ( var_id_lit_c )
      cube_land = cube_from_array(lit_c)
      cube = map_from_land(cube_land)

    CASE ( var_id_npp_dr_out )
      cube_land = cube_from_array(npp_dr_out)
      cube = map_from_land(cube_land)

    CASE ( var_id_npp )
      cube_land = cube_from_array(npp_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_o3_exp_fac )
      cube_land = cube_from_array(fo3_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_p )
      cube_land = cube_from_array(resp_p_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_w_dr_out )
      cube_land = cube_from_array(resp_w_dr_out)
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_w )
      cube_land = cube_from_array(resp_w_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_cs )
      cube_land = cube_from_array(cs)
      cube = map_from_land(cube_land)

    CASE ( var_id_resp_s )
      cube_land = cube_from_array(resp_s)
      cube = map_from_land(cube_land)

    CASE ( var_id_con_rain )
      cube = cube_from_array(con_rain)

    CASE ( var_id_con_snow )
      cube = cube_from_array(con_snow)

    CASE ( var_id_cosz )
      cube = cube_from_array(RESHAPE(cosz, (/ t_i_length, t_j_length /)))

    CASE ( var_id_diff_frac )
      cube = cube_from_array(RESHAPE(diff_frac, (/ t_i_length, t_j_length /)))

    CASE ( var_id_ecan_gb )
      cube = cube_from_array(ecan)

    CASE ( var_id_ei_gb )
      cube = cube_from_array(ei)

    CASE ( var_id_esoil_gb )
      cube = cube_from_array(esoil)

    CASE ( var_id_fqw_gb )
      cube = cube_from_array(fqw_1)

    CASE ( var_id_ftl_gb )
      cube = cube_from_array(ftl_1)

    CASE ( var_id_land_albedo_1 )
      cube = cube_from_array(land_albedo(:,:,1))

    CASE ( var_id_land_albedo_2 )
      cube = cube_from_array(land_albedo(:,:,2))

    CASE ( var_id_land_albedo_3 )
      cube = cube_from_array(land_albedo(:,:,3))

    CASE ( var_id_land_albedo_4 )
      cube = cube_from_array(land_albedo(:,:,4))

    CASE ( var_id_latent_heat )
      cube = cube_from_array(latent_heat)

    CASE ( var_id_ls_rain )
      cube = cube_from_array(ls_rain)

    CASE ( var_id_ls_snow )
      cube = cube_from_array(ls_snow)

    CASE ( var_id_lw_down )
      cube = cube_from_array(lw_down)

    CASE ( var_id_precip )
      cube = cube_from_array(ls_rain + con_rain + ls_snow + con_snow)

    CASE ( var_id_pstar )
      cube = cube_from_array(pstar)

    CASE ( var_id_q1p5m_gb )
      cube = cube_from_array(q1p5m)

    CASE ( var_id_qw1 )
      cube = cube_from_array(qw_1)

    CASE ( var_id_rainfall )
      cube = cube_from_array(ls_rain + con_rain)

    CASE ( var_id_snomlt_surf_htf )
      cube = cube_from_array(snomlt_surf_htf)

    CASE ( var_id_snowfall )
      cube = cube_from_array(ls_snow + con_snow)

    CASE ( var_id_snow_mass_gb )
!     Don't use the snow_mass variable as that is calculated under control.
      cube_land = cube_from_array(tiles_to_gbm(snow_grnd + snow_tile))
      cube = map_from_land(cube_land)

    CASE ( var_id_surf_ht_flux_gb )
      cube = cube_from_array(surf_ht_flux)

    CASE ( var_id_sw_down )
      cube = cube_from_array(sw_down)

    CASE ( var_id_t1p5m_gb )
      cube = cube_from_array(t1p5m)

    CASE ( var_id_taux1 )
      cube = cube_from_array(taux_1)

    CASE ( var_id_tauy1 )
      cube = cube_from_array(tauy_1)

    CASE ( var_id_tl1 )
      cube = cube_from_array(tl_1)

    CASE ( var_id_tstar_gb )
!     Don't use tstar as that is calculated part-way through a timestep.
      cube_land = cube_from_array(tiles_to_gbm( tstar_tile ))
      cube = map_from_land(cube_land)

    CASE ( var_id_u1 )
      cube = cube_from_array(u_1)

    CASE ( var_id_u10m )
      cube = cube_from_array(u10m)

    CASE ( var_id_v1 )
      cube = cube_from_array(v_1)

    CASE ( var_id_v10m )
      cube = cube_from_array(v10m)

    CASE ( var_id_wind )
      cube = cube_from_array(SQRT(u_1**2 + v_1**2))

    CASE ( var_id_ext )
      cube_land = cube_from_array(ext)
      cube = map_from_land(cube_land)

    CASE ( var_id_smcl )
      cube_land = cube_from_array(smcl)
      cube = map_from_land(cube_land)

    CASE ( var_id_soil_wet )
      cube_land = cube_from_array(sthu + sthf)
      cube = map_from_land(cube_land)

    CASE ( var_id_sthf )
      cube_land = cube_from_array(sthf)
      cube = map_from_land(cube_land)

    CASE ( var_id_sthu )
      cube_land = cube_from_array(sthu)
      cube = map_from_land(cube_land)

    CASE ( var_id_t_soil )
      cube_land = cube_from_array(t_soil)
      cube = map_from_land(cube_land)

    CASE ( var_id_alb_tile_1 )
      cube_land = cube_from_array(alb_tile(:,:,1))
      cube = map_from_land(cube_land)

    CASE ( var_id_alb_tile_2 )
      cube_land = cube_from_array(alb_tile(:,:,2))
      cube = map_from_land(cube_land)

    CASE ( var_id_alb_tile_3 )
      cube_land = cube_from_array(alb_tile(:,:,3))
      cube = map_from_land(cube_land)

    CASE ( var_id_alb_tile_4 )
      cube_land = cube_from_array(alb_tile(:,:,4))
      cube = map_from_land(cube_land)

    CASE ( var_id_anthrop_heat )
      cube_land = cube_from_array(anthrop_heat)
      cube = map_from_land(cube_land)

    CASE ( var_id_canopy )
      cube_land = cube_from_array(canopy)
      cube = map_from_land(cube_land)

    CASE ( var_id_catch )
      cube_land = cube_from_array(catch)
      cube = map_from_land(cube_land)

    CASE ( var_id_ecan )
      cube_land = cube_from_array(ecan_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_ei )
      cube_land = cube_from_array(ei_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_emis )
      cube_land = cube_from_array(emis_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_esoil )
      cube_land = cube_from_array(esoil_tile)
      cube = map_from_land(cube_land)

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
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_ftl )
      cube_land = cube_from_array(ftl_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_gc )
      cube_land = cube_from_array(gc)
      cube = map_from_land(cube_land)

    CASE ( var_id_le )
      cube_land = cube_from_array(le_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_nsnow )
      cube_land = cube_from_array(REAL(nsnow))
      cube = map_from_land(cube_land)

    CASE ( var_id_q1p5m )
      cube_land = cube_from_array(q1p5m_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_rad_net_tile )
      cube_land = cube_from_array(radnet_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_rgrain )
      cube_land = cube_from_array(rgrain)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_can_melt )
! Only include tiles where canopy snow model is used
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) workspace_tile(:,i) = melt_tile(:,i)
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_can )
! Only include tiles where canopy snow model is used
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) workspace_tile(:,i) = snow_tile(:,i)
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_depth )
      cube_land = cube_from_array(snowDepth)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_grnd_rho )
      cube_land = cube_from_array(rho_snow_grnd)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_grnd )
! Only include tiles where canopy snow model is used
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) workspace_tile(:,i) = snow_grnd(:,i)
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_ground )
      DO i = 1,ntiles
        IF ( canSnowTile(i) ) THEN
          workspace_tile(:,i) = snow_grnd(:,i)
        ELSE
          workspace_tile(:,i) = snow_tile(:,i)
        ENDIF
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_ice_tile )
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sice(i,n,1:nsnow(i,n)))
        ENDDO
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_liq_tile )
      DO n = 1,ntiles
        DO j = 1,tile_pts(n)
          i = tile_index(j,n)
          workspace_tile(i,n) = SUM(sliq(i,n,1:nsnow(i,n)))
        ENDDO
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_mass )
      workspace_tile(:,:) = snow_tile(:,:)
! Add snow below canopy
      DO n = 1,ntiles
        IF ( canSnowTile(n) )                                                 &
          workspace_tile(:,n) = workspace_tile(:,n) + snow_grnd(:,n)
      ENDDO
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_melt )
      workspace_tile(:,:) = melt_tile(:,:)
! Add melting of snow below canopy.
!     IF ( can_model == 4 ) workspace_tile(:,:) = workspace_tile(:,:) + snowGMeltDiag(:,:)
      cube_land = cube_from_array(workspace_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_surf_ht_flux )
      cube_land = cube_from_array(surf_htf_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_surf_ht_store )
      cube_land = cube_from_array(surf_ht_store)
      cube = map_from_land(cube_land)

    CASE ( var_id_t1p5m )
      cube_land = cube_from_array(t1p5m_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_tstar )
      cube_land = cube_from_array(tstar_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_z0 )
      cube_land = cube_from_array(z0_tile)
      cube = map_from_land(cube_land)

    CASE ( var_id_tile_index )
      cube_land = cube_from_array(REAL(tile_index))
      cube = map_from_land(cube_land)

    CASE ( var_id_isoprene_gb )
      cube_land = cube_from_array(isoprene)
      cube = map_from_land(cube_land)

    CASE ( var_id_isoprene )
      cube_land = cube_from_array(isoprene_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_terpene_gb )
      cube_land = cube_from_array(terpene)
      cube = map_from_land(cube_land)

    CASE ( var_id_terpene )
      cube_land = cube_from_array(terpene_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_methanol_gb )
      cube_land = cube_from_array(methanol)
      cube = map_from_land(cube_land)

    CASE ( var_id_methanol )
      cube_land = cube_from_array(methanol_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_acetone_gb )
      cube_land = cube_from_array(acetone)
      cube = map_from_land(cube_land)

    CASE ( var_id_acetone )
      cube_land = cube_from_array(acetone_ft)
      cube = map_from_land(cube_land)

    CASE ( var_id_rgrainl )
      cube_land = cube_from_array(rgrainl)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_ds )
      cube_land = cube_from_array(ds)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_ice )
      cube_land = cube_from_array(sice)
      cube = map_from_land(cube_land)

    CASE ( var_id_snow_liq )
      cube_land = cube_from_array(sliq)
      cube = map_from_land(cube_land)

    CASE ( var_id_tsnow )
      cube_land = cube_from_array(tsnow)
      cube = map_from_land(cube_land)

    CASE DEFAULT
      CALL log_fatal("extract_var",                                           &
                     "Unrecognised variable for output - '" //                &
                     TRIM(get_string_identifier(var_id)) // "'. " //          &
                     "See docs for available variables")
  END SELECT

!-----------------------------------------------------------------------------
! Free the land cube
!
! Note that is safe to call this routine even if the cube has not been
! allocated - it just won't do anything
!-----------------------------------------------------------------------------
  CALL cube_free(cube_land)


  RETURN

END FUNCTION extract_var
#endif
