!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/model_interface_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-09-27 09:33:34 +0100 (Thu, 27 Sep 2012) $
!!   $LastChangedRevision: 609 $
!!
!!****************************************************************************

MODULE cable_vars_mod

! vars to pass to CABLE which are not ordinarily in JULES subrs are 
! inherited from here. clean this up to include df of above USESs

!#define CABLE_from_control T
#ifndef CABLE_from_control
  USE model_time_mod, ONLY : a_step => timestep
  USE snow_param, ONLY : rho_snow => rho_snow_const
#endif
  
#ifdef CABLE_from_control
  ! exists in sf_exch
  USE soil_param, ONLY : dzsoil
#endif
  
  USE ancil_info, ONLY : row_length, rows, sm_levels,DIM_CS1,DIM_CS2

  USE model_time_mod, ONLY : timestep_len, timestep_number => timestep

  USE forcing, ONLY : sw_down, lw_down, ls_rain, ls_snow 

  USE latlon_mod, ONLY : latitude, longitude !, & 

  USE p_s_parms, ONLY : hcon, satcon, sathh, smvcwt,     &
                        smvccl, albsoil, sthf,sthu,sthf, bexp => b 
  USE prognostics, ONLY : canht_ft, lai, t_soil, smcl, gs, canopy_gb

  USE aero, ONLY : co2_mmr

  !jhan: can use others here and eliminate passing
  ! some exists in sf_exch
  !USE fluxes, ONLY : ftl_1, fqw_1, melt_tile,&
  USE fluxes, ONLY : melt_tile,&
               ECAN_TILE,ESOIL_TILE,EI_TILE 
  USE screen, ONLY : T1P5M_TILE, Q1P5M_TILE
  USE trifctl 

  USE coastal, ONLY : surf_ht_flux_land

  USE nstypes, ONLY : npft 

  IMPLICIT NONE
   
   ! CABLE uses these in diag routines, dummies here
   INTEGER :: endstep=0, mype=0

   INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE ::                         & 
      isnow_flg3l   ! 3 layer snow scheme flag (make logical)

   ! CABLE-JULES vars to be initialized from file thru JULES i/o
   
   ! suffix number refers to layer. JULES i/o can't handle tiled 
   ! soil/snow multi layered vars
   
   !snow vars
   REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::                            &
      snow_rho1l,  & ! snow density for 1 layer scheme
      snage_tile,  & ! snow age 
      snow_rho1,   & ! 
      snow_rho2,   & !
      snow_rho3,   & !
      snow_mass1,  & !
      snow_mass2,  & !
      snow_mass3,  & !
      snow_depth1, & !
      snow_depth2, & !
      snow_depth3, & !
      snow_tmp1,   & !
      snow_tmp2,   & !
      snow_tmp3      !
   
   !soil vars
   REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::                            &
      sthu_tile1, &
      sthu_tile2, &
      sthu_tile3, &
      sthu_tile4, &
      sthu_tile5, &
      sthu_tile6, &
      sthf_tile1, &
      sthf_tile2, &
      sthf_tile3, &
      sthf_tile4, &
      sthf_tile5, &
      sthf_tile6, &
      smcl_tile1, &
      smcl_tile2, &
      smcl_tile3, &
      smcl_tile4, &
      smcl_tile5, &
      smcl_tile6, &
      tsoil_tile1, &
      tsoil_tile2, &
      tsoil_tile3, &
      tsoil_tile4, &
      tsoil_tile5, &
      tsoil_tile6
   
   ! End - CABLE-JULES vars to be initialized from file thru JULES i/o

   ! Re-packed CABLE-JULES vars from above for passing to CABLE
   ! consistency with UM coupling here
   
   REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::                          &
      snow_rho3l,   & !
      snow_mass3l,  & !
      snow_depth3l, & !
      snow_tmp3l      !
   
   REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::                          &
      sthu_tile, &
      sthf_tile, &
      smcl_tile, &
      tsoil_tile

   ! End - Re-packed CABLE-JULES vars from above for passing to CABLE

   ! CABLE specific var BUT doesnt need initializing 
   !jhan:  - value taken from .... ??
   REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::                            &
      snow_cond      !
   
   REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::                            &
      snow_cond1,  & ! snow conductance for layer 1 
      snow_cond2,  & ! snow conductance for layer 2 
      snow_cond3     ! snow conductance for layer 3 
   
   ! Re-pack in CABLE-JULES format
   REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: cosz_cable

   !jhan: initialize??(row_length,rows,4)
!   REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
!      surf_down_sw

!   REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
!     LAND_ALBEDO_CABLE
   
   !(land_points,NTILES)                                
   REAL, DIMENSION(:,:), ALLOCATABLE :: &
     RTSOIL_TILE, & !
     GFLUX_TILE,  & !
     SGFLUX_TILE, & !
     TOT_ALB,     & !
     conv_rain,   & ! convective rain
     conv_snow      ! convective snow

CONTAINS
   
SUBROUTINE setup_cable(timestep_number)
   
  USE ancil_info, ONLY : ntiles, land_pts, sm_levels, msn => nsmax 
   
  USE p_s_parms, ONLY : cosz

   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: timestep_number 
   INTEGER :: n,k

   IF( timestep_number == 1 ) THEN

      ALLOCATE(snow_rho3l(land_pts, ntiles, msn) )
      ALLOCATE(snow_mass3l(land_pts, ntiles, msn) )
      ALLOCATE(snow_depth3l(land_pts, ntiles, msn) )
      ALLOCATE(snow_tmp3l(land_pts, ntiles, msn) )
   
      ALLOCATE(sthu_tile(land_pts, ntiles, sm_levels) )
      ALLOCATE(sthf_tile(land_pts, ntiles, sm_levels) )
      ALLOCATE(smcl_tile(land_pts, ntiles, sm_levels) )
      ALLOCATE(tsoil_tile(land_pts, ntiles, sm_levels) )
   
      ALLOCATE(snow_cond(land_pts, ntiles, msn) )
      ALLOCATE(snow_cond1(land_pts, ntiles) )
      ALLOCATE(snow_cond2(land_pts, ntiles) )
      ALLOCATE(snow_cond3(land_pts, ntiles) )
   
      !ALLOCATE( surf_down_sw(row_length, rows, 4) )
      ALLOCATE( cosz_cable(row_length, rows) )
      ALLOCATE( conv_rain(row_length, rows) )
      ALLOCATE( conv_snow(row_length, rows) )
      
      !ALLOCATE( LAND_ALBEDO_CABLE(row_length,rows,4) )
      ALLOCATE( RTSOIL_TILE(land_pts,NTILES) )
      ALLOCATE( GFLUX_TILE(land_pts,NTILES) )
      ALLOCATE( SGFLUX_TILE(land_pts,NTILES) )
      ALLOCATE( TOT_ALB(land_pts,NTILES) )   


      ! Re-pcking here is done largely so inrerface between JULES/UM
      ! can remain as close as possible
      ! re-pack snow layerd vars into single array
      snow_cond1(:,:) = 0.
      snow_cond2(:,:) = 0.
      snow_cond3(:,:) = 0.
      snow_cond(:,:,1) = snow_cond1(:,:)
      snow_cond(:,:,2) = snow_cond2(:,:)
      snow_cond(:,:,3) = snow_cond3(:,:)
      
      snow_depth3l(:,:,1) = snow_depth1
      snow_depth3l(:,:,2) = snow_depth2
      snow_depth3l(:,:,3) = snow_depth3
      
      snow_mass3l(:,:,1) = snow_mass1
      snow_mass3l(:,:,2) = snow_mass2
      snow_mass3l(:,:,3) = snow_mass3
      
      snow_rho3l(:,:,1) = snow_rho1
      snow_rho3l(:,:,2) = snow_rho2
      snow_rho3l(:,:,3) = snow_rho3
      
      snow_tmp3l(:,:,1) = snow_tmp1
      snow_tmp3l(:,:,2) = snow_tmp2
      snow_tmp3l(:,:,3) = snow_tmp3
      
      sthu_tile(:,:,1) = sthu_tile1
      sthu_tile(:,:,2) = sthu_tile2
      sthu_tile(:,:,3) = sthu_tile3
      sthu_tile(:,:,4) = sthu_tile4
      sthu_tile(:,:,5) = sthu_tile5
      sthu_tile(:,:,6) = sthu_tile6
      
      sthf_tile(:,:,1) = sthf_tile1
      sthf_tile(:,:,2) = sthf_tile2
      sthf_tile(:,:,3) = sthf_tile3
      sthf_tile(:,:,4) = sthf_tile4
      sthf_tile(:,:,5) = sthf_tile5
      sthf_tile(:,:,6) = sthf_tile6
      
      smcl_tile(:,:,1) = smcl_tile1
      smcl_tile(:,:,2) = smcl_tile2
      smcl_tile(:,:,3) = smcl_tile3
      smcl_tile(:,:,4) = smcl_tile4
      smcl_tile(:,:,5) = smcl_tile5
      smcl_tile(:,:,6) = smcl_tile6

      tsoil_tile(:,:,1) = tsoil_tile1
      tsoil_tile(:,:,2) = tsoil_tile2
      tsoil_tile(:,:,3) = tsoil_tile3
      tsoil_tile(:,:,4) = tsoil_tile4
      tsoil_tile(:,:,5) = tsoil_tile5
      tsoil_tile(:,:,6) = tsoil_tile6

    ENDIF 
    
   ! Re-pack cosine zenith angle for CABLE 
   ! (UM processes as dimension row_length, rows)
    DO k=1,rows
      DO n=1,row_length
        cosz_cable(n,k) = cosz( n + ( (k-1)*row_length ) )
      END DO
    END DO

   !jhan:TODO initialize surf_sw_down?
   !jhan:TODO Eva can these be removed?
   !EAK  SW_DOWN calculations for checking only !!!
   !SW_DOWN = 0.
   !SW_DOWN = (surf_down_sw(:,:,1)+surf_down_sw(:,:,2) + &
   !   surf_down_sw(:,:,3)+surf_down_sw(:,:,4))*cosz_cable(:,:)

   conv_rain = 0. ! convective rain and snow calculated after 
                  ! call to ni_bl_ctl
   conv_snow = 0.

END SUBROUTINE setup_cable 


!jhan
!! existed in 7.3? in r2_swrad. 8.2 uses r2_swrad3z which also existed in 7.3,
!! BUT 8.2 does not have r2_swrad. thf. cannot init - 
!! surf_down_sw = surf_vis_dir_g.....
!
!   REAL, DIMENSION(row_length, rows, 4) ::                          &
!      surf_down_sw
!   

END MODULE cable_vars_mod


