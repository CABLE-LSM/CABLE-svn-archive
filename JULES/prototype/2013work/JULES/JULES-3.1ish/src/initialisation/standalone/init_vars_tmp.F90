#if !defined(UM_JULES)
!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_vars_tmp.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_vars_tmp()

  USE Ancil_info
  USE Switches
  USE Top_pdm, ONLY : inlandout_atm
  USE Forcing, ONLY : u_0,v_0
  USE Prognostics
  USE Aero
  USE Orog
  USE Trifctl
  USE Coastal
  USE C_Densty
  USE Sea_ice
  USE C_kappai
  USE C_Rough
  USE p_s_parms, ONLY : satcon,soil_clay
  USE c_elevate, ONLY : z_land

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises various variables that may change their initialisation in
!   future versions
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: i,l  ! Loop counters


!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Initialise accumulated fluxes for TRIFFID and phenology.
! This is not necessary if these are read from a restart file - but at
! present they're not.
!-----------------------------------------------------------------------
  g_leaf_acc(:,:)      = 0.0
  g_leaf_phen_acc(:,:) = 0.0
  npp_ft_acc(:,:)      = 0.0
  resp_s_acc(:,:)      = 0.0
  resp_w_ft_acc(:,:)   = 0.0

!-----------------------------------------------------------------------
! Set saturated hydraulic conductivity to zero at land ice points.
!-----------------------------------------------------------------------
  IF ( lice_pts > 0 ) THEN
    CALL log_info("init_vars_tmp",                                            &
                  "Setting satcon to zero at land ice points")
    DO i=1,lice_pts
      l = lice_index(i)
      satcon(l,:) = 0.0
    ENDDO
  ENDIF

!-----------------------------------------------------------------------
! Set surface velocity to be zero
!-----------------------------------------------------------------------
  u_0(:,:) = 0.0
  v_0(:,:) = 0.0

!-----------------------------------------------------------------------
! Set CO2 variables
!-----------------------------------------------------------------------
  co2_3d(:,:) = 0.0

!-----------------------------------------------------------------------
! Set coastal tiling variables
!-----------------------------------------------------------------------
  tstar_sea(:,:)  = 280.0
  tstar_sice_ncat(:,:,:) = 270.0

!-----------------------------------------------------------------------
! Set orographic roughness variables
!-----------------------------------------------------------------------
  h_blend_orog(:,:) = 0.0
  sil_orog_land(:)  = 0.0
  ho2r2_orog(:)     = 0.0

!-----------------------------------------------------------------------
! Set up prognostics which are not currently in dump
!-----------------------------------------------------------------------
  ti(:,:)                   = 270.0
  z0msea(:,:)               = z0hsea
  snow_mass(:,:)            = 0.0
  ice_fract_ncat(:,:,:)     = 0.0
  di_ncat(:,:,:)            = 0.0
  snow_mass_sea_ncat(:,:,:) = 0.0
  k_sice(:,:,:)             = 2.0 * kappai / de

!-----------------------------------------------------------------------
! Set up sea-ice parameter variables
!-----------------------------------------------------------------------
  IF(L_SSICE_ALBEDO) THEN
    SW_alpham=0.65
    SW_alphac=0.80
    SW_alphab=0.57
    SW_dtice=2.00
  ELSE
    SW_alpham=0.50
    SW_alphac=0.80
    SW_alphab=-1.00
    SW_dtice=10.00
  ENDIF

  soot(:) = 0.0

!------------------------------------------------------------------
! Temporary initialisation of variables added during UM integration
!------------------------------------------------------------------
  z_land(:,:)     = 0.0
  soil_clay(:,:)  = 0.0
  
  INLANDOUT_ATM(:) = 0.0

  RETURN

END SUBROUTINE init_vars_tmp
#endif
