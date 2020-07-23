!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
! 
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located 
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Converts selected CABLE hydrology variables to UM variables for 
!          UM hydrology code
!
! Called from: UM code hydrol
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Written for CABLE v1.8. No change for CABLE v2.0
!          In future could be combined with standard unpacking of cable 
!          variables at end of implicit call
!
!          2018 WB_LAKE extraction removed and placed in a separate routine
!          as needed for interation with ACCESS rivers' scheme
!
! ==============================================================================
module cable_hyd_driv_mod
  
contains

SUBROUTINE cable_hyd_driver( land_pts, ntiles, L_tile_pts, lying_snow, SNOW_TILE, SURF_ROFF,&
                             SUB_SURF_ROFF, TOT_TFALL, ssnow, canopy, veg )

  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
  USE cable_data_module, ONLY : cable
  !from old version
  !USE cable_common_module!, only : cable_runtime, cable_user
  USE cable_data_module,   ONLY : PHYS, OTHER
  USE cable_um_tech_mod, only : um1 

USE cable_canopy_type_mod,    ONLY : canopy_type
USE cable_soil_snow_type_mod, ONLY : soil_snow_type
USE cable_params_mod,         ONLY : veg_parameter_type

  USE cable_common_module, ONLY : cable_runtime, cable_user, l_casacnp,       &
                                  l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl
  
  implicit none
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

  !___ re-decl input args
  integer :: land_pts, ntiles
  LOGICAL :: L_tile_pts(land_pts,ntiles)     ! packing mask 
  
  REAL, INTENT(OUT), DIMENSION(LAND_PTS,NTILES) ::                    &
    SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        

  REAL, INTENT(OUT), DIMENSION(LAND_PTS) ::                               &
    LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
    SUB_SURF_ROFF, & !
    SURF_ROFF,     & !
    TOT_TFALL        !

  !___ local vars

  REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
    SURF_CAB_ROFF,    &
    TOT_TFALL_TILE                

  REAL :: miss =0. 
  REAL, POINTER :: TFRZ
 
  ! std template args 
  character(len=*), parameter :: subr_name = "cable_explicit_main"

  !-------- Unique subroutine body -----------
 
  TFRZ => PHYS%TFRZ
   
  SNOW_TILE= UNPACK(ssnow%snowd, L_TILE_PTS, miss) 

  LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass

  SURF_CAB_ROFF  = UNPACK(ssnow%rnof1, L_TILE_PTS, miss)
  SURF_ROFF      = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

  SURF_CAB_ROFF  = UNPACK(ssnow%rnof2, L_TILE_PTS, miss)
  SUB_SURF_ROFF  = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

  ! %through is /dels in UM app. for STASH output  
  canopy%through = canopy%through / kwidth_gl
  TOT_TFALL_TILE = UNPACK(canopy%through, L_TILE_PTS, miss)
  TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)

  !-------- End Unique subroutine body -----------
  
return

END SUBROUTINE cable_hyd_driver

!H!SUBROUTINE cable_lakesrivers(TOT_WB_LAKE)
!H!    
!H!    USE cable_um_tech_mod, ONLY : um1, ssnow
!H!    !jhan : thisversion of L_tile_pts is currently inctive
!H!    USE cbl_masks_mod, ONLY : L_tile_pts
!H!    
!H!    IMPLICIT NONE
!H!    
!H!    !routine extracts daily integrated ssnow%totwblake - water added to keep 
!H!    !lake tiles saturated - and grid cell averages (over land fraction)
!H!    !for use in river flow scaling routines
!H!    
!H!    REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) :: TOT_WB_LAKE
!H!    
!H!    !working variables
!H!    REAL :: miss = 0.
!H!    REAL, DIMENSION(um1%LAND_PTS, um1%ntiles) :: TOT_WB_LAKE_TILE
!H!    
!H!    TOT_WB_LAKE_TILE = UNPACK(ssnow%totwblake, um1%L_TILE_PTS, miss)
!H!    TOT_WB_LAKE = SUM(um1%TILE_FRAC * TOT_WB_LAKE_TILE,2)
!H!      
!H!    !zero the current integration
!H!    ssnow%totwblake = 0.
!H!
!H!END SUBROUTINE cable_lakesrivers

End module cable_hyd_driv_mod



