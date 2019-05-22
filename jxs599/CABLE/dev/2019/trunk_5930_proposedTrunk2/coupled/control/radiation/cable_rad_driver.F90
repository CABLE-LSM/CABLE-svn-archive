module cable_rad_driv_mod
  
contains

subroutine cable_rad_driver(  &
#                           include "cbl_rad_driver_args.inc"
                             )

!subrs:HaC1.3
USE cbl_hruff_mod,              ONLY : HgtAboveSnow
USE cbl_LAI_eff_mod,            ONLY : LAI_eff
USE cbl_albedo_mod,             ONLY : albedo
USE cbl_init_radiation_module,  ONLY : init_radiation

implicit none

!___ re-decl input args
#include "cbl_rad_driver_decs.inc"
 
character(len=*), parameter :: subr_name = "cable_rad_driver"
    
!set Height of Canopy above snow (rough%hruff) 
call HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft_cbl, & 
                   SnowDepth, SnowDensity )

!set Effective LAI considering ipotentail snow coverage
call LAI_eff( mp, LAI_pft_cbl, HGT_pft_cbl, HeightAboveSnow, &
                reducedLAIdue2snow)

!Defines Extinction Coefficients to use in calculation of Canopy 
!Reflectance/Transmitance. 
CALL init_radiation(  &
#                   include "cbl_init_radiation_args.inc"                  
                   )

!Finally call albedo to get what we really need to fill contract with JULES
!Defines 4-"band" albedos [VIS/NIR bands. direct beam/diffuse components] from 
!considering albedo of Ground (snow?) and Canopy Reflectance/Transmitance. 
call Albedo(  &
#             include "cbl_albedo_args.inc" 
           )

return
 
END SUBROUTINE cable_rad_driver
 
End module cable_rad_driv_mod

