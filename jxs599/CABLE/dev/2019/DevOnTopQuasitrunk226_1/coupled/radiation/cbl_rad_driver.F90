module cable_rad_driv_mod
  
contains

subroutine cable_rad_driver(  &
#                           include "cbl_rad_driver_args.inc"
                             )

!subrs:HaC1.3
use cbl_masks_mod,              ONLY : fveg_mask, fsunlit_mask
use cbl_masks_mod,              ONLY : fsunlit_veg_mask
USE cbl_hruff_mod,              ONLY : HgtAboveSnow
USE cbl_LAI_eff_mod,            ONLY : LAI_eff
USE cbl_albedo_mod,             ONLY : albedo
USE cbl_init_radiation_module,  ONLY : init_radiation

implicit none

!___ re-decl input args
#include "cbl_rad_driver_decs.inc"
 
character(len=*), parameter :: subr_name = "cable_rad_driver"

logical :: cbl_standalone = .FALSE.    

!set Height of Canopy above snow (rough%hruff) 
call HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft_cbl, & 
                   SnowDepth, SnowDensity )

!set Effective LAI considering ipotentail snow coverage
call LAI_eff( mp, LAI_pft_cbl, HGT_pft_cbl, HeightAboveSnow, &
                reducedLAIdue2snow)

!define logical masks that are used throughout
call fveg_mask( veg_mask, mp, Clai_thresh, reducedLAIdue2snow )
call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, coszen )
call fsunlit_veg_mask( sunlit_veg_mask, mp, veg_mask, sunlit_mask )

!Defines Extinction Coefficients to use in calculation of Canopy 
!Reflectance/Transmitance. 
CALL init_radiation( &!rad%extkb, rad%extkd,                                     &
                     ExtCoeff_beam, ExtCoeff_dif,                              &
                     !rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
                     EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,              &
                     c1, rhoch, xk,                                            &
                     mp,nrb,                                                   &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,         &
                     cbl_standalone, jls_standalone, jls_radiation,            &
                     subr_name,                                                &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                   &
                     !veg%Xfang, veg%taul, veg%refl,                            &
                     VegXfang, VegTaul, VegRefl,                               & 
                     !met%coszen, int(met%DoY), met%fsd,                        &
                     coszen, metDoY, SW_down,                                  &
                     !canopy%vlaiw                                              &
                     reducedLAIdue2snow )
 
!Finally call albedo to get what we really need to fill contract with JULES
!Defines 4-"band" albedos [VIS/NIR bands. direct beam/diffuse components] from 
!considering albedo of Ground (snow?) and Canopy Reflectance/Transmitance. 

  call Albedo( &!ssnow%AlbSoilsn, soil%AlbSoil,                                &
               AlbSnow, AlbSoil,                                             &
               mp, nrb,                                                      &
               jls_radiation,                                                &
               veg_mask, sunlit_mask, sunlit_veg_mask,                       &  
               Ccoszen_tols, CGAUSS_W,                                       & 
               !veg%iveg, veg%refl, veg%taul,                                 & 
               surface_type, VegRefl, VegTaul,                               &
               !met%tk, met%coszen, canopy%vlaiw,                             &
               metTk, coszen, reducedLAIdue2snow,                            &
               !ssnow%snowd, ssnow%osnowd, ssnow%isflag,                      & 
               SnowDepth, SnowODepth, SnowFlag_3L,                           &
               !ssnow%ssdnn, ssnow%tgg(:,1), ssnow%snage,                     & 
               SnowDensity, SoilTemp, SnowAge,                               &
               xk, c1, rhoch,                                                & 
               !rad%fbeam, rad%albedo,                                        &
               RadFbeam, RadAlbedo,                                          & 
               !rad%extkd, rad%extkb,                                         & 
               ExtCoeff_dif, ExtCoeff_beam,                                  &   
               !rad%extkdm, rad%extkbm,                                       & 
               EffExtCoeff_dif, EffExtCoeff_beam,                             & 
               !rad%rhocdf, rad%rhocbm,                                       &
               CanopyRefl_dif,CanopyRefl_beam,                               &
               !rad%cexpkdm, rad%cexpkbm,                                     & 
               CanopyTransmit_dif, CanopyTransmit_beam,                      & 
               !rad%reffdf, rad%reffbm                                        &
               EffSurfRefl_dif, EffSurfRefl_beam )


return
 
END SUBROUTINE cable_rad_driver
 
End module cable_rad_driv_mod

