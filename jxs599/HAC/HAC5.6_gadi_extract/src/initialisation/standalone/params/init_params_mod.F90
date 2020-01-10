#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE init_params_mod

USE logging_mod, ONLY: log_info, log_debug, log_warn, log_error, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC init_params
PUBLIC init_veg_from_vegin_JAC

CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "init_params.inc"
#include "init_pftparm.inc"
#include "init_nvegparm.inc"
#include "init_cropparm.inc"
#include "init_triffid.inc"
#include "init_soilparm_cable.inc"
#include "init_deposition_species.inc"

SUBROUTINE init_veg_from_vegin_JAC(VegXfang, VegTaul, VegRefl, fmp, fnrb, SurfaceType, VeginXfang, VeginTaul, VeginRefl ) 
USE ancil_info,                 ONLY : nsurft
USE ancil_info,                 ONLY : lp => land_pts 
USE cable_types_mod,            ONLY : L_tile_pts
implicit none
integer ::  fmp
integer ::  fnrb
integer ::  SurfaceType(nsurft)
real :: VegXfang(fmp)
real :: VegTaul(fmp, fnrb)
real :: VegRefl(fmp, fnrb)
real :: VeginXfang(nsurft)
real :: VeginTaul(fnrb, nsurft)
real :: VeginRefl( fnrb, nsurft)
real :: VeginXfang_lp(lp, nsurft)
real :: VeginTaul1_lp(lp, nsurft)
real :: VeginTaul2_lp(lp, nsurft)
real :: VeginRefl1_lp(lp, nsurft)
real :: VeginRefl2_lp(lp, nsurft)
integer :: i, h

DO i = 1, nsurft
  VeginXfang_lp(:,i) = VeginXfang(i)
  VeginTaul1_lp(:,i) = VeginTaul(1,i)
  VeginTaul2_lp(:,i) = VeginTaul(2,i)
  VeginRefl1_lp(:,i) = VeginRefl(1,i)
  VeginRefl2_lp(:,i) = VeginRefl(2,i)
EndDo

VegXfang(:)  = PACK( VeginXfang_lp, L_tile_pts ) 
VegTaul(:,1) = PACK( VeginTaul1_lp, L_tile_pts ) 
VegTaul(:,2) = PACK( VeginTaul2_lp, L_tile_pts ) 
VegRefl(:,1) = PACK( VeginRefl1_lp, L_tile_pts )   
VegRefl(:,2) = PACK( VeginRefl2_lp, L_tile_pts )     
  
END SUBROUTINE init_veg_from_vegin_JAC 

END MODULE init_params_mod

#endif
