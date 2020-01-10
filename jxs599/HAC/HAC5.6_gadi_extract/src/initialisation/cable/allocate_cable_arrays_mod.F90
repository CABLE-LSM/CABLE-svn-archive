MODULE allocate_cable_arrays_mod

!Common Non-science modules

USE yomhook,                  ONLY: lhook, dr_hook
USE ereport_mod,              ONLY: ereport
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE
PUBLIC :: allocate_cable_arrays

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_CABLE_ARRAYS_MOD'

CONTAINS

SUBROUTINE allocate_cable_arrays()

USE cable_types_mod,      ONLY : mp
USE cable_types_mod,      ONLY : nrb
USE cable_types_mod,      ONLY : L_tile_pts
USE cable_types_mod,      ONLY : vegin    
USE cable_types_mod,      ONLY : SurfaceTypeID_cbl
USE cable_types_mod,      ONLY : VegXfang
USE cable_types_mod,      ONLY : VegTaul
USE cable_types_mod,      ONLY : VegRefl

USE ancil_info,               ONLY: land_pts, surft_pts, frac_surft
USE jules_surface_types_mod,  ONLY: ntype

!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------
INTEGER ::  error = 0

! Determine the number of active tiles
mp = SUM(surft_pts)

ALLOCATE( l_tile_pts(land_pts, ntype),  stat = error )
ALLOCATE( SurfaceTypeID_cbl(mp),     stat = error ) 
ALLOCATE( VegXfang(mp),                 stat = error )                
ALLOCATE( VegTaul(mp,nrb),              stat = error )             
ALLOCATE( VegRefl(mp,nrb),              stat = error )             

END SUBROUTINE allocate_cable_arrays

END MODULE allocate_cable_arrays_mod
