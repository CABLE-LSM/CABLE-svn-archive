
!model dimensions
!-------------------------------------------------------------------------------
!JaC:todo:ultimatelty get this from JaC~
integer :: mp                       !total number of "tiles"  
integer :: nrb                      !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]
integer :: row_length                       !grid cell x
integer :: rows                             !grid cell y
integer :: land_pts                         !grid cell land points on the x,y grid
integer :: nsurft                           !grid cell number of surface types 
integer :: sm_levels                        !grid cell number of soil levels 
integer :: tile_pts(nsurft)                 !Number of land points per PFT 
integer :: tile_index(land_pts,nsurft)      !Index of land point in (land_pts) array
integer :: land_index(land_pts)             !Index of land points in (x,y) array - see below
!-------------------------------------------------------------------------------

!Return from CABLE to vulvill contract with JULES
!-------------------------------------------------------------------------------
real :: land_albedo(row_length,rows,4)      
real :: alb_surft(Land_pts,nsurft,4)        
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. 
!-------------------------------------------------------------------------------
real :: tile_frac(land_pts,nsurft)          !fraction of each surface type per land point 
logical :: L_tile_pts( land_pts, nsurft )   !mask:=TRUE where tile_frac>0, else FALSE. pack mp according to this mask
!-------------------------------------------------------------------------------

!recieved as spatial maps from the UM. remapped to "mp"
!-------------------------------------------------------------------------------
integer:: surface_type(mp)          ! Integer index of Surface type (veg%iveg)
real :: LAI_pft_cbl(mp)             !LAI -  "limited" and remapped
real :: HGT_pft_cbl(mp)             !canopy height -  "limited" and remapped
!-------------------------------------------------------------------------------

!Prognostics
!-------------------------------------------------------------------------------
integer:: SnowFlag_3L(mp)           !Flag to treat snow as 3 layer  - if enough present. Updated depending on total depth (ssnow%isflag)
real :: snow_flag_cable(land_pts,nsurft) !Flag to treat snow as 3 layer  - REALized and in JULES dimensioonal format
!-------------------------------------------------------------------------------
                                                                                           
! Albedos
!-------------------------------------------------------------------------------
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!-------------------------------------------------------------------------------

