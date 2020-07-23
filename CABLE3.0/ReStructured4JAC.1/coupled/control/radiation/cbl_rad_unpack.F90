module cable_rad_unpack_mod
  
contains
 
subroutine cable_rad_unpack( &
land_albedo,            & 
alb_surft,              &
mp,                     &
nrb,                    &
row_length,             &
rows,                   &
land_pts,               &
nsurft,                 &
sm_levels,              &
tile_pts,               &
tile_index,             &
land_index,             &
tile_frac,              &
L_tile_pts,             &
EffSurfRefl_dif,        &
EffSurfRefl_beam,       &
snow_flag_cable,        &
SnowFlag_3L             & 


                           )

implicit none

!___ re-decl input args

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


!___ local vars
INTEGER :: i,J,K,L,N
REAL :: miss = 0.0
 
! std template args 
character(len=*), parameter :: subr_name = "cable_rad_unpack"
real :: Sumreffbm(mp) 
real :: Sumreffdf(mp) 
  
! only for land points, at present do not have a method for treating 
! mixed land/sea or land/seaice points as yet.
!   Albedo for surface tiles
!     (:,:,1) direct beam visible
!     (:,:,2) diffuse visible
!     (:,:,3) direct beam near-IR
!     (:,:,4) diffuse near-IR
alb_surft(:,:,:) = 0.0
alb_surft(:,:,1) = UNPACK(EffSurfRefl_beam(:,1),L_TILE_PTS, miss)
alb_surft(:,:,3) = UNPACK(EffSurfRefl_beam(:,2),L_TILE_PTS, miss)
alb_surft(:,:,2) = UNPACK(EffSurfRefl_dif(:,1),L_TILE_PTS, miss)
alb_surft(:,:,4) = UNPACK(EffSurfRefl_dif(:,2),L_TILE_PTS, miss)

do i=1,land_pts
  do j=1,nsurft
              
  if( alb_surft(i,j,1)> 1.) then
    print *, 'alb > 1',alb_surft(i,j,1)
    STOP
  ENDIF
 
  ENDDO
ENDDO

LAND_ALBEDO=0

DO N=1,nsurft
  DO K=1,TILE_PTS(N)
    L = TILE_INDEX(K,N)
    J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
    I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
    
    ! direct beam visible
    LAND_ALBEDO(I,J,1) = LAND_ALBEDO(I,J,1) +              &
                               TILE_FRAC(L,N)*ALB_surft(L,N,1)

    ! diffuse beam visible
    LAND_ALBEDO(I,J,2) = LAND_ALBEDO(I,J,2) +               &
                               TILE_FRAC(L,N)*ALB_surft(L,N,2)

    ! direct beam nearinfrared 
    LAND_ALBEDO(I,J,3) = LAND_ALBEDO(I,J,3) +               &
                               TILE_FRAC(L,N)*ALB_surft(L,N,3)

    ! diffuse beam nearinfrared
    LAND_ALBEDO(I,J,4) = LAND_ALBEDO(I,J,4) +               &
                               TILE_FRAC(L,N)*ALB_surft(L,N,4)
  ENDDO
ENDDO

return
     
End subroutine cable_rad_unpack
End module cable_rad_unpack_mod
