!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: cable_explicit_driver
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
MODULE cable_expl_unpack_mod

implicit none

contains

SUBROUTINE cable_expl_unpack( latitude, longitude, FTL_TILE, FQW_TILE,       &
                           TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow_snowd, ssnow_cls, air_rlam, air_rho,          &
                           canopy_fe, canopy_fh, canopy_us, canopy_cdtq,       &
                           canopy_fwet, canopy_wetfac_cs, canopy_rnet,         &
                           canopy_zetar, canopy_epot, met_ua, rad_trad,        &
                           rad_transd, rough_z0m, rough_zref_tq )

   USE cable_def_types_mod, ONLY : mp, NITER 
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1, basic_diag
   USE cable_common_module!, ONLY : cable_runtime, cable_user, &
                          !         ktau_gl, knode_gl, kend_gl, &
                          !         fudge_out, L_fudge,         &
                          !         fprintf_dir_root, fprintf_dir
   
  USE cable_unpack_checks_mod, ONLY : cable_unpack_check
!fprintf
  USE cable_common_module, ONLY :    fprintf_dir_root, fprintf_dir
  USE cable_diag_module

use arraydiag_m

  IMPLICIT NONE         


  character(len=*), parameter :: subr_name = "cable_explicit_unpack"
  !-------------------------------------------------------------------------- 
  !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
  !-------------------------------------------------------------------------- 


  !___ UM variables to recieve unpacked CABLE vars

  REAL,  DIMENSION(um1%row_length,um1%rows) ::                             &
     latitude,   &
     longitude

  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
     TSTAR_TILE
       
  !___return fluxes
  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
     FTL_TILE,   &  ! Surface FTL for land tiles     
     FQW_TILE      ! Surface FQW for land tiles     

  !___return temp and roughness
  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
     Z0H_TILE,         &
     Z0M_TILE

  !___return friction velocities/drags/ etc
  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
     CD_TILE,    &     ! Drag coefficient
     CH_TILE,    &     ! Transfer coefficient for heat & moisture
     U_S_STD_TILE      ! Surface friction velocity
  REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  :: &
     U_S               ! Surface friction velocity (m/s)

  !___return miscelaneous 
  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
     RADNET_TILE,   & ! Surface net radiation
     RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                      ! factor for fraction (1-FRACA) of snow-free land tiles
     RESFT,         & ! Total resistance factor.
                      ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                      ! 1 for snow.    
     FRACA,         & ! Fraction of surface moisture
     RECIP_L_MO_TILE,&! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
     EPOT_TILE

     
  LOGICAL,DIMENSION(um1%land_pts,um1%ntiles) :: l_tile_pts

  !___UM vars used but NOT returned 
  REAL, INTENT(IN), DIMENSION(um1%land_pts) ::   &
     FLAND(um1%land_pts)              ! IN Land fraction on land tiles.


  !___ decs of intent(in) CABLE variables to be unpacked

  ! snow depth (liquid water), factor for latent heat
  REAL, INTENT(IN), DIMENSION(mp) :: ssnow_snowd, ssnow_cls
  
  ! surface wind speed (m/s)
  REAL, INTENT(IN), DIMENSION(mp) :: met_ua 
  
  ! latent heat for water (j/kg), dry air density (kg m-3)
  REAL, INTENT(IN), DIMENSION(mp) :: air_rlam, air_rho 
  
  ! frac SW diffuse transmitted thru canopy, rad. temp. (soil and veg)
  REAL, INTENT(IN), DIMENSION(mp) :: rad_trad,rad_transd 
  
  ! total latent heat (W/m2), total sensible heat (W/m2)
  REAL, INTENT(IN), DIMENSION(mp) :: canopy_fe, canopy_fh  
  
  ! fraction of canopy wet
  REAL, INTENT(IN), DIMENSION(mp) :: canopy_fwet, canopy_wetfac_cs
  
  ! friction velocity, drag coefficient for momentum
  REAL, INTENT(IN), DIMENSION(mp) :: canopy_us, canopy_cdtq
  
  ! net rad. absorbed by surface (W/m2), total potential evaporation 
  REAL, INTENT(IN), DIMENSION(mp) :: canopy_rnet, canopy_epot        
  
  ! stability correction
  REAL, INTENT(IN), DIMENSION(mp,niter) :: canopy_zetar
  
  ! roughness length, Reference height for met forcing
  REAL, INTENT(IN), DIMENSION(mp) :: rough_z0m, rough_zref_tq 

  !-------------------------------------------------------------------------- 
  !--- end INPUT ARGS FROM cable_explicit_driver() ------------------------------
  !-------------------------------------------------------------------------- 
  
  !___vars in local calc. of latent heat fluxes
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
     LE_TILE

  !___vars in local calc of Surface friction velocities
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
     U_S_TILE
  REAL, DIMENSION(mp)  :: &
     CDCAB
  
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
     Lpts_nTILE


  !___local miscelaneous
  REAL, DIMENSION(mp)  :: &
    THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
  INTEGER :: i,j,k,N,L
  REAL :: miss = 0.0
  LOGICAL, SAVE :: first_cable_call = .true.
  REAL, POINTER :: CAPP 
  LOGICAL :: Lunpack = .false.
  LOGICAL :: checks = .false.
  real :: tols = 0.25

  !fprintf{_____________________________________________________________________  
  !USE cable_common_module, ONLY :    fprintf_dir_root, fprintf_dir
  !USE cable_diag_module
  character(len=25) :: vname
  character(len=70) :: dir
 
  INTEGER, SAVE ::                                                            &
   cDiag00=0, cDiag0=0, cDiag1=0, cDiag2=0, cDiag3=0, cDiag4=0,                &
   cDiag5=0, cDiag6=0, cDiag7=0, cDiag8=0, cDiag9=0, cDiag10=0,    cDiag11=0,  &  
   cDiag12=0, cDiag13=0, cDiag14=0, cDiag15=0, cDiag16=0, cDiag17=0, cDiag18=0,& 
   cDiag19=0,cDiag20=0, cDiag21=0, cDiag22=0, cDiag23=0, cDiag24=0, cDiag25=0, &
   cDiag26=0, cDiag27=0, cDiag28=0, cDiag29=0, cDiag30=0, cDiag31=0, cDiag32=0,&  
   cDiag33=0, cDiag34=0, cDiag35=0, cDiag36=0, cDiag37=0, cDiag38=0, cDiag39=0,& 
   cDiag40=0, cDiag41=0, cDiag42=0, cDiag43=0, cDiag44=0, cDiag45=0, cDiag46=0,& 
   cDiag47=0, cDiag48=0, cDiag49=0, cDiag50=0, cDiag51=0, cDiag52=0, cDiag53=0,& 
   cDiag54=0, cDiag55=0, cDiag56=0, cDiag57=0, cDiag58=0, cDiag59=0, cDiag60=0,& 
   cDiag61=0, cDiag62=0, cDiag63=0, cDiag64=0, cDiag65=0, cDiag66=0, cDiag67=0,& 
   cDiag68=0, cDiag69=0

  logical, parameter :: L_fprint_HW = .false.
  logical :: L_fprint

  L_fprint = .false. !default

  !if( L_fprint_HW ) then
  !  if ( ktau_gl==1 .OR. ktau_gl==54 .OR. ktau_gl==154 .OR. &
  !       ktau_gl==154 .OR. ktau_gl==154 .OR. mod(ktau_gl,10)==0. ) then
  !    L_fprint = .true.
  !  endif  
  !endif  

  fprintf_dir=trim(fprintf_dir_root)//trim("expl_unpack")//"/"
  
  !vname='' 
  !call cable_fprintf( cDiagX, vname, %, mp, L_fprint )
  !fprintf_____________________________________________________________________} 
        
  CAPP => PHYS%CAPP
     
  !___return fluxes

  if( L_fprint_HW ) L_fprint = .true.
  
  vname='canopy_fh' ! FTL_TILE 
  call cable_fprintf( cDiag0, vname, canopy_fh/CAPP, mp, L_fprint )
  FTL_TILE = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
  FTL_TILE = FTL_TILE / CAPP

  vname='canopy_fe' ! FQW_TILE 
  fe_dlh = canopy_fe/(air_rlam*ssnow_cls)
  call cable_fprintf( cDiag1, vname, fe_dlh, mp, L_fprint )
  FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

  !___return temp and roughness
  vname='rad_trad' ! FTL_TILE 
  call cable_fprintf( cDiag2, vname, rad_trad, mp, L_fprint )
    !testing TSTAR_TILE v miss: test if makes difference outside CABLE
    !Lpts_nTILE = tstar_tile
    !TSTAR_TILE = UNPACK(rad_trad,  um1%l_tile_pts, Lpts_nTILE )
    !call arraydiag("jhA_Tstar",tstar_tile)
  TSTAR_TILE = UNPACK(rad_trad,  um1%l_tile_pts, miss )
    !call arraydiag("jhB_Tstar",tstar_tile)




Z0M_TILE = UNPACK(rough_z0m,  um1%l_tile_pts, miss)
Z0H_TILE = Z0M_TILE
      
      !___return friction velocities/drags/ etc
U_S_TILE  =  UNPACK(canopy_us, um1%l_tile_pts, miss)
      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
      ! for Cable CD*
CD_TILE =  UNPACK(CDCAB,um1%l_tile_pts, miss)
      ! for Cable CH*
CH_TILE =  UNPACK(canopy_cdtq,um1%l_tile_pts, miss)

U_S_STD_TILE=U_S_TILE

U_S = 0.
      DO N=1,um1%ntiles
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
            I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
U_S(I,J) = U_S(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*U_S_TILE(L,N)
         ENDDO
      ENDDO


      !___return miscelaneous 
      fraca_cab = canopy_fwet * (1.-rad_transd)
      WHERE( ssnow_snowd > 1.0 ) fraca_cab = 1.0
      rfsfs_cab = MIN( 1., MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /         &
                  MAX( 0.01,1. - fraca_cab ) )
FRACA = UNPACK( fraca_cab, um1%l_tile_pts, miss )
RESFT = UNPACK( canopy_wetfac_cs,um1%l_tile_pts, miss )
RESFS = UNPACK( rfsfs_cab , um1%l_tile_pts, miss )

RADNET_TILE = UNPACK( canopy_rnet , um1%l_tile_pts, miss )

      THETAST = ABS( canopy_fh ) / ( air_rho * capp*canopy_us )
      RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
RECIP_L_MO_TILE = UNPACK( RECIPLMOTILE, um1%l_tile_pts, miss )
EPOT_TILE = UNPACK( canopy_epot, um1%l_tile_pts, miss )

      IF(first_cable_call) THEN 
         l_tile_pts = um1%l_tile_pts
         first_cable_call = .FALSE.
      ENDIF


END SUBROUTINE cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

END MODULE cable_expl_unpack_mod
