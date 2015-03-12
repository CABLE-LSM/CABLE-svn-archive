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
! Purpose: Converts CABLE variables into JULES variables for the call to
!          surf_couple_extra
!
! Called from: cable_couple_extra
!
! Author: Matt Pryor (Met Office)
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Written for CABLE v1.8. No change for CABLE v2.0
!          In future could be combined with standard unpacking of cable 
!          variables at end of implicit call
!
!
! ==============================================================================

SUBROUTINE cable_extra_driver( SMVCST, TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,    &
                               STHF, STHF_TILE, STHU, STHU_TILE, SNOW_TILE,   &
                               SNOW_RHO1L, ISNOW_FLG3L, SNOW_DEPTH3L,         &
                               SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L,           &
                               SNOW_COND, SNAGE_TILE, CANOPY_TILE, CANOPY_GB, &
                               LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,          &
                               TOT_TFALL )

  USE cable_data_module,   ONLY : PHYS, OTHER
  USE cable_common_module!, only : cable_runtime, cable_user
  USE cable_um_tech_mod, only : um1, ssnow, canopy, veg, soil

  IMPLICIT NONE

  REAL, INTENT(IN), DIMENSION(um1%LAND_PTS) ::                                &
    SMVCST

  REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%SM_LEVELS) ::                 &
    SMCL,                                                                     &
    STHF,                                                                     &
    STHU,                                                                     &
    TSOIL

  REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%NTILES,um1%SM_LEVELS) ::      &
    SMCL_TILE,                                                                &
    STHU_TILE,                                                                &
    TSOIL_TILE,                                                               &
    STHF_TILE

   !___flag for 3 layer snow pack
  INTEGER, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%NTILES) ::                 &
    ISNOW_FLG3L

  REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%NTILES) ::                    &
    SNOW_TILE,                                                                &
    SNOW_RHO1L,                                                               &
    SNAGE_TILE,                                                               &
    CANOPY_TILE

   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
  REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles,3) ::                  &
    SNOW_DEPTH3L,                                                             &
    SNOW_MASS3L,                                                              &
    SNOW_RHO3L,                                                               &
    SNOW_TMP3L,                                                               &
    SNOW_COND

  REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) ::                               &
    CANOPY_GB,                                                                &
    LYING_SNOW,                                                               &
    SUB_SURF_ROFF,                                                            &
    SURF_ROFF,                                                                &
    TOT_TFALL

! Local variables
  REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
    SURF_CAB_ROFF,                                                            &
    TOT_TFALL_TILE

  REAL :: miss =0.
  INTEGER:: i_miss = 0
  REAL, POINTER :: TFRZ
  INTEGER :: i, j, k, n

!-----------------------------------------------------------------------------
      
  TFRZ => PHYS%TFRZ

  TSOIL_TILE = 0.0
  SMCL_TILE  = 0.0
  STHF_TILE  = 0.0
  STHU_TILE  = 0.0

  DO j = 1,um1%SM_LEVELS
    TSOIL_TILE(:,:,j) = UNPACK(ssnow%tgg(:,j), um1%L_TILE_PTS, miss)
    SMCL_TILE(:,:,j) = UNPACK(REAL(ssnow%wb(:,j)), um1%L_TILE_PTS, miss)
    SMCL_TILE(:,:,j) = SMCL_TILE(:,:,j)*soil%zse(j)*um1%RHO_WATER
    STHF_TILE(:,:,j) = UNPACK(REAL(ssnow%wbice(:,j)), um1%L_TILE_PTS, miss)
    SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
    TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)

    DO N=1,um1%NTILES
      DO K=1,um1%TILE_PTS(N)
        I = um1%TILE_INDEX(K,N)
        IF ( SMVCST(I) > 0. ) THEN ! Exclude permanent ice - mrd
          STHF_TILE(I,N,J)= STHF_TILE(I,N,J)/SMVCST(I)
          STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) -                &
                            STHF_TILE(I,N,J) * SMVCST(I) * soil%zse(J) &
                            * um1%RHO_WATER ) / ( soil%zse(J) *        &
                            um1%RHO_WATER * SMVCST(I) )
        END IF
      END DO
    END DO

    STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
    STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)
  END DO
   
  SNOW_TILE  = UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss)
  LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass

  SNOW_RHO1L  = UNPACK(ssnow%ssdnn, um1%L_TILE_PTS, miss)
  ISNOW_FLG3L = UNPACK(ssnow%isflag, um1%L_TILE_PTS, i_miss)
  !--- unpack layered snow vars
  DO k = 1,3
    SNOW_TMP3L(:,:,k)   = UNPACK(ssnow%tggsn(:,k), um1%L_TILE_PTS, miss)
    SNOW_MASS3L(:,:,k)  = UNPACK(ssnow%smass(:,k), um1%L_TILE_PTS, miss)
    SNOW_RHO3L(:,:,k)   = UNPACK(ssnow%ssdn(:,k), um1%L_TILE_PTS, miss)
    SNOW_COND(:,:,k)    = UNPACK(ssnow%sconds(:,k),um1%L_TILE_PTS,miss)
    SNOW_DEPTH3L(:,:,k) = UNPACK(ssnow%sdepth(:,k),um1%L_TILE_PTS,miss)
  END DO
  SNAGE_TILE = UNPACK(ssnow%snage, um1%L_TILE_PTS, miss)

  CANOPY_TILE = UNPACK(canopy%cansto, um1%L_TILE_PTS, miss)
  CANOPY_GB   = SUM(um1%TILE_FRAC * CANOPY_TILE,2)


  SURF_CAB_ROFF = UNPACK(ssnow%rnof1, um1%L_TILE_PTS, miss)
  SURF_ROFF     = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)
      
  ! Don't include sub-soil drainage for lakes
  ! NB: Hard-wired type to be removed in future version
  WHERE( veg%iveg == 16 ) ssnow%rnof2 = 0.0
  
  SURF_CAB_ROFF = UNPACK(ssnow%rnof2, um1%L_TILE_PTS, miss)
  SUB_SURF_ROFF = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

  TOT_TFALL_TILE = UNPACK(canopy%through, um1%L_TILE_PTS, miss)
  TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)
      
END SUBROUTINE cable_extra_driver
