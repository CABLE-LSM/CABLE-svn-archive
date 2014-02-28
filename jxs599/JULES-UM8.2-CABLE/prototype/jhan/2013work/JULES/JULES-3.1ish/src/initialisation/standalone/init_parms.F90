#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_parms.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_parms()

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE switches, ONLY : can_model, l_aggregate

  USE ancil_info, ONLY : ssi_index, fssi, ice_fract, ice_fract_ncat,          &
                         sea_frac, sice_frac, sice_frac_ncat, sea_index,      &
                         sice_index, sice_frac, sea_frac, sice_pts_ncat,      &
                         sice_index_ncat, land_pts, nice, ntiles, sea_pts,    &
                         sice_pts, ssi_pts, tile_pts, tile_index, frac,       &
                         land_index

  USE coastal

  USE prognostics, ONLY : di, di_ncat, tstar_tile, canht_ft, lai

  USE p_s_parms, ONLY : catch, catch_snow, infil_tile, satcon, z0_tile

  USE fluxes, ONLY : tstar

  USE u_v_grid, ONLY : dtrdz_charney_grid_1

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
  INTEGER :: error, error_sum  ! Error indicator

  INTEGER :: i,j,l,n  ! Loop counters


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Calculate surface parameters.
!-----------------------------------------------------------------------------
  CALL SPARM(LAND_PTS, NTILES, CAN_MODEL, L_AGGREGATE, TILE_PTS, TILE_INDEX,  &
             FRAC, CANHT_FT, LAI, SATCON, CATCH_SNOW, CATCH, INFIL_TILE,      &
             Z0_TILE)

!-----------------------------------------------------------------------
! Set up index for sea and sea-ice
!-----------------------------------------------------------------------
  SSI_PTS = 0
  SSI_INDEX(:) = 0
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      IF ( FLANDG(i,j) < 1.0 ) THEN
        SSI_PTS = SSI_PTS + 1
        SSI_INDEX(SSI_PTS) = (j - 1) * t_i_length + i
      ENDIF
      FSSI(i,j)=1.0 - FLANDG(i,j)
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set sea ice fraction.
!-----------------------------------------------------------------------
  DO i = 1,t_i_length
    DO j = 1,t_j_length
      ICE_FRACT(I,J) = 0.0
      DI(I,J) = 0.0
      TSTAR_SICE(I,J) = 0.0
      DO N=1,NICE
        ICE_FRACT(I,J) = ICE_FRACT(I,J) + ICE_FRACT_NCAT(I,J,N)
        DI(I,J) = DI(I,J) + ICE_FRACT_NCAT(I,J,N) * DI_NCAT(I,J,N)
      ENDDO
      IF (ICE_FRACT(I,J) > 0.0) THEN
        DO N=1,NICE  !assuming nice=nice_use here
          TSTAR_SICE(I,J) = TSTAR_SICE(I,J)                                   &
                          + ICE_FRACT_NCAT(I,J,N) * TSTAR_SICE_NCAT(I,J,N) /  &
                                                                ICE_FRACT(I,J)
        END DO
      END IF
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Allocate space for sea and sea-ice indices.
!-----------------------------------------------------------------------
  error_sum = 0
  ALLOCATE(SEA_FRAC(SSI_PTS), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(SICE_FRAC(SSI_PTS), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(SICE_INDEX_NCAT(ssi_pts,nice), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(SICE_FRAC_NCAT(SSI_PTS,NICE), STAT=error )
  error_sum = error_sum + error
! Check for error.
  IF ( error_sum /= 0 )                                                       &
    CALL log_fatal("init_parms", "Error allocating sea and sea-ice arrays")

!-----------------------------------------------------------------------
! Initialise sea and sea-ice indices
!-----------------------------------------------------------------------
  SEA_PTS  = 0
  SICE_PTS = 0
  SEA_INDEX(:)  = 0
  SICE_INDEX(:) = 0
  SICE_FRAC(:) = 0.0
  SEA_FRAC(:)  = 0.0
  DO L = 1,SSI_PTS
    J = (SSI_INDEX(L) - 1) / t_i_length + 1
    I = SSI_INDEX(L) - (J - 1) * t_i_length
    IF ( SSI_INDEX(L) > 0 ) THEN
      IF ( ICE_FRACT(I,J) > 0.0 ) THEN
        SICE_PTS = SICE_PTS + 1
        SICE_INDEX(SICE_PTS) = L
        SICE_FRAC(L) = ICE_FRACT(I,J)
      ENDIF
      IF( ICE_FRACT(I,J) < 1.0 ) THEN
        SEA_PTS = SEA_PTS + 1
        SEA_INDEX(SEA_PTS) = L
        SEA_FRAC(L) = 1.0 - SICE_FRAC(L)
      ENDIF
    ENDIF
  ENDDO

  SICE_PTS_NCAT(:) = 0
  SICE_INDEX_NCAT(:,:) = 0
  SICE_FRAC_NCAT(:,:) = 0.0
  DO N = 1,NICE
    DO L = 1,SSI_PTS
      J = (SSI_INDEX(L) - 1) / t_i_length + 1
      I = SSI_INDEX(L) - (J - 1) * t_i_length
      IF ( SSI_INDEX(L) > 0 ) THEN
        IF ( ICE_FRACT_NCAT(I,J,N) > 0.0 ) THEN
          SICE_PTS_NCAT(N) = SICE_PTS_NCAT(N) + 1
          SICE_INDEX_NCAT(SICE_PTS_NCAT(N),N) = L
          SICE_FRAC_NCAT(L,N) = ICE_FRACT_NCAT(I,J,N)
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set up gridbox "prognostics".
!-----------------------------------------------------------------------

  tstar(:,:)      = 0.0
  tstar_land(:,:) = 0.0
  tstar_ssi(:,:)  = 0.0

  DO L = 1,LAND_PTS
    J = (LAND_INDEX(L) - 1) / t_i_length + 1
    I = LAND_INDEX(L) - (J - 1) * t_i_length
    IF ( L_AGGREGATE ) THEN
      TSTAR_LAND(I,J) = TSTAR_TILE(L,1)
    ELSE
      DO N = 1,NTILES
        TSTAR_LAND(I,J) = TSTAR_LAND(I,J) + FRAC(L,N) * TSTAR_TILE(L,N)
      ENDDO
    ENDIF
  ENDDO

  TSTAR_SSI(:,:) = (1.0 - ICE_FRACT(:,:)) * TSTAR_SEA(:,:)                    &
                 + ICE_FRACT(:,:) * TSTAR_SICE(:,:)
  TSTAR(:,:) = FLANDG(:,:) * TSTAR_LAND(:,:)                                  &
             + (1.0 - FLANDG(:,:)) * TSTAR_SSI(:,:)

!-----------------------------------------------------------------------------
! Set up information on U, V and T grids (assume that att grids are the same)
!-----------------------------------------------------------------------------
  DTRDZ_CHARNEY_GRID_1(:,:) = 0.0

  RETURN

END SUBROUTINE init_parms
#endif
