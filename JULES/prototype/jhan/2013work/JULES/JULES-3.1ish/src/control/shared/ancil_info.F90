! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing size and dimension parameters, indexing variables
! ...and more.
!
! Most of these are not required in the UM implementation
!

MODULE ancil_info

IMPLICIT NONE

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:

!-----------------------------------------------------------------------
! Define variables that are needed by both standalone and UM
!-----------------------------------------------------------------------
  INTEGER ::                                                            &
    nsmax = 0                                                           &
                      !  Maximum number of snow layers
   ,ssi_pts                                                             &
                      !  Number of sea or sea-ice points
   ,sea_pts                                                             &
                      !  Number of sea points
   ,sice_pts
                      !  Number of sea-ice points

  INTEGER, ALLOCATABLE ::                                               &
    ssi_index(:)                                                        &
                      !  index of sea and sea-ice points
   ,sea_index(:)                                                        &
                      !  index of sea points
   ,sice_index(:)                                                       &
                      !  index of sea-ice points
   ,sice_pts_ncat(:)                                                    &
                      !  Number of points for each sea-ice category
   ,sice_index_ncat(:,:)  !  index of points for each sea-ice category


  REAL, ALLOCATABLE ::                                                  &
    fssi(:,:)                                                           &
                      !  Fraction of gridbox covered by sea
!                     !  or sea-ice
   ,sea_frac(:)                                                         &
                      !  Fraction of gridbox covered by sea
!                     !  (converted to single vector array)
   ,sice_frac(:)                                                        &
                      !  Fraction of gridbox covered by sea-ice
!                     !  (converted to single vector array)
   ,sice_frac_ncat(:,:)     !  Fraction of gridbox covered by each
!                           !  sea-ice category
!                           !  (converted to single vector array)


!-----------------------------------------------------------------------
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------
#if !defined(UM_JULES)

! The following block declares variables that are being removed
! from the UM. They have not been removed from here because this
! code is never used in the UM. When the standalone code is
! supplied with an equivalent to atm_fields_bounds_mod, they
! may be mapped to the variables in that equivalent module.

  INTEGER ::                                                            &
    halo_i = 0                                                          &
                      !  Size of halo in i direction
   ,halo_j = 0                                                          &
                      !  Size of halo in j direction
   ,n_rows                                                              &
                      !  Number of rows in a v field
   ,off_x = 0                                                           &
                      !  Size of small halo in i
   ,off_y = 0                                                           &
                      !  Size of small halo in j
   ,row_length                                                          &
                      !  Number of points on a row
   ,rows              !  Number of rows in a theta field

  INTEGER ::                                                            &
    co2_dim_len                                                         &
                      !  Length of a CO2 field row
   ,co2_dim_row                                                         &
                      !  Number of CO2 field rows
   ,land_pts = 0                                                        &
                      !  No. of land points
   ,land_pts_trif                                                       &
                      !  For dimensioning land fields in TRIFFID
   ,lice_pts                                                            &
                      !  Number of land ice points
   ,npft_trif                                                           &
                      !  For dimensioning pft fields in TRIFFID
!                     !   =npft when TRIFFID on, otherwise =1
   ,ntiles                                                              &
                      !  Number of surface tiles
   ,sm_levels = 4                                                       &
                      !  Number of soil layers
   ,soil_pts                                                            &
                      !  Number of soil points
   ,nice                                                                &
                      !  Number of sea ice categories
   ,nice_use                                                            &
                      !  Number of sea ice cats used in radiation and
                      !   the explicit part of surface exchange
   ,dim_cs1                                                             &
                      !  size of second dimension in soil carbon (cs)
!                        and related respiration variables
   ,dim_cs2           !  size used for some variables that are only
!                        used with TRIFFID. If not using TRIFFID these
!                        variables are set to be smaller to save some space.

  INTEGER, ALLOCATABLE ::                                               &
    land_index(:)                                                       &
                          !  index of land points
   ,tile_index(:,:)                                                     &
                          !  indices of land points which include the
!                         !  nth surface type
   ,soil_index(:)                                                       &
                          !  index of soil points (i.e. land point
!                         !  number for each soil point)
   ,lice_index(:)                                                       &
                          !  index of land ice points (i.e. land point
!                         !  number for each land ice point)
   ,tile_pts(:)
                          !  Number of land points which include the
!                         !  nth surface type

  REAL, ALLOCATABLE ::                                                  &
    frac(:,:)                                                           &
                            !  fractional cover of each surface type
   ,z1_tq(:,:)                                                          &
                            !  height of temperature data
   ,z1_uv(:,:)                                                          &
                            !  height of wind data
   ,ice_fract(:,:)                                                      &
                            !  fraction of gridbox covered by sea-ice
!                           !  (decimal fraction)
   ,ice_fract_ncat(:,:,:)
                            !  fraction of gridbox covered by sea-ice
!                           !  on catagories


  LOGICAL, ALLOCATABLE ::                                               &
    land_mask(:,:)          !  T if land, F elsewhere

#endif

!-----------------------------------------------------------------------
! Namelist for UM only - sm_levels and nsmax are read in the namelist
! jules_model_levels in standalone
!-----------------------------------------------------------------------
  NAMELIST /jules_ancil_info/ nsmax

END MODULE ancil_info
