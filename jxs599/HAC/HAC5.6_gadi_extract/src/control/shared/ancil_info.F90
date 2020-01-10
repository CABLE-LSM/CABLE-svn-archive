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

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:

!-----------------------------------------------------------------------
! Define variables that are needed by both standalone and UM
!-----------------------------------------------------------------------
INTEGER ::                                                                    &
  land_pts = 0                                                                &
                    !  No. of land points
  ,ssi_pts                                                                    &
                    !  Number of sea or sea-ice points
                    !  NOTE: Currently, this is set to the total number
                    !        of grid boxes
                    !        This probably needs to be looked at in
                    !        the future...
 ,sea_pts                                                                     &
                    !  Number of sea points
 ,sice_pts                                                                    &
                    !  Number of sea-ice points
 ,nsoilt = 1
                    !  Number of soil tiles, defaulting to 1

INTEGER, ALLOCATABLE ::                                                       &
 land_index(:)                                                                &
                        !  index of land points
 ,ssi_index(:)                                                                &
                    !  index of sea and sea-ice points
 ,sea_index(:)                                                                &
                    !  index of sea points
 ,sice_index(:)                                                               &
                    !  index of sea-ice points
 ,sice_pts_ncat(:)                                                            &
                    !  Number of points for each sea-ice category
 ,sice_index_ncat(:,:)                                                        &
                    !  index of points for each sea-ice category
 ,soilt_pts(:)                                                                &
                    !  Number of points for each soil tile
 ,soilt_index(:,:)
                    !  Number of points that include the nth soil tile


LOGICAL, ALLOCATABLE ::                                                       &
  l_soil_point(:),                                                            &
                        !  TRUE if a soil point, FALSE otherwise
  l_lice_point(:),                                                            &
                        !  TRUE if a land ice point, FALSE otherwise
                        !  Used as a test for ice points, replacing the
                        !  old test using sm_sat <= 0.0, which will not
                        !  work (or be rather inflexible) with tiled soils.
  l_lice_surft(:)
                        !  TRUE if a land ice (surface) tile, FALSE otherwise

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  fssi_ij(:,:)                                                                &
                    !  Fraction of gridbox covered by sea
!                     !  or sea-ice
   ,sea_frac(:)                                                               &
                      !  Fraction of gridbox covered by sea
!                     !  (converted to single vector array)
   ,sice_frac(:)                                                              &
                      !  Fraction of gridbox covered by sea-ice
!                     !  (converted to single vector array)
   ,sice_frac_ncat(:,:)                                                       &
                      !  Fraction of gridbox covered by each
                      !  sea-ice category
                      !  (converted to single vector array)
   ,frac_soilt(:,:)
                      !  Fraction of gridbox for each soil tile


!-----------------------------------------------------------------------
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------

! The following block declares variables that are being removed
! from the UM. They have not been removed from here because this
! code is never used in the UM. When the standalone code is
! supplied with an equivalent to atm_fields_bounds_mod, they
! may be mapped to the variables in that equivalent module.

INTEGER ::                                                                    &
  halo_i = 0                                                                  &
                    !  Size of halo in i direction
 ,halo_j = 0                                                                  &
                    !  Size of halo in j direction
 ,n_rows                                                                      &
                    !  Number of rows in a v field
 ,off_x = 0                                                                   &
                    !  Size of small halo in i
 ,off_y = 0
                    !  Size of small halo in j

INTEGER ::                                                                    &
  co2_dim_len                                                                 &
                    !  Length of a CO2 field row
 ,co2_dim_row                                                                 &
                    !  Number of CO2 field rows
 ,land_pts_trif                                                               &
                    !  For dimensioning land fields in TRIFFID
 ,lice_pts                                                                    &
                    !  Number of land ice points
 ,npft_trif                                                                   &
                    !  For dimensioning pft fields in TRIFFID
                    !   =npft when TRIFFID on, otherwise =1
 ,nsurft                                                                      &
                    !  Number of surface tiles
 ,soil_pts                                                                    &
                    !  Number of soil points
 ,dim_cs1                                                                     &
                    !  size of pool dimension in soil carbon (cs)
                    !  and related respiration variables
 ,dim_cslayer = 1                                                             &
                    !  size of depth dimension in soil carbon (cs)
                    !  and related respiration variables
                    !  Initialised to 1 to ensure that configuations that
                    !  do not use science connected with it are not broken
                    !  by zero-sized array dimensions
 ,dim_cs2                                                                     &
                    !  size used for some variables that are only
                    !  used with TRIFFID. If not using TRIFFID these
                    !  variables are set to be smaller to save some space
 ,row_length                                                                  &
                    !  Number of points on a row
 ,rows              !  Number of rows in a theta field

INTEGER, ALLOCATABLE ::                                                       &
 surft_index(:,:)                                                             &
                        !  indices of land points which include the
                        !  nth surface type
 ,soil_index(:)                                                               &
                        !  index of soil points (i.e. land point
                        !  number for each soil point)
 ,lice_index(:)                                                               &
                        !  index of land ice points (i.e. land point
                        !  number for each land ice point)
 ,surft_pts(:)
                        !  Number of land points which include the
                        !  nth surface type

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  frac_surft(:,:)                                                             &
                          !  fractional cover of each surface type
 ,z1_tq_ij(:,:)                                                               &
                          !  height of temperature data
 ,z1_uv_ij(:,:)                                                               &
                          !  height of wind data
 ,ice_fract_ij(:,:)                                                           &
                          !  fraction of gridbox covered by sea-ice
                         !  (decimal fraction)
 ,ice_fract_ncat_sicat(:,:,:)                                                 &
                          !  fraction of gridbox covered by sea-ice
                         !  on catagories
 ,ti_cat_sicat(:,:,:)                                                         &
                          ! sea ice surface temperature on categories
 ,pond_frac_cat_sicat(:,:,:)                                                  &
                         ! Meltpond fraction on sea ice categories
 ,pond_depth_cat_sicat(:,:,:)                                                 &
                         ! Meltpond depth on sea ice categories (m)
 ,sstfrz_ij(:,:)
                         ! Salinity-dependent sea surface freezing
                         ! temperature (K)

LOGICAL, ALLOCATABLE ::                                                       &
  land_mask(:,:)          !  T if land, F elsewhere

END MODULE ancil_info
