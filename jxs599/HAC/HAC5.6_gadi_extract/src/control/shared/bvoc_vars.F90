! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Module containing BVOC diagnostics

MODULE bvoc_vars

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  isoprene_gb(:),                                                             &
          ! Gridbox mean isoprene emission flux (kgC/m2/s)
  isoprene_pft(:,:),                                                          &
          ! Isoprene emission flux on PFTs (kgC/m2/s)
  terpene_gb(:),                                                              &
          ! Gridbox mean (mono-)terpene emission flux (kgC/m2/s)
  terpene_pft(:,:),                                                           &
          ! (Mono-)Terpene emission flux on PFTs (kgC/m2/s)
  methanol_gb(:),                                                             &
          ! Gridbox mean methanol emission flux (kgC/m2/s)
  methanol_pft(:,:),                                                          &
          ! Methanol emission flux on PFTs (kgC/m2/s)
  acetone_gb(:),                                                              &
          ! Gridbox mean acetone emission flux (kgC/m2/s)
  acetone_pft(:,:)
          ! Acetone emission flux on PFTs (kgC/m2/s)

END MODULE bvoc_vars
