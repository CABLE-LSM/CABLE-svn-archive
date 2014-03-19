! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains variables used for reading in triffid data
! and initialisations


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE trif_io

  USE max_dimensions, ONLY:                                           &
    npft_max

  IMPLICIT NONE

!---------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in trif that we want to initialise).
!---------------------------------------------------------------------
  INTEGER ::                                                          &
    crop_io(npft_max)

  REAL ::                                                             &
    g_area_io(npft_max),                                              &
    g_grow_io(npft_max),                                              &
    g_root_io(npft_max),                                              &
    g_wood_io(npft_max),                                              &
    lai_max_io(npft_max),                                             &
    lai_min_io(npft_max)

!---------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!---------------------------------------------------------------------
  NAMELIST /jules_triffid/ crop_io,g_area_io,g_grow_io,g_root_io,     &
                           g_wood_io,lai_max_io,lai_min_io

CONTAINS
  SUBROUTINE print_nlist_jules_triffid()
    USE jules_print_mgr, ONLY : jules_print
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer

    CALL jules_print('trif_io', &
        'Contents of namelist jules_triffid')

    WRITE(lineBuffer,*)' crop_io = ',crop_io
    CALL jules_print('trif_io',lineBuffer)
    WRITE(lineBuffer,*)' g_area_io = ',g_area_io
    CALL jules_print('trif_io',lineBuffer)
    WRITE(lineBuffer,*)' g_grow_io = ',g_grow_io
    CALL jules_print('trif_io',lineBuffer)
    WRITE(lineBuffer,*)' g_root_io = ',g_root_io
    CALL jules_print('trif_io',lineBuffer)
    WRITE(lineBuffer,*)' g_wood_io = ',g_wood_io
    CALL jules_print('trif_io',lineBuffer)
    WRITE(lineBuffer,*)' lai_max_io = ',lai_max_io
    CALL jules_print('trif_io',lineBuffer)
    WRITE(lineBuffer,*)' lai_min_io = ',lai_min_io
    CALL jules_print('trif_io',lineBuffer)

    CALL jules_print('trif_io', &
        '- - - - - - end of namelist - - - - - -')

  END SUBROUTINE print_nlist_jules_triffid
END MODULE trif_io
