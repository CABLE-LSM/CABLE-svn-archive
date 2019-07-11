#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
     SUBROUTINE OASIS4_ADVANCE_DATE64(ierr64)
!
! Description: This subroutine is a 64 bit wrapper interface 
!              around the call to OASIS4_ADVANCE_DATE32 
!              code needed in order to hide the 32 bit 
!
! History:
! UM Version        Date        Description
! ----------     -----------    ----------------------------------
!    6.4          Jan 2007      Original Code. R. Hill
!=================================================================

! DEPENDS ON: OASIS4_grid32_mod
      USE OASIS4_grid32_mod

      IMPLICIT NONE

      INTEGER :: ierr64

! DEPENDS ON: OASIS4_grid32_mod
      CALL OASIS4_ADVANCE_DATE32(ierr64)

      END SUBROUTINE OASIS4_ADVANCE_DATE64
#endif
