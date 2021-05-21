#if ! defined(MO_GRIB)

      !**********************************************************************
      ! Modified version of coder.F90 from UM version 7.8 to suit
      ! version 7.3
      ! Change made by SAW 20111201
      !**********************************************************************

      ! *****************************COPYRIGHT*******************************
      ! (C) Crown copyright Met Office. All rights reserved.
      ! For further details please refer to the file COPYRIGHT.txt
      ! which you should have received as part of this distribution.
      ! *****************************COPYRIGHT*******************************
      ! Dummy routine to resolve externals

      ! Code Owner: See Unified Model Code Owner's HTML page
      ! This file belongs in section: Control

      SUBROUTINE coder()

#if defined(RECON)
      USE ereport_mod, ONLY: ereport
#endif

      IMPLICIT NONE

      INTEGER                       ::  icode
      CHARACTER (len=80)            ::  cmessage
      CHARACTER (len=* ), PARAMETER ::  routinename='CODER'

      cmessage = 'Routine not available.  Please check options.'
      icode = 1
#if ! defined(RECON)
      ! DEPENDS ON: ereport
#endif
      CALL ereport(routinename,icode,cmessage)

      RETURN
      END SUBROUTINE coder
#endif

