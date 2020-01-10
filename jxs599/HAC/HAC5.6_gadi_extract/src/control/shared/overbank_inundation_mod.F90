!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE overbank_inundation_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains river overbank inundation variables and switches
!
!  Code Owner: Please refer to ModuleLeaders.txt
!  This file belongs in section: Hydrology
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: rmdi

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Overbank inundation parameters
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm) ::                                                     &
   riv_a = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn1)
   riv_b = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn1)
   riv_c = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn2)
   riv_f = rmdi,                                                              &
        ! Parameter in Leopold & Maddock (1953: eqn2)
   coef_b = rmdi,                                                             &
        ! Parameter to calculate the bankflow discharge power-law
        ! relationship from "Flood Modeling, Prediction and Mitigation"
        ! by Şen 2018.
   exp_c = rmdi,                                                              &
        ! Parameter to calculate the bankflow discharge power-law
        ! relationship from "Flood Modeling, Prediction and Mitigation"
        ! by Şen 2018.
   ent_ratio = rmdi
        ! Entrenchment ratio (Rosgen 1994). Ratio of the width of
        ! flood-prone area to surface width of bankfull channel.
        ! The flood-prone area width is measured at the elevation that
        ! corresponds to twice the maximum depth of the bankfull channel.

LOGICAL ::                                                                    &
   l_riv_overbank = .FALSE.,                                                  &
                            ! Logical to control overbank inundation
   l_riv_hypsometry = .TRUE.,                                                 &
                            ! TRUE indicates to calculate inundated area
                            ! using hypsometry, which requires additional
                            ! ancillaries logn_mean and logn_stdev
                            ! FALSE indicates to use a simpler width scaling
                            ! method (generally only to be used for testing
                            ! when those ancillaries are not available).
   use_rosgen = .FALSE.
                            ! Modify floodplain width using the Rosgen
                            ! entrenchment ratio (n.b. only for local use)

!-----------------------------------------------------------------------------
! Array variables defined on land points updated in overbank inundation
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  frac_fplain_lp(:)
                          ! fraction of gridcell that is temporarily-
                          ! inundated open water surface (for overbank
                          ! inundation)

!-----------------------------------------------------------------------------
! Array variables defined on full 2D rivers grid (as read in from ancillary)
!----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  logn_mean(:,:),                                                             &
                          ! ln(mean(elevation - elev_min)) for each gridcell
                          ! (elevation is relative to minimum elevation)
                          !in grid box)
  logn_stdev(:,:)
                          ! ln(SD(elevation - elev_min)) for each gridcell

!-----------------------------------------------------------------------------
! Other array variables
!----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  logn_mean_rp(:),                                                            &
           ! as logn_mean but on river points
  logn_stdev_rp(:),                                                           &
           ! as logn_stdev but on river points
  qbf(:),                                                                     &
           ! bankfull discharge rate in m3/s
  dbf(:),                                                                     &
           ! channel depth when at bankfull discharge rate (m)
  wbf(:),                                                                     &
           ! channel width when at bankfull discharge rate (m)
  frac_fplain_rp(:)
           ! Fraction of gridcell predicted to have overbank inundation
           ! (on rivers grid)

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_overbank/                                                   &
  l_riv_overbank, l_riv_hypsometry, use_rosgen,                               &
  riv_a, riv_b, riv_c, riv_f, coef_b, exp_c, ent_ratio

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OVERBANK_INUNDATION_MOD'

CONTAINS

SUBROUTINE check_jules_overbank()

USE ereport_mod, ONLY: ereport

USE jules_rivers_mod, ONLY: l_rivers

USE jules_print_mgr, ONLY: jules_print

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_OVERBANK switches for consistency
!
! Current Code Owner: Toby Marthews (CEH)
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer
INTEGER :: errcode

WRITE(lineBuffer,*) ' l_rivers ', l_rivers
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*) ' l_riv_overbank ', l_riv_overbank
CALL jules_print('jules_overbank',lineBuffer)

IF ( ( .NOT. l_rivers) .AND. l_riv_overbank) THEN
  errcode = 101
  CALL ereport("check_jules_overbank", errcode,                               &
               'l_riv_overbank=T requires l_rivers=T ')
END IF

!! Further sanity checks (e.g. parameters within range) possible here....
!! see check_jules_rivers()

END SUBROUTINE check_jules_overbank

SUBROUTINE print_nlist_jules_overbank()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_JULES_OVERBANK'

CALL jules_print('overbank_inundation_mod',                                   &
   'Contents of namelist jules_overbank')

WRITE(lineBuffer,*)' l_riv_overbank = ',l_riv_overbank
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' l_riv_hypsometry = ',l_riv_hypsometry
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' use_rosgen = ',use_rosgen
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_a = ',riv_a
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_b = ',riv_b
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_c = ',riv_c
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' riv_f = ',riv_f
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' coef_b = ',coef_b
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' exp_c = ',exp_c
CALL jules_print('jules_overbank',lineBuffer)
WRITE(lineBuffer,*)' ent_ratio = ',ent_ratio
CALL jules_print('jules_overbank',lineBuffer)

CALL jules_print('overbank_inundation_mod',                                   &
   '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_overbank


#if defined(UM_JULES) && !defined(LFRIC)
SUBROUTINE read_nml_jules_overbank(unit_in)

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_OVERBANK'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_real = 7
INTEGER, PARAMETER :: n_log = 3

TYPE my_namelist
  SEQUENCE
  LOGICAL :: l_riv_overbank
  LOGICAL :: l_riv_hypsometry
  LOGICAL :: use_rosgen
  REAL(KIND=real_jlslsm) :: riv_a
  REAL(KIND=real_jlslsm) :: riv_b
  REAL(KIND=real_jlslsm) :: riv_c
  REAL(KIND=real_jlslsm) :: riv_f
  REAL(KIND=real_jlslsm) :: coef_b
  REAL(KIND=real_jlslsm) :: exp_c
  REAL(KIND=real_jlslsm) :: ent_ratio
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type,                                &
                    n_real_in = n_real, n_log_in = n_log)

IF (mype == 0) THEN

  READ(UNIT = unit_in, NML = jules_overbank, IOSTAT = ErrorStatus, IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist JULES_OVERBANK", iomessage)

  my_nml % l_riv_overbank   = l_riv_overbank
  my_nml % l_riv_hypsometry = l_riv_hypsometry
  my_nml % use_rosgen       = use_rosgen
  my_nml % riv_a            = riv_a
  my_nml % riv_b            = riv_b
  my_nml % riv_c            = riv_c
  my_nml % riv_f            = riv_f
  my_nml % coef_b           = coef_b
  my_nml % exp_c            = exp_c
  my_nml % ent_ratio        = ent_ratio

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  l_riv_overbank = my_nml % l_riv_overbank
  l_riv_hypsometry = my_nml % l_riv_hypsometry
  use_rosgen = my_nml % use_rosgen
  riv_a = my_nml % riv_a
  riv_b = my_nml % riv_b
  riv_c = my_nml % riv_c
  riv_f = my_nml % riv_f
  coef_b = my_nml % coef_b
  exp_c = my_nml % exp_c
  ent_ratio = my_nml % ent_ratio

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_jules_overbank
#endif

END MODULE overbank_inundation_mod

!-----------------------------------------------------------------------------
