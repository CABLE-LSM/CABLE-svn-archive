SUBROUTINE POP_IO ( POP, casamet, YEAR, ACTION, CF )


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! POP        : POP structure containing all specific parameter 
  ! casamet    : structure containing met and grid specific parameters from CASA
  ! YEAR       : Current year <YYYY>
  ! ACTION     : What do you want?
  !              "READ_RST"  : Read a restart file (will be looking for a file 
  !                            either with given name or from YEAR-1
  !              "WRITE_RST" : Write a restart file for YEAR+1
  !              "WRITE_EPI" : Write data at the end of each year 
  ! CLOSE_FILE : Flag to close file at the end of Episode (Episode only)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE netcdf
  USE POP_constants
  USE POP_types
  USE CASAVARIABLE
  USE CABLE_COMMON_MODULE
  USE cable_def_types_mod

  IMPLICIT NONE

  TYPE(POP_TYPE),  INTENT(INOUT) :: POP
  TYPE (casa_met), INTENT(IN)    :: casamet
  INTEGER       ,  INTENT(IN)    :: YEAR
  CHARACTER(LEN=9),INTENT(IN)    :: ACTION
  LOGICAL,OPTIONAL,INTENT(IN)    :: CF

  INTEGER*4            :: STATUS,i,m,p,l,land_ID,patch_ID,ndis_ID
  INTEGER*4            :: nlay_ID,hgtb_ID,ncoh_ID,t_ID
  INTEGER*4            :: nlayer_dim, ndisturb_dim, ndisturb1_dim,land_dim
  INTEGER*4            :: HEIGHT_BINS_dim,npatch2d_dim,NCOHORT_MAX_dim
  INTEGER*4            :: dID, t_dim, tx = -1, ntile
  CHARACTER(len=3)   :: typ = 'rst'
  CHARACTER          :: dum*9,fname*120, RUNPATH*100
  LOGICAL            :: CLOSE_FILE, EXISTFILE
  INTEGER*4:: np
   integer*4, parameter :: fi4=1
   integer*4, parameter :: ui4 = kind(fi4)
   integer*4, parameter :: onei4=1_ui4
   real*4, parameter :: fr4=1.
   integer*4, parameter :: ur4 = kind(fr4)
   real*4, parameter :: oneR4=real(1.0,ur4)
    
  !   ! 1 dim arrays (np)
  !   CHARACTER(len=40),DIMENSION( 2), PARAMETER :: AR0 = (/'latitude','longitude'/)

  !   ! LANDSCAPE STRUCTURE
  !   ! 2 dim arrays (np,t)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI1 = (/ 'npatch_active' /)
  !   CHARACTER(len=40),DIMENSION(11), PARAMETER :: AR1 = (/ 'cmass_sum',           &
  !        'densindiv','height_mean','height_max','basal_area','stress_mortality',  &
  !        'fire_mortality','growth','crown_cover','crown_area','crown_volume' /)
  !   ! 3 dim arrays (np,nlayer,t)
  !   CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR2 = (/ 'biomass','density',   &
  !        'hmean','hmax' /)
  !   ! 3 dim arrays (np,height_bins,t)
  !   CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR3 = (/ 'cmass_stem_bin',      &
  !        'densindiv_bin','height_bin','diameter_bin' /)
  !   ! 3 dim arrays (np,ndisturb,t)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI4 = (/ 'n_age' /)

  !   ! PATCH STRUCTURE
  !   ! 3 dim arrays (np,npatch2d,t)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI5 = (/ 'patch_id' /)
  !   CHARACTER(len=40),DIMENSION(10), PARAMETER :: AR5 = (/ 'patch_freq',          &
  !        'patch_freq_old','patch_freq_old2','patch_factor_recruit',               &
  !        'patch_biomass','patch_biomass_old','patch_biomass_old2',                &
  !        'patch_stress_mortality','patch_fire_mortality','patch_growth' /)
  !   ! 4 dim arrays (np,npatch2d,ndisturb+1,t)
  !   ! 4 dim arrays (np,npatch2d,ndisturb,t)
  !   CHARACTER(len=40),DIMENSION( 5), PARAMETER :: AI7 =                           &
  !        (/ 'patch_disturbance_interval','patch_first_disturbance_year',          &
  !        'patch_age','patch_age_old','patch_ranked_age_unique' /)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AR7 =                           &
  !        (/ 'patch_freq_ranked_age_unique' /)

  !   ! LAYER STRUCTURE
  !   ! 4 dim arrays (np,npatch2d,nlayer,t)
  !   CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR8 = (/ 'layer_biomass',       &
  !        'layer_density','layer_hmean','layer_hmax' /)
  !   CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI8 = (/ 'layer_ncohort' /)

  !   ! COHORT STRUCTURE
  !   ! 5 dim arrays (np,npatch2d,nlayer,ncohort_max,t)
  !   CHARACTER(len=40),DIMENSION( 5), PARAMETER :: AR9 = (/ 'cohort_biomass',      &
  !        'cohort_density','cohort_frac_resource_uptake','cohort_height',          &
  !        'cohort_diameter' /)
  !   CHARACTER(len=40),DIMENSION( 2), PARAMETER :: AI9 = (/'cohort_age','cohort_id'/)

  ! 1 dim arrays (np)
  CHARACTER(len=40),DIMENSION( 2) :: AR0

  ! LANDSCAPE STRUCTURE
  ! 2 dim arrays (np,t)
  CHARACTER(len=40),DIMENSION( 1) :: AI1
  CHARACTER(len=40),DIMENSION(11) :: AR1
  ! 3 dim arrays (np,nlayer,t)
  CHARACTER(len=40),DIMENSION( 4) :: AR2
  ! 3 dim arrays (np,height_bins,t)
  CHARACTER(len=40),DIMENSION( 4) :: AR3
  ! 3 dim arrays (np,ndisturb,t)
  CHARACTER(len=40),DIMENSION( 1) :: AI4

  ! PATCH STRUCTURE
  ! 3 dim arrays (np,npatch2d,t)
  CHARACTER(len=40),DIMENSION( 1) :: AI5
  CHARACTER(len=40),DIMENSION( 8) :: AR5
  ! 4 dim arrays (np,npatch2d,ndisturb,t)
  CHARACTER(len=40),DIMENSION( 4) :: AI7
  CHARACTER(len=40),DIMENSION( 1) :: AR7

  ! LAYER STRUCTURE
  ! 4 dim arrays (np,npatch2d,nlayer,t)
  CHARACTER(len=40),DIMENSION( 4) :: AR8
  CHARACTER(len=40),DIMENSION( 1) :: AI8

  ! COHORT STRUCTURE
  ! 5 dim arrays (np,npatch2d,nlayer,ncohort_max,t)
  CHARACTER(len=40),DIMENSION( 5) :: AR9
  CHARACTER(len=40),DIMENSION( 2) :: AI9


  INTEGER, SAVE :: VIDtime, VIDR0(SIZE(AR0)),VIDR1(SIZE(AR1)),VIDI1(SIZE(AI1)), &
       VIDR2(SIZE(AR2)),VIDR3(SIZE(AR3)),VIDI4(SIZE(AI4)),VIDR5(SIZE(AR5)),     &
       VIDI5(SIZE(AI5)),VIDI7(SIZE(AI7)),VIDR7(SIZE(AR7)),VIDR8(SIZE(AR8)),     &
       VIDI8(SIZE(AI8)),VIDR9(SIZE(AR9)),VIDI9(SIZE(AI9))
  INTEGER*4, SAVE :: FILE_ID, CNT = 0 

  ! TEMPORARY ARRAYS
  INTEGER, ALLOCATABLE :: I1(:), I2(:,:),I3(:,:,:),I4(:,:,:,:)
  REAL   , ALLOCATABLE :: R1(:), R2(:,:),R3(:,:,:),R4(:,:,:,:)

  AR0(1) = 'latitude'
  AR0(2) = 'longitude'

  AI1(1) = 'npatch_active'

  AR1(1) = 'cmass_sum'
  AR1(2) = 'densindiv'
  AR1(3) = 'height_mean'
  AR1(4) = 'height_max'
  AR1(5) = 'basal_area'
  AR1(6) = 'stress_mortality'
  AR1(7) = 'fire_mortality'
  AR1(8) = 'growth'
  AR1(9) = 'crown_cover'
  AR1(10) = 'crown_area'
  AR1(11) = 'crown_volume'

  AR2(1) = 'biomass'
  AR2(2) = 'density'
  AR2(3) = 'hmean'
  AR2(4) = 'hmax'

  AR3(1) = 'cmass_stem_bin'
  AR3(2) = 'densindiv_bin'
  AR3(3) = 'height_bin'
  AR3(4) = 'diameter_bin'

  AI4(1) = 'n_age'

  AI5(1) = 'patch_id'

  AR5(1) = 'patch_freq'
  AR5(2) = 'patch_freq_old'
  AR5(3) = 'patch_factor_recruit'
  AR5(4) = 'patch_biomass'
  AR5(5) = 'patch_biomass_old'
  AR5(6) = 'patch_stress_mortality'
  AR5(7) = 'patch_fire_mortality'
  AR5(8) = 'patch_growth'

  AI7(1) = 'patch_disturbance_interval'
  AI7(2) = 'patch_first_disturbance_year'
  AI7(3) = 'patch_age'
  AI7(4) = 'patch_ranked_age_unique'

  AR7(1) = 'patch_freq_ranked_age_unique'

  AR8(1) = 'layer_biomass'
  AR8(2) = 'layer_density'
  AR8(3) = 'layer_hmean'
  AR8(4) = 'layer_hmax'

  AI8(1) = 'layer_ncohort'

  AR9(1) = 'cohort_biomass'
  AR9(2) = 'cohort_density'
  AR9(3) = 'cohort_frac_resource_uptake'
  AR9(4) = 'cohort_height'
  AR9(5) = 'cohort_diameter'

  AI9(1) = 'cohort_age'
  AI9(2) = 'cohort_id'

  !ntile  = NINT(REAL(mp)/REAL(np)) 
ntile = mp ! warning: not currently accoutning for > 1 patch per grid cell
np = mp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF ( PRESENT( CF ) ) THEN 
     CLOSE_FILE = CF
  ELSE
     CLOSE_FILE = .FALSE.
  END IF

  ! Check for valid ACTION
  IF ( INDEX(ACTION,"WRITE_EPI") .GT. 0 ) THEN
     typ = 'out'
  ELSE IF ( INDEX(ACTION,"WRITE_RST") .GT. 0 ) THEN
     typ = 'rst'
  ELSE IF ( INDEX(ACTION,"WRITE_INI") .GT. 0 ) THEN
     typ = 'ini'
  ELSE IF ( .NOT. INDEX(ACTION,"READ_RST") .GT. 0 ) THEN
     WRITE(*,*)  "WRONG ACTION:'",TRIM(ACTION),"' in call to pop_bios_io!"
     STOP -1
  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! WRITE POP VALUES TO OUTPUT FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF ( INDEX(ACTION,"WRITE") .GT. 0 ) THEN

     CNT = CNT + onei4

     IF ( CNT .EQ. onei4 ) THEN
        ! Get File-Name
        IF( typ .EQ. 'rst' ) THEN
           WRITE( dum, FMT="(I4)")YEAR
           fname = TRIM(cable_user%POP_rst)//'/'//TRIM(dum)//'_pop_'//TRIM(cable_user%RunIDEN)//'_'//typ//'.nc'
        ELSE 
           fname = TRIM(cable_user%POP_rst)//'/pop_'//TRIM(cable_user%RunIDEN)//'_'//typ//'.nc'
        ENDIF

        write(*,*) "pop_io", fname
        ! Create NetCDF file:
        STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        ! Put the file in define mode:
        STATUS = NF90_redef(FILE_ID)

        ! GLOBAL ATTRIBUTES
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle" , icycle             )
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Year"   , YEAR               )
        STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden", CABLE_USER%RunIden )

        ! Define dimensions:
        STATUS = NF90_def_dim(FILE_ID, 'land'       , np         , land_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'NPATCH2D'   , NPATCH2D   , patch_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'NDISTURB'   , NDISTURB   , ndis_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'NLAYER'     , NLAYER     , nlay_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'HEIGHT_BINS', HEIGHT_BINS, hgtb_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'NCOHORT_MAX', NCOHORT_MAX, ncoh_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_def_dim(FILE_ID, 'time'   ,  NF90_UNLIMITED, t_ID    )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        ! Define variables
        STATUS = NF90_def_var(1,'Time' ,NF90_INT,(/1/),VIDtime )
        !STATUS = NF90_def_var(FILE_ID,'Time' ,NF90_INT,(/t_ID/),VIDtime )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        DO i = 1, SIZE(AR0)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR0(i)), nf90_double,(/land_ID/),VIDR0(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AI1)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)), NF90_INT  ,(/land_ID,t_ID/),VIDI1(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR1)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR1(i)), nf90_double,(/land_ID,t_ID/),VIDR1(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR2)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR2(i)), nf90_double,(/land_ID,nlay_ID,t_ID/),VIDR2(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR3)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR3(i)), nf90_double,(/land_ID,hgtb_ID,t_ID/),VIDR3(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AI4)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI4(i)), NF90_INT  ,(/land_ID,ndis_ID,t_ID/),VIDI4(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AI5)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI5(i)), NF90_INT  ,(/land_ID,patch_ID,t_ID/),VIDI5(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR5)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR5(i)), nf90_double,(/land_ID,patch_ID,t_ID/),VIDR5(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AI7)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI7(i)), NF90_INT,   &
                (/land_ID,patch_ID,ndis_ID,t_ID/),VIDI7(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR7)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR7(i)), nf90_double, &
                (/land_ID,patch_ID,ndis_ID,t_ID/),VIDR7(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AI8)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI8(i)), NF90_INT  , &
                (/land_ID,patch_ID,nlay_ID,t_ID/),VIDI8(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR8)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR8(i)), nf90_double, &
                (/land_ID,patch_ID,nlay_ID,t_ID/),VIDR8(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AI9)
           STATUS = NF90_def_var(FILE_ID,TRIM(AI9(i)), NF90_INT  , &
                (/land_ID,patch_ID,nlay_ID,ncoh_ID,t_ID/),VIDI9(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO
        DO i = 1, SIZE(AR9)
           STATUS = NF90_def_var(FILE_ID,TRIM(AR9(i)), nf90_double, &
                (/land_ID,patch_ID,nlay_ID,ncoh_ID,t_ID/),VIDR9(i))
           IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        END DO

        ! End define mode:
        STATUS = NF90_enddef(FILE_ID)
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

        ! PUT LAT / LON ( np )
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR0(1), REAL(casamet%lat(1:mp:ntile)),&
             start=(/ onei4 /), count=(/ np /)  )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        STATUS = NF90_PUT_VAR(FILE_ID, VIDR0(2), REAL(casamet%lon(1:mp:ntile)),&
             start=(/ onei4 /), count=(/ np /)  )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     END IF

     ! WRITE CURRENT STATE
     ! TIME  ( t )
     STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, YEAR, start=(/ CNT /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     ! PUT 2D VARS ( np, t )
     STATUS = NF90_PUT_VAR(FILE_ID, VIDI1( 1), POP%pop_grid(:)%npatch_active,      &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 1), POP%pop_grid(:)%cmass_sum,          &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 2), POP%pop_grid(:)%densindiv,          &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 3), POP%pop_grid(:)%height_mean,        &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 4), POP%pop_grid(:)%height_max,         &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 5), POP%pop_grid(:)%basal_area,         &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 6), POP%pop_grid(:)%stress_mortality,   &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 7), POP%pop_grid(:)%fire_mortality,     &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 8), POP%pop_grid(:)%growth,             &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 9), POP%pop_grid(:)%crown_cover,        &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(10), POP%pop_grid(:)%crown_area,         &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(11), POP%pop_grid(:)%crown_volume,       &
          start=(/ onei4, CNT /), count=(/ np, onei4 /) )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

     ! PUT 3D VARS ( np,nlayer, t )
     MPS:DO m = onei4, np
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 1), POP%pop_grid(m)%biomass,         &
             start=(/ m, onei4, CNT /), count=(/ onei4, nlayer, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 2), POP%pop_grid(m)%density,         &
             start=(/ m, onei4, CNT /), count=(/ onei4, nlayer, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 3), POP%pop_grid(m)%hmean,           &
             start=(/ m, onei4, CNT /), count=(/ onei4, nlayer, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 4), POP%pop_grid(m)%hmax,            &
             start=(/ m, onei4, CNT /), count=(/ onei4, nlayer, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        ! PUT 3D VARS ( np,height_bins, t )
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 1), POP%pop_grid(m)%cmass_stem_bin,  &
             start=(/ m, onei4, CNT /), count=(/ onei4, HEIGHT_BINS, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 2), POP%pop_grid(m)%densindiv_bin,   &
             start=(/ m, onei4, CNT /), count=(/ onei4, HEIGHT_BINS, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 3), POP%pop_grid(m)%height_bin,      &
             start=(/ m, onei4, CNT /), count=(/ onei4, HEIGHT_BINS, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 4), POP%pop_grid(m)%diameter_bin,    &
             start=(/ m, onei4, CNT /), count=(/ onei4, HEIGHT_BINS, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        ! PUT 3D VARS ( np,ndisturb, t )
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI4( 1), POP%pop_grid(m)%n_age,           &
             start=(/ m, onei4, CNT /), count=(/ onei4, ndisturb, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        ! PATCH STRUCTURE
        ! PUT 3D VARS ( np,npatch2d, t )
        STATUS = NF90_PUT_VAR(FILE_ID, VIDI5( 1), POP%pop_grid(m)%patch(:)%id,     &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 1), POP%pop_grid(m)%freq,            &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 2), POP%pop_grid(m)%freq_old,        &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 3), POP%pop_grid(m)%patch(:)%factor_recruit,  &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 4), POP%pop_grid(m)%patch(:)%biomass,         &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 5), POP%pop_grid(m)%patch(:)%biomass_old,     &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 6), POP%pop_grid(m)%patch(:)%stress_mortality,&
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 7), POP%pop_grid(m)%patch(:)%fire_mortality,  &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
        STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 8), POP%pop_grid(m)%patch(:)%growth,          &
             start=(/ m, onei4, CNT /), count=(/ onei4, npatch2d, onei4 /) )
        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

        PAT:DO p = onei4, npatch2d
           ! PUT 4D VARS ( np,npatch2d, ndisturb,t )
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 1), POP%pop_grid(m)%patch(p)%disturbance_interval,&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, NDISTURB, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 2), POP%pop_grid(m)%patch(p)%first_disturbance_year,&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, NDISTURB, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 3), POP%pop_grid(m)%patch(p)%age,&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, NDISTURB, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 4), POP%pop_grid(m)%ranked_age_unique(p,:),&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, NDISTURB, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           STATUS = NF90_PUT_VAR(FILE_ID, VIDR7( 1), POP%pop_grid(m)%freq_ranked_age_unique(p,:),&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, NDISTURB, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           ! LAYER STRUCTURE
           ! PUT 4D VARS ( np,npatch2d, nlayer,t )
           STATUS = NF90_PUT_VAR(FILE_ID, VIDI8( 1), POP%pop_grid(m)%patch(p)%layer(:)%ncohort,&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, nlayer, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 1), POP%pop_grid(m)%patch(p)%layer(:)%biomass,&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, nlayer, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 2), POP%pop_grid(m)%patch(p)%layer(:)%density,&
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, nlayer, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 3), POP%pop_grid(m)%patch(p)%layer(:)%hmean,  &
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, nlayer, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
           STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 4), POP%pop_grid(m)%patch(p)%layer(:)%hmax,   &
                start=(/ m, p, onei4, CNT /), count=(/ onei4, onei4, nlayer, onei4 /) )
           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           LAY:DO l = onei4, nlayer
              ! COHORT STRUCTURE
              ! PUT 5D VARS ( np,npatch2d, nlayer,ncohort_max,t )
              STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( 1), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%age,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
              STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( 2), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%id,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
              STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 1), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%biomass,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
              STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 2), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%density,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
              STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 3), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_resource_uptake,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
              STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 4), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%height,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
              STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 5), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%diameter,&
                   start=(/ m, p, l, onei4, CNT /), count=(/ onei4, onei4, onei4, ncohort_max, onei4 /) )
              IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

           END DO LAY
        END DO PAT
     END DO MPS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! READ POP VALUES AS RESTART VALUES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELSE IF ( INDEX(ACTION,'READ') .GT. 0 ) THEN

     WRITE( dum, FMT="(I4)")YEAR-1
     fname = TRIM(cable_user%POP_rst)//'/'//TRIM(dum)//'_pop_'//TRIM(cable_user%RunIDEN)//'_'//typ//'.nc'
     INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
     ! If suitable restart-file, try ini-restart
     IF ( .NOT. EXISTFILE ) THEN
        WRITE(*,*) "Restart file not found: ",TRIM(fname)
        WRITE(*,*) "Looking for initialization file..."        
        fname = TRIM(cable_user%POP_rst)//'/'//'pop_'//TRIM(cable_user%RunIDEN)//'_ini.nc'
        INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
        IF (.NOT. EXISTFILE) THEN
           WRITE(*,*) " No ini-restart file found either! ", TRIM(fname)
           STOP -1
        ENDIF
     END IF


     STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
     IF (STATUS /= NF90_noerr)THEN
        WRITE(*,*)"Error opening file (pop_bios_io.f90) ",TRIM(fname)
        CALL handle_err(STATUS)
     ENDIF
     ! DIMS
     STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NPATCH2D', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=npatch2d_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NDISTURB', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NDISTURB_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NLAYER', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NLAYER_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'HEIGHT_BINS', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=HEIGHT_BINS_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_DIMID( FILE_ID, 'NCOHORT_MAX', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=NCOHORT_MAX_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     IF ( land_dim .NE. np .OR.  npatch2d_dim .NE. NPATCH2D .OR.  &
          HEIGHT_BINS_dim .NE. HEIGHT_BINS .OR. NCOHORT_MAX_dim .NE. NCOHORT_MAX &
          .OR. NLAYER_dim .NE. NLAYER .OR. NDISTURB_dim .NE. NDISTURB ) THEN
        WRITE(*,*)"Dimension misfit in pop_bios_io.f90!"
        WRITE(*,*)"Restart file  | Current Run"
        WRITE(*,*)"# points   ",land_dim,"     ",np
        WRITE(*,*)"# patches  ",NPATCH2D_dim,"     ",NPATCH2D
        WRITE(*,*)"# HGT_BINS ",HEIGHT_BINS_dim,"     ",HEIGHT_BINS
        WRITE(*,*)"NCOHORT_MAX",NCOHORT_MAX_dim,"     ",NCOHORT_MAX
        WRITE(*,*)"# NLAYER   ",NLAYER_dim,"     ",NLAYER
        WRITE(*,*)"# NDISTURB ",NDISTURB_dim,"     ",NDISTURB
        STOP
     ENDIF

     ! TIME
     STATUS = NF90_INQ_DIMID( FILE_ID, 'time', dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=t_dim )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     STATUS = NF90_INQ_VARID( FILE_ID, 'Time', t_ID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     ALLOCATE( I1 (t_dim) )
     STATUS = NF90_GET_VAR( FILE_ID, t_ID, I1 )
     IF ( STATUS /= NF90_noerr ) CALL handle_err(STATUS)
     DO i = 1, t_dim
        IF ( YEAR .EQ. I1(i)+1 ) THEN ! DATA FROM PRECEDING YEAR !
           tx = i
           EXIT
        END IF
     END DO
     DEALLOCATE( I1 )

     IF ( tx .LE. 0 ) THEN
        WRITE(*,*) 'FILE '//TRIM(fname)//" doesn't contain data for ",YEAR
        WRITE(*,*) 'Resetting  tx to 1!'
        tx = 1
        IF ( typ .NE. "ini" ) THEN
           WRITE(*,*) "Wrong date in input pop restart-file! ",TRIM(fname)
        ENDIF
     ENDIF

     ! CHECK LAT 'N LON
     ALLOCATE( R1( np ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR0(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR ( FILE_ID, dID, R1 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     IF ( ANY ( casamet%lat(1:mp:ntile) .NE. R1 ) ) THEN
        WRITE(*,*)"INPUT LATs don't match casamet! pop_bios_io.f90" &
             , TRIM(fname)
        STOP
     ENDIF
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR0(2)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR ( FILE_ID, dID, R1 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     IF ( ANY ( casamet%lon(1:mp:ntile) .NE. R1 ) ) THEN
        WRITE(*,*)"INPUT LONs don't match casamet! pop_bios_io.f90" &
             , TRIM(fname)
        STOP
     ENDIF
     DEALLOCATE ( R1 )

     ! GET 1D VARS ( np )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI1(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, POP%pop_grid(:)%npatch_active, &
          start=(/ 1,tx /), count=(/ np, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     ALLOCATE ( R1( np ) )
     DO i = 1, SIZE(AR1)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR1(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R1, start=(/1,tx/), count=(/np,1/) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        SELECT CASE ( i )
        CASE( 1); POP%pop_grid(:)%cmass_sum        = R1
        CASE( 2); POP%pop_grid(:)%densindiv        = R1
        CASE( 3); POP%pop_grid(:)%height_mean      = R1
        CASE( 4); POP%pop_grid(:)%height_max       = R1
        CASE( 5); POP%pop_grid(:)%basal_area       = R1
        CASE( 6); POP%pop_grid(:)%stress_mortality = R1
        CASE( 7); POP%pop_grid(:)%fire_mortality   = R1
        CASE( 8); POP%pop_grid(:)%growth           = R1
        CASE( 9); POP%pop_grid(:)%crown_cover      = R1
        CASE(10); POP%pop_grid(:)%crown_area       = R1
        CASE(11); POP%pop_grid(:)%crown_volume     = R1
        CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
        END SELECT
     END DO
     DEALLOCATE ( R1 )

     ! GET 2D VARS ( np,nlayer )
     ALLOCATE( R2( np, nlayer ) )
     DO i = 1, SIZE(AR2)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR2(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R2, &
             start=(/ 1, 1, tx /), count=(/ np, nlayer, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO m = 1, np
           SELECT CASE ( i )
           CASE( 1); POP%pop_grid(m)%biomass = R2(m,:)
           CASE( 2); POP%pop_grid(m)%density = R2(m,:)
           CASE( 3); POP%pop_grid(m)%hmean   = R2(m,:)
           CASE( 4); POP%pop_grid(m)%hmax    = R2(m,:)
           CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
     END DO
     DEALLOCATE( R2 )

     ! GET 2D VARS ( np,height_bins )
     ALLOCATE( R2( np, height_bins ) )
     DO i = 1, SIZE(AR3)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR3(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R2, &
             start=(/ 1, 1, tx /), count=(/ np, height_bins, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO m = 1, np
           SELECT CASE ( i )
           CASE( 1); POP%pop_grid(m)%cmass_stem_bin = R2(m,:)
           CASE( 2); POP%pop_grid(m)%densindiv_bin  = R2(m,:)
           CASE( 3); POP%pop_grid(m)%height_bin     = R2(m,:)
           CASE( 4); POP%pop_grid(m)%diameter_bin   = R2(m,:)
           CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
     END DO
     DEALLOCATE( R2 )

     ! GET 2D VARS ( np,ndisturb)
     ALLOCATE ( I2( np, ndisturb ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI4(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, I2, &
          start=(/ 1, 1, tx /), count=(/ np, ndisturb, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     DO m = 1, np
        POP%pop_grid(m)%n_age = I2( m,: )
     END DO
     DEALLOCATE( I2 )

     ! GET 2D VARS ( np,npatch2d)
     ALLOCATE ( I2( np,npatch2d ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI5(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, I2, &
          start=(/ 1, 1, tx /), count=(/ np, npatch2d, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     DO m = 1, np
        POP%pop_grid(m)%patch(:)%id = I2( m,: )
     END DO
     DEALLOCATE( I2 )

     ! PATCH STRUCTURE
     ! GET 2D VARS ( np,npatch2d )
     ALLOCATE( R2( np, npatch2d ) )
     DO i = 1, SIZE(AR5)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR5(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R2, &
             start=(/ 1, 1, tx /), count=(/ np, npatch2d, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO m = 1, np
           SELECT CASE ( i )
           CASE( 1); POP%pop_grid(m)%freq                      = R2(m,:)
           CASE( 2); POP%pop_grid(m)%freq_old                  = R2(m,:)
           CASE( 3); POP%pop_grid(m)%patch(:)%factor_recruit   = R2(m,:)
           CASE( 4); POP%pop_grid(m)%patch(:)%biomass          = R2(m,:)
           CASE( 5); POP%pop_grid(m)%patch(:)%biomass_old      = R2(m,:)
           CASE( 6); POP%pop_grid(m)%patch(:)%stress_mortality = R2(m,:)
           CASE( 7); POP%pop_grid(m)%patch(:)%fire_mortality   = R2(m,:)
           CASE( 8); POP%pop_grid(m)%patch(:)%growth           = R2(m,:)
           CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
           END SELECT
        END DO
     END DO
     DEALLOCATE( R2 )
     
     ! GET 3D VARS ( np,npatch2d,ndisturb )
     ALLOCATE( I3( np, npatch2d, ndisturb ) )
     DO i = 1, SIZE(AI7)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI7(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, I3, &
             start=(/ 1, 1, 1, tx /), count=(/ np, npatch2d, ndisturb, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO p = 1, npatch2d
           DO m = 1, np
              SELECT CASE ( i )
              CASE( 1); POP%pop_grid(m)%patch(p)%disturbance_interval   = I3(m,p,:)
              CASE( 2); POP%pop_grid(m)%patch(p)%first_disturbance_year = I3(m,p,:)
              CASE( 3); POP%pop_grid(m)%patch(p)%age                    = I3(m,p,:)
              CASE( 4); POP%pop_grid(m)%ranked_age_unique(p,:)          = I3(m,p,:)
              CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE( I3 )

     ALLOCATE( R3( np, npatch2d,ndisturb ) )
     STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR7(1)), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR  ( FILE_ID, dID, R3, &
          start=(/ 1, 1, 1, tx /), count=(/ np, npatch2d, ndisturb, 1 /) )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     DO m = 1, np
        POP%pop_grid(m)%freq_ranked_age_unique = R3(m,:,:)
     END DO
     DEALLOCATE( R3 )

     ! LAYER STRUCTURE
     ! GET 3D VARS ( np,npatch2d,nlayer )
     ALLOCATE( I3( np, npatch2d, nlayer ) )
     DO i = 1, SIZE(AI8)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI8(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, I3, &
             start=(/ 1, 1, 1, tx /), count=(/ np, npatch2d, nlayer, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO p = 1, npatch2d
           DO m = 1, np
              SELECT CASE ( i )
              CASE( 1); POP%pop_grid(m)%patch(p)%layer(:)%ncohort = I3(m,p,:)
              CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE( I3 )

     ALLOCATE( R3( np, npatch2d, nlayer ) )
     DO i = 1, SIZE(AR8)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR8(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R3, &
             start=(/ 1, 1, 1, tx /), count=(/ np, npatch2d, nlayer, 1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO p = 1, npatch2d
           DO m = 1, np
              SELECT CASE ( i )
              CASE( 1); POP%pop_grid(m)%patch(p)%layer(:)%biomass = R3(m,p,:)
              CASE( 2); POP%pop_grid(m)%patch(p)%layer(:)%density = R3(m,p,:)
              CASE( 3); POP%pop_grid(m)%patch(p)%layer(:)%hmean   = R3(m,p,:)
              CASE( 4); POP%pop_grid(m)%patch(p)%layer(:)%hmax    = R3(m,p,:)
              CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
              END SELECT
           END DO
        END DO
     END DO
     DEALLOCATE( R3 )

     ! COHORT STRUCTURE
     ! GET 4D VARS ( np,npatch2d,nlayer,ncohort_max )
     ALLOCATE( I4( np, npatch2d, nlayer, ncohort_max ) )
     DO i = 1, SIZE(AI9)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AI9(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, I4, start=(/ 1, 1, 1, 1, tx /), &
             count=(/ np,npatch2d,nlayer,ncohort_max,1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO l = 1, nlayer
           DO p = 1, npatch2d
              DO m = 1, np
                 SELECT CASE ( i )
                 CASE( 1); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%age = I4(m,p,l,:)
                 CASE( 2); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%id  = I4(m,p,l,:)
                 CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
                 END SELECT
              END DO
           END DO
        END DO
     END DO
     DEALLOCATE( I4 )

     ALLOCATE( R4( np, npatch2d, nlayer, ncohort_max ) )
     DO i = 1, SIZE(AR9)
        STATUS = NF90_INQ_VARID( FILE_ID, TRIM(AR9(i)), dID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        STATUS = NF90_GET_VAR  ( FILE_ID, dID, R4, start=(/ 1, 1, 1, 1, tx /), &
             count=(/ np,npatch2d,nlayer,ncohort_max,1 /) )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        DO l = 1, nlayer
           DO p = 1, npatch2d
              DO m = 1, np
                 SELECT CASE ( i )
                 CASE( 1); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%biomass = R4(m,p,l,:)
                 CASE( 2); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%density = R4(m,p,l,:)
                 CASE( 3); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_resource_uptake&
                      = R4(m,p,l,:)
                 CASE( 4); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%height  = R4(m,p,l,:)
                 CASE( 5); POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%diameter= R4(m,p,l,:)
                 CASE default; STOP "Parameter not assigned in pop_bios_io.f90!"
                 END SELECT
              END DO
           END DO
        END DO
     END DO
     DEALLOCATE( R4 )

  ELSE
     WRITE(*,*) 'ACTION = ',TRIM(ACTION)
     STOP 'Please, enter either "READ" or "WRITE" when calling pop_bios_io.f90!'
  END IF

  IF ( CLOSE_FILE .OR. INDEX(ACTION,"READ").GT.0 .OR. typ .EQ. 'rst') THEN
     ! Close NetCDF file:
     STATUS = NF90_close(FILE_ID)
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     ! Reset CNT for next Write
     CNT = 0
  ENDIF

END SUBROUTINE POP_IO

