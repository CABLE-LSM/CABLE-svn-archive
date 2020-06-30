!
! ==============================================================================
! Purpose: adjustment of Ci-based Jmax and Vcmax to their Cc-based values
!          (accounting for a finite mesophyll conductance) using a nonlinear
!          curve fitting routine as described in Knauer et al. 2019 GCB.
!
! Called from: SUBROUTINE bgcdriver in casa_cable.F90
!
! History: Juergen Knauer July/August 2019
! ==============================================================================
!
MODULE cable_adjust_JV_gm_module

  USE cable_def_types_mod,      ONLY: dp => r_2
  USE cable_def_types_mod,      ONLY: mp, veg_parameter_type
  USE cable_data_module,        ONLY: icanopy_type, point2constants
  USE cable_abort_module,       ONLY: nc_abort
  USE netcdf
  USE minpack

  type(icanopy_type) :: C

  integer, parameter :: nrci=3000
  integer, parameter :: nrcic4=1200
  real(dp) :: gmmax25, Vcmax25Ci, Jmax25Ci, Vcmax25Cc, Jmax25Cc, k25Ci, k25Cc
  real(dp) :: Rd
  real(dp) :: Kc_ci, Ko_ci, gammastar_ci, Km_ci
  real(dp) :: Kc_cc, Ko_cc, gammastar_cc, Km_cc
  
CONTAINS 

  SUBROUTINE adjust_JV_gm(veg)  

    IMPLICIT NONE

    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg  ! vegetation parameters

    ! local variables
    LOGICAL  :: Cc_based_OK, sw ! sw = stability switch
    INTEGER  :: p,i,k,z
    INTEGER  :: kmax=20  ! maximum nr of iterations (inner loop)
    INTEGER  :: zmax=8   ! maximum nr of iterations (outer loop)
    INTEGER  :: lAn, cntr
    REAL(dp) :: vstart, v
    REAL(dp) :: Vcmax25Cct1  ! Vcmax25Cc of previous iteration
    REAL(dp) :: Vcmax_diff
    REAL(dp) :: maxdiff=0.002e-6_dp
    REAL(dp), DIMENSION(nrci)  :: An1, Ci1
    REAL(dp), DIMENSION(:), ALLOCATABLE :: An, Ci, Cc, An_Cc
    
    ! MINPACK params
    INTEGER, PARAMETER      :: N=2     ! Number of variables
    REAL(dp), DIMENSION(N)  :: X
    REAL(dp), ALLOCATABLE   :: FVEC(:)
    INTEGER                 :: info
    REAL(dp)                :: tol=0.00001_dp

    
    ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)
    
    DO p=1,mp

      Ci1          = (/(real(i,dp),i=1,nrci,1)/) / 2.0_dp * 1.0e-6_dp  ! 1-1500 umol mol-1
      Rd           = real(veg%cfrd(p) * veg%vcmax(p),dp)
      gmmax25      = real(veg%gm(p),dp)
      Vcmax25Ci    = real(veg%vcmax(p),dp)
      Jmax25Ci     = real(veg%ejmax(p),dp)
      Kc_ci        = real(C%conkc0,dp)
      Ko_ci        = real(C%conko0,dp)
      gammastar_ci = real(C%gam0,dp)
      Kc_cc        = real(C%conkc0cc,dp)
      Ko_cc        = real(C%conko0cc,dp)
      gammastar_cc = real(C%gam0cc,dp)

      Km_ci        = Kc_ci * (1.0_dp + 0.21_dp / Ko_ci)
      Km_cc        = Kc_cc * (1.0_dp + 0.21_dp / Ko_cc)

      IF (veg%frac4(p) .lt. 0.001) THEN ! not C4

        !! 1) Calculate An-Ci curve
        CALL PHOTOSYN25(Ci1,nrci,Vcmax25Ci,Jmax25Ci,Rd,Km_ci,gammastar_ci,An1)

        !! 2) Exclude negative parts of the An-Ci curve
        lAn = count(An1 .GT. 0.0_dp)

        ALLOCATE(An(lAn))
        ALLOCATE(An_Cc(lAn))
        ALLOCATE(Ci(lAn))
        ALLOCATE(Cc(lAn))
        ALLOCATE(FVEC(lAn))

        An = pack(An1, An1 > 0.0_dp)
        Ci = pack(Ci1, An1 > 0.0_dp)
        !cntr = 0
        !DO i = 1, 3000
        !  IF (An1(i) .GT. 0.0_dp) THEN
        !    cntr = cntr + 1
        !    An(cntr) = An1(i)
        !    Ci(cntr) = Ci1(i)
        !  END IF
        !END DO

        Cc_based_OK = .FALSE.
        z = 0
        !! 3) calculate Cc based on gm and An
        DO WHILE(.NOT. Cc_based_OK .AND. z < zmax) ! if it iterates more than once, check gm and Vcmax, Jmax

           z = z + 1
           k = 0
           X(:) = [Vcmax25Ci,Jmax25Ci]
           sw = .FALSE.
           vstart = 1.0_dp
           Vcmax25Cct1 = Vcmax25Ci
           Vcmax_diff = 1.0e-6_dp
           An_Cc = An

           DO WHILE (Vcmax_diff > maxdiff .AND. k < kmax)

              k = k + 1
              Cc = Ci - An_Cc / gmmax25

              CALL LMDIF1(PHOTOSYN25_f,lAn,N,X,FVEC,tol,info,An_Cc,Cc,Rd,Km_cc,gammastar_cc)
              Vcmax25Cc = X(1)
              Jmax25Cc  = X(2)

              Vcmax_diff = ABS(Vcmax25Cc - Vcmax25Cct1)
              Vcmax25Cct1 = Vcmax25Cc

              CALL PHOTOSYN25(Cc,lAn,Vcmax25Cc,Jmax25Cc,Rd,Km_cc,gammastar_cc,An_Cc)

              ! safety switch ensuring stability
              IF (MINVAL(An_Cc) < 0.0_dp .AND. (.NOT. sw)) THEN
                 sw = .TRUE.
                 v  = vstart
              ENDIF

              IF (sw) THEN
                 v = MAX(v - (vstart/(0.8_dp*kmax)),0.0_dp)
                 An_Cc = v * An + (1.0_dp-v) * An_Cc
              ENDIF
           END DO   
           !! Avoid unrealistic Vcmax and Jmax values
           IF (Vcmax25Cc < 0.9_dp*Vcmax25Ci .OR. Vcmax25Cc > 2.5_dp*Vcmax25Ci &
               .OR. Jmax25Cc < 0.9_dp*Jmax25Ci .OR. Jmax25Cc > 1.5_dp*Jmax25Ci) THEN
              gmmax25 = 1.2_dp * gmmax25  ! If no solution, try again with higher gmmax25          
           ELSE       
              Cc_based_OK = .TRUE.
              veg%vcmaxcc(p) = real(Vcmax25Cc)
              veg%ejmaxcc(p) = real(Jmax25Cc)
           ENDIF

        END DO
     
        DEALLOCATE(An)
        DEALLOCATE(An_Cc)
        DEALLOCATE(Ci)
        DEALLOCATE(Cc)
        DEALLOCATE(FVEC)
        
   if (p == 1) then
      write(89,*) 'gmmax25:', gmmax25
      write(89,*) 'Vcmax25Ci:', Vcmax25Ci
      write(89,*) 'Rd:', Rd
      write(89,*) 'veg%vcmaxcc(p):', veg%vcmaxcc(p)
      write(89,*) 'veg%ejmaxcc(p):', veg%ejmaxcc(p)
   endif
   
      ELSE  ! C4 (Vcmax and Jmax do not change with gm in C4 plants)

        veg%vcmaxcc(p) = real(Vcmax25Ci)
        veg%ejmaxcc(p) = real(Jmax25Ci)
         
      ENDIF ! C4 flag
    END DO ! tile loop
          
  END SUBROUTINE adjust_JV_gm


  ! Function to use within LMDIF1
  SUBROUTINE PHOTOSYN25_f(M,N,X,FVEC,IFLAG,Anx,Cix,Rd,Km,gammastar)

    INTEGER,                INTENT(IN)    :: M,N,IFLAG
    REAL(dp), DIMENSION(N), INTENT(INOUT) :: X
    REAL(dp), DIMENSION(M), INTENT(OUT)   :: FVEC
    REAL(dp), DIMENSION(M), INTENT(IN)    :: Anx,Cix
    REAL(dp),               INTENT(IN)    :: Rd,Km,gammastar
    ! local
    REAL(dp), DIMENSION(M) :: Ac, Aj

    
    Ac = (X(1) * (Cix - gammastar) / (Cix + Km))
    Aj = (X(2) * (Cix - gammastar) / 4.0 / (Cix + 2.0 * gammastar))
    
    ! avoid discontinuity (e.g. Duursma 2015, PLOS ONE)
    FVEC = Anx  - ( (Ac + Aj - SQRT((Ac + Aj)**2 - 4.0*0.99999999_dp*Ac*Aj)) / &
                    (2.0*0.99999999_dp) - Rd  &
                  )
  
  END SUBROUTINE PHOTOSYN25_f


  
  ! Function to calculate An-Ci curve under standard conditions
  SUBROUTINE PHOTOSYN25(Ciz,nrci,Vcmax25,Jmax25,Rd,Km,gammastar,Anz)

    INTEGER,                   INTENT(IN)  :: nrci
    REAL(dp), DIMENSION(nrci), INTENT(IN)  :: Ciz
    REAL(dp),                  INTENT(IN)  :: Vcmax25,Jmax25,Rd,Km,gammastar
    REAL(dp), DIMENSION(nrci), INTENT(OUT) :: Anz
    ! local
    REAL(dp), DIMENSION(nrci)  :: Wc,We


    ! Rubisco-limited
    Wc =  Vcmax25 * (Ciz - gammastar) / (Ciz + Km)

    ! RuBP regeneration-limited
    We =  Jmax25 * (Ciz - gammastar) / 4.0 / (Ciz + 2.0 * gammastar)

    ! Net photosynthesis
    Anz = Min(Wc,We) - Rd

  END SUBROUTINE PHOTOSYN25



  
  
  SUBROUTINE read_gm_LUT(gm_LUT_file,veg)
    ! Read lookup table needed for parameter conversion of photosynthetic parameters
    ! from Ci- to Cc-based values (the latter considering gm explicitly)
    implicit none
    
    character(len=200), intent(in) :: gm_LUT_file
    TYPE(veg_parameter_type), intent(inout) :: veg  ! vegetation parameters
    
    ! local
    integer  :: ncid_gmlut                      ! netcdf ID
    integer  :: ok                              ! netcdf error status
    integer  :: gm_dimid, vcmax_dimid, Rd_dimid ! dimension IDs
    integer  :: gm_len, vcmax_len, Rd_len       ! dimensions of LUT
    integer  :: vcmax_id, jmax_id

    ok = nf90_open(trim(gm_LUT_file),0,ncid_gmlut)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error opening gm lookup table.')
    ok = nf90_inq_dimid(ncid_gmlut,'gm',gm_dimid)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring dimension gm from LUT.')
    ok = nf90_inq_dimid(ncid_gmlut,'Vcmax_Ci',vcmax_dimid)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring dimension Vcmax_Ci from LUT.')
    ok = nf90_inq_dimid(ncid_gmlut,'Rd',Rd_dimid)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring dimension Rd from LUT.')
    
    ok = nf90_inquire_dimension(ncid_gmlut,gm_dimid,len=gm_len)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring length of dimension gm from LUT.')
    ok = nf90_inquire_dimension(ncid_gmlut,vcmax_dimid,len=vcmax_len)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring length of dimension Vcmax_Ci from LUT.')
    ok = nf90_inquire_dimension(ncid_gmlut,Rd_dimid,len=Rd_len)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring length of dimension Rd from LUT.')

    write(*,*) 'gm LUT dimensions:'
    write(*,*) 'gm_len:', gm_len
    write(*,*) 'vcmax_len:', vcmax_len
    write(*,*) 'Rd_len:', Rd_len

    ! allocate variables in veg structure
    allocate(veg%LUT_VcmaxJmax(2,Rd_len,vcmax_len,gm_len))
    allocate(veg%LUT_gm(gm_len))
    allocate(veg%LUT_vcmax(vcmax_len))
    allocate(veg%LUT_Rd(Rd_len))
    
    ok = nf90_inq_varid(ncid_gmlut,'Vcmax',vcmax_id)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring variable Vcmax from LUT.')
    ok = nf90_inq_varid(ncid_gmlut,'Jmax',jmax_id)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error inquiring variable Jmax from LUT.')

    ok = nf90_get_var(ncid_gmlut,vcmax_id,veg%LUT_VcmaxJmax(1,:,:,:))
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error getting variable Vcmax from LUT.')
    ok = nf90_get_var(ncid_gmlut,jmax_id,veg%LUT_VcmaxJmax(2,:,:,:))
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error getting variable Jmax from LUT.')
    ok = nf90_get_var(ncid_gmlut,gm_dimid,veg%LUT_gm)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error getting dimension gm from LUT.')
    ok = nf90_get_var(ncid_gmlut,vcmax_dimid,veg%LUT_vcmax)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error getting dimension Vcmax_Ci from LUT.')
    ok = nf90_get_var(ncid_gmlut,Rd_dimid,veg%LUT_Rd)
    if (ok /= NF90_NOERR) call nc_abort(ok,'Error getting dimension Rd from LUT.')

    ! convert values to those used in CABLE
    veg%LUT_VcmaxJmax = veg%LUT_VcmaxJmax * 1.0e-06
    veg%LUT_vcmax     = veg%LUT_vcmax     * 1.0e-06
    veg%LUT_Rd        = veg%LUT_Rd        * 1.0e-06

    veg%is_read_gmLUT = .true.

  END SUBROUTINE read_gm_LUT



     

     
  
  SUBROUTINE find_Vcmax_Jmax_LUT(veg,p)
    ! JK: note that LUT is single precision at the moment
    
    implicit none

    type (veg_parameter_type), intent(inout) :: veg      ! vegetation parameters
    integer :: p             ! veg type (tile)

    ! local
    logical :: val_ok        ! check for NAs
    integer :: maxit, i      ! maximum nr of iterations in while loop, loop counter
    integer :: igm, ivc, ird ! indices for LUT

    
    if (veg%frac4(p) .lt. 0.001) then ! not C4
       ! determine current Ci-based values
       Rd        = veg%cfrd(p) * veg%vcmax(p)
       gmmax25   = veg%gm(p)
       Vcmax25Ci = veg%vcmax(p)  ! LUT assumes a given Jmax/Vcmax ratio! see details in nc LUT

       ! determine right indices of LUT
       igm = minloc(abs(gmmax25 - veg%LUT_gm),1)
       ivc = minloc(abs(Vcmax25Ci - veg%LUT_vcmax),1)
       ird = minloc(abs(Rd - veg%LUT_Rd),1)

       i = 0
       maxit = 50
       val_ok = .false. 

       do while (.not. val_ok .and. i .le. maxit)
          i = i + 1
          veg%vcmaxcc(p) = veg%LUT_VcmaxJmax(1,ird,ivc,igm)
          veg%ejmaxcc(p) = veg%LUT_VcmaxJmax(2,ird,ivc,igm)

          ! check for implausible parameter combinations that result in NAs
          if (veg%vcmaxcc(p) .gt. 0.0 .and. veg%ejmaxcc(p) .gt. 0.0) then
             val_ok = .true.  
          else
write(84,*) 'unrealistic Vcmax_ci:', veg%LUT_vcmax(ivc)
write(84,*) 'iteration:', i
             igm = igm + 1
          endif   
       end do
         
       
    else ! C4 (Vcmax and Jmax do not change with gm in C4 plants)

       veg%vcmaxcc(p) = veg%vcmax(p)
       veg%ejmaxcc(p) = veg%ejmax(p)
       
    endif  


    if (p == 1) then
  
      write(90,*) 'Rd:', Rd
      write(90,*) 'gmmax25:', gmmax25
      write(90,*) 'Vcmax25Ci:', Vcmax25Ci

      write(90,*) 'veg%vcmaxcc(p):', veg%vcmaxcc(p)
      write(90,*) 'veg%ejmaxcc(p):', veg%ejmaxcc(p)
      
    endif 

    
  END SUBROUTINE find_Vcmax_Jmax_LUT


  

  ! conversion of k Parameter in Collatz et al. 1992 from implicit gm model
  ! to explicit gm model.
  Subroutine adjust_k_Collatz(veg,p)

    implicit none

    type (veg_parameter_type), intent(inout) :: veg   ! vegetation parameters
    integer, intent(in) :: p   ! vegetation type
    
    ! local
    integer  :: i,k
    integer  :: kmax=1000
    integer  :: lAn
    real(dp) :: diff, diffx
    real(dp), dimension(nrcic4) :: An_Ci1, Ci1, Aj_Ci, Ae_Ci
    real(dp), dimension(:), allocatable :: An_Ci, An_Cc, Ci, Cc
    real(dp) :: kinc=0.001_dp  ! increment of k
    
    if (veg%frac4(p) .gt. 0.001) then ! C4
       Ci1       = (/(real(i,dp),i=1,nrcic4,1)/) / 4.0_dp * 1.0e-6_dp 
       Rd        = real(veg%cfrd(p) * veg%vcmax(p),dp)
       gmmax25   = real(veg%gm(p),dp)
       Vcmax25Ci = real(veg%vcmax(p),dp)
       k25Ci     = real(veg%c4kci(p),dp)


       ! 1) calculate An-ci curves (no light limitation)
       Aj_Ci = Vcmax25Ci - Rd
       Ae_Ci = k25Ci * Ci1 - Rd

       An_Ci1 = Min(Aj_Ci,Ae_Ci)
    
    
       ! 2) exclude negative An values and those not limited by Ci
       lAn = count(An_Ci1 .GT. 0.0_dp .AND. An_Ci1 .EQ. Ae_Ci)

       allocate(An_Ci(lAn))
       allocate(An_Cc(lAn))
       allocate(Ci(lAn))
       allocate(Cc(lAn))

       An_Ci = pack(An_Ci1, An_Ci1 .GT. 0.0_dp .AND. An_Ci1 .EQ. Ae_Ci)
       Ci    = pack(Ci1, An_Ci1 .GT. 0.0_dp .AND. An_Ci1 .EQ. Ae_Ci)

    
       ! 3) calculate Cc
       Cc = Ci - An_Ci / gmmax25

    
       ! 4) fit k to Cc-based model (a poor man's optimisation...)
       k     = 0
       diffx = 1.0e6_dp
       diff  = 0.0_dp
       k25Cc = k25Ci
       Do while (diff < diffx .AND. k < kmax)
          if (k > 0) then
            diffx = diff
          endif
          An_Cc = k25Cc * Cc - Rd
          diff = sqrt(sum((An_Cc - An_Ci)**2)/lAn)
          k25Cc = k25Cc + kinc
          k = k + 1
       End Do

       veg%c4kcc(p) = real(k25Cc - 2.0*kinc) ! single precision
       ! subtract 2x the increment to get the right value

    else   ! C3 (not used)
      veg%c4kcc(p) = 0.0
    endif

  End Subroutine adjust_k_Collatz



END MODULE cable_adjust_JV_gm_module
