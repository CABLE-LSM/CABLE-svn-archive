!*******************************************************************************

MODULE TypeDef
  !-------------------------------------------------------------------------------
  ! PRB: 18-06-2000
  ! This module explicitly defines the sizes of variable types
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
  ! Define integer kind parameters to accommodate the range of numbers usually
  ! associated with 4, 2, and 1 byte integers.
  INTEGER,PARAMETER :: i4b = SELECTED_INT_KIND(9)
  INTEGER,PARAMETER :: i2b = SELECTED_INT_KIND(4)
  INTEGER,PARAMETER :: i1b = SELECTED_INT_KIND(2)
  ! Define single and double precision real kind parameters:
  ! * Kind(1.0)   defines sp as the machine's default size for single precision
  ! * Kind(1.0d0) defines dp as the machine's default size for double precision
  !INTEGER,PARAMETER :: sp  = KIND(1.0)
  !INTEGER,PARAMETER :: dp  = KIND(1.0d0)
  !INTEGER,PARAMETER :: sp = SELECTED_REAL_KIND(6,37)    ! IEEE single Precision Real
  INTEGER,PARAMETER :: dp  = SELECTED_REAL_KIND(12, 50) ! IEEE double Precision Real
  ! lgt is set to the default kind required for representing logical values.
  INTEGER,PARAMETER :: lgt = KIND(.TRUE.)

END MODULE TypeDef


!*******************************************************************************
MODULE POP_Constants
  USE TYPEdef, ONLY: dp, i4b

  !REAL(dp),PARAMETER:: FULTON_ALPHA=9.61 ! recruitment scalar alpha in Fulton (1991)
  REAL(dp),PARAMETER:: FULTON_ALPHA=5.6 ! recruitment scalar alpha in Fulton (1991)
  REAL(dp),PARAMETER:: DENSINDIV_MAX=2.0  !  Maximum density of individuals within a cohort indiv/m2
  REAL(dp),PARAMETER:: DENSINDIV_MIN=1e-12 !
  REAL(dp),PARAMETER:: Kbiometric=50.0 ! Constant in height-diameter relationship
  REAL(dp),PARAMETER:: WD= 300.0 ! Wood density kgC/m3
  REAL(dp),PARAMETER:: GROWTH_EFFICIENCY_MIN=0.015 ! threshold growth efficiency for enhanced mortality
  !REAL(dp),PARAMETER:: GROWTH_EFFICIENCY_MIN=0.0076 ! threshold growth efficiency for enhanced mortality
  REAL(dp),PARAMETER:: Pmort=5.0 ! exponent in mortality formula
  REAL(dp),PARAMETER:: MORT_MAX=0.3 ! upper asymptote for enhanced mortality
  REAL(dp),PARAMETER:: THETA_recruit=0.95 ! shape parameter in recruitment equation
  REAL(dp),PARAMETER:: CMASS_STEM_INIT= 1e-6 ! initial biomass kgC/m2
  REAL(dp),PARAMETER:: POWERbiomass=0.75 ! exponent for biomass in proportion to which cohorts preempt resources
  REAL(dp),PARAMETER:: POWERGrowthEfficiency = 0.75
  REAL(dp),PARAMETER:: CrowdingFactor = 0.016
  !REAL(dp),PARAMETER:: CrowdingFactor = 0.049
  REAL(dp),PARAMETER:: ALPHA_CPC = 10.0
  !REAL(dp),PARAMETER:: ALPHA_CPC = 3.33
  REAL(dp),PARAMETER:: k_allom1 = 200.0 ! crownarea =  k_allom1 * diam ** k_rp
  REAL(dp),PARAMETER:: k_rp = 1.67


  REAL(dp),PARAMETER:: Q=7.0 ! governs rate of increase of mortality with age (2=exponential)
  REAL(dp),PARAMETER:: F_MORT=0.0001 ! proportion of cohort surviving to age Amax
  REAL(dp),PARAMETER:: AMAX=250.0 ! age to which proportion F of cohort will survive in absence of resource stress
  REAL(dp),PARAMETER:: EPS=1e-12;
  INTEGER(i4b),PARAMETER :: NLAYER = 1
  INTEGER(i4b),PARAMETER :: NCOHORT_MAX = 20
  INTEGER(i4b),PARAMETER :: NDISTURB=2 ! number of disturbance regimes (1 or 2)
  INTEGER(i4b),PARAMETER :: PATCH_REPS=7 ! higher number reduces 'noise'
  INTEGER(i4b),PARAMETER :: NAGE_MAX = 6 ! number of maxium ages
  INTEGER(i4b),PARAMETER :: NPATCH=(NAGE_MAX*PATCH_REPS)**NDISTURB ! number of patches with different disturbance intervals to simulate
  INTEGER(i4b),PARAMETER :: NPATCH1D=(NAGE_MAX*PATCH_REPS)
  INTEGER(i4b),PARAMETER :: NPATCH2D= NPATCH  ! + NPATCH1D*(NDISTURB-1) ! number of patches to be simulated, including those correponding to a1>a2
  INTEGER(i4b),PARAMETER ::  HEIGHT_BINS=12 ! number of height categories to keep track of for diagnostics
  REAL(dp),PARAMETER:: BIN_POWER=1.4 ! bins have muscles
  INTEGER(i4b),PARAMETER :: TIMEBASE_FACTOR=50 ! Time base factor (to be multiplied by mean dist interval to give TIMEBASE)  for sampling disturbance probabilities from Poisson distribution
  REAL(dp),PARAMETER:: PI=3.14159265358979323846264
  INTEGER(i4b),PARAMETER :: ALLOM_SWITCH = 0 ! 0 == default; 1 = top-end allometry (requires precip as input to POPSTEP)
  INTEGER(i4b),PARAMETER :: MAX_HEIGHT_SWITCH = 2 ! 0 == binnned max height variable; 1 = continuous (needs lots of memory); 2 = binned by integer heights

END MODULE POP_Constants
!*******************************************************************************
MODULE POP_Types
  USE TYPEdef, ONLY: dp, i4b
  USE POP_Constants, ONLY: NCOHORT_MAX, NLAYER, HEIGHT_BINS, NDISTURB, NPATCH, NPATCH2D
  !USE BGCModule, ONLY: np

  TYPE Cohort
     INTEGER(i4b) :: id
     INTEGER(i4b) :: age ! cohort age
     REAL(dp)     :: biomass ! cohort biomass
     REAL(dp)     :: density ! landscape tree density (weighted mean over patches)
     REAL(dp)     :: frac_resource_uptake
     REAL(dp)     :: height
     REAL(dp)     :: diameter

  END TYPE Cohort

  TYPE Layer
     TYPE (Cohort), DIMENSION(NCOHORT_MAX) :: Cohort
     INTEGER(i4b) :: ncohort ! number of cohorts with density >0
     REAL(dp)    :: biomass ! layer biomass
     REAL(dp)    :: density ! layer tree density
     REAL(dp)     :: hmean ! layer mean tree height (weighted mean over patches)
     REAL(dp)     :: hmax  ! layer max tree height
  END TYPE Layer

  TYPE Patch
     TYPE (Layer), DIMENSION(NLAYER) :: Layer
     REAL(dp)     :: factor_recruit
     REAL(dp) :: biomass ! total biomass in patch
     REAL(dp) :: biomass_old ! total biomass in patch
     REAL(dp) :: biomass_old2 ! total biomass in patch
     REAL(dp) :: stress_mortality ! biomass lost in each patch due to stress
     REAL(dp) :: fire_mortality ! biomass lost in each patch due partial fire disturbance
     REAL(dp) :: crowding_mortality ! biomass lost to crowding mortality
     REAL(dp) :: mortality !
     REAL(dp) :: growth ! biomass growth in each patch due to stem increment
     INTEGER(i4b) :: disturbance_interval(NDISTURB)  ! prescribed disturbance(s) interval for this patch
     INTEGER(i4b) :: first_disturbance_year(NDISTURB)
     INTEGER(i4b) :: age(NDISTURB) ! number of years since last disturbance(s)
     INTEGER(i4b) :: age_old(NDISTURB) ! number of years since last disturbance(s)
     INTEGER(i4b) :: id
  END TYPE Patch

  TYPE Landscape
     TYPE (Patch), DIMENSION(NPATCH2D) :: patch
     REAL(dp), DIMENSION(NPATCH2D)     :: freq ! patch weighting
     REAL(dp), DIMENSION(NPATCH2D)     :: freq_old ! patch weighting
     REAL(dp), DIMENSION(NPATCH2D)     :: freq_old2 ! patch weighting
     REAL(dp), DIMENSION(NPATCH2D,NDISTURB)     :: freq_ranked_age_unique ! unique age weighting
     INTEGER(i4b), DIMENSION(NPATCH2D, NDISTURB)     :: ranked_age_unique ! unique age
     INTEGER(i4b), DIMENSION(NDISTURB)     :: n_age ! number of unique ages
     REAL(dp), DIMENSION(NLAYER)     :: biomass ! landscape stem biomass (weighted mean over patches)
     REAL(dp), DIMENSION(NLAYER)     :: density ! landscape tree density (weighted mean over patches)
     REAL(dp), DIMENSION(NLAYER)     :: hmean ! landscape mean treen height (weighted mean over patches)
     REAL(dp), DIMENSION(NLAYER)     :: hmax  ! landscape max tree height
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: cmass_stem_bin ! biomass by height bin
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: densindiv_bin ! density by height bin
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: height_bin ! mean height in each bin
     REAL(dp), DIMENSION(HEIGHT_BINS)     :: diameter_bin ! mean diameter in each bin
     CHARACTER(100), DIMENSION(HEIGHT_BINS) :: bin_labels ! text strings for bin bounds
     REAL(dp) :: cmass_sum ! landscape biomass
     REAL(dp) :: densindiv ! landscape density of individuals
     REAL(dp) :: height_mean
     REAL(dp) :: height_max
     REAL(dp) :: basal_area
     REAL(dp) :: stress_mortality ! (kg C m-2 y-1)
     REAL(dp) :: fire_mortality ! (kg C m-2 y-1)
     REAL(dp) :: growth
     REAL(dp) :: crown_cover
     REAL(dp) :: crown_area
     REAL(dp) :: crown_volume
     INTEGER(i4b) :: npatch_active
  END TYPE Landscape

  TYPE POP_TYPE
     TYPE(Landscape), DIMENSION(:), ALLOCATABLE :: pop_grid
     INTEGER                                    :: it_pop
  END TYPE POP_TYPE

END MODULE POP_Types
!*******************************************************************************

MODULE POPModule
  !-------------------------------------------------------------------------------
  ! * This module contains all subroutines for POP calcs at a single time step.
  !-------------------------------------------------------------------------------
  USE TYPEdef, ONLY: i4b
  USE POP_Types
  USE POP_Constants
  !USE BGCModule, ONLY: np

CONTAINS


  !*******************************************************************************
  SUBROUTINE ZeroPOP(POP)
    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER:: g,k,l,c, np

    PRINT*,"ZEROPOP"

    np = SIZE(pop%pop_grid)

    DO g=1,np
       POP%pop_grid(g)%freq = 0 ! patch weighting
       POP%pop_grid(g)%freq_old = 0 ! patch weighting
       POP%pop_grid(g)%freq_old2 = 0 ! patch weighting
       POP%pop_grid(g)%biomass = 0 ! landscape stem biomass (weighted mean over patches)
       POP%pop_grid(g)%density= 0  ! landscape tree density (weighted mean over patches)
       POP%pop_grid(g)%hmean= 0  ! landscape mean treen height (weighted mean over patches)
       POP%pop_grid(g)%hmax= 0   ! landscape max tree height
       POP%pop_grid(g)%cmass_stem_bin= 0  ! biomass by height bin
       POP%pop_grid(g)%densindiv_bin= 0  ! density by height bin
       POP%pop_grid(g)%height_bin= 0  ! mean height in each bin
       POP%pop_grid(g)%diameter_bin= 0  ! mean diameter in each bin
       POP%pop_grid(g)%bin_labels= ' '  ! text strings for bin bounds
       POP%pop_grid(g)%cmass_sum= 0  ! landscape biomass
       POP%pop_grid(g)%densindiv= 0  ! landscape density of individuals
       POP%pop_grid(g)%height_mean= 0
       POP%pop_grid(g)%height_max= 0
       POP%pop_grid(g)%basal_area= 0
       POP%pop_grid(g)%stress_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%fire_mortality = 0 ! (kg C m-2 y-1)
       POP%pop_grid(g)%growth= 0
       POP%pop_grid(g)%crown_cover = 0
       POP%pop_grid(g)%crown_area = 0
       POP%pop_grid(g)%crown_volume = 0

       DO k=1,NPATCH2D
          POP%pop_grid(g)%patch(k)%factor_recruit= 0
          POP%pop_grid(g)%patch(k)%biomass= 0  ! total biomass in patch
          POP%pop_grid(g)%patch(k)%biomass_old= 0
          POP%pop_grid(g)%patch(k)%biomass_old2 = 0
          POP%pop_grid(g)%patch(k)%stress_mortality = 0  ! biomass lost in each patch due to stress
          POP%pop_grid(g)%patch(k)%fire_mortality= 0 ! biomass lost in each patch due to fire partial dist
          POP%pop_grid(g)%patch(k)%crowding_mortality = 0
          POP%pop_grid(g)%patch(k)%mortality = 0
          POP%pop_grid(g)%patch(k)%growth= 0  ! biomass growth in each patch due stem increment
          POP%pop_grid(g)%patch(k)%disturbance_interval= 0   ! prescribed disturbance(s) interval for this patch
          POP%pop_grid(g)%patch(k)%first_disturbance_year = 0
          POP%pop_grid(g)%patch(k)%age= 0  ! number of years since last disturbance(s)
          POP%pop_grid(g)%patch(k)%age_old= 0  ! number of years since last disturbance(s)
          POP%pop_grid(g)%patch(k)%id = 0
          DO l=1,NLAYER
             POP%pop_grid(g)%patch(k)%Layer(L)%ncohort = 0 ! number of cohorts with density >0
             POP%pop_grid(g)%patch(k)%Layer(L)%biomass = 0 ! layer biomass
             POP%pop_grid(g)%patch(k)%Layer(L)%density= 0  ! layer tree density
             POP%pop_grid(g)%patch(k)%Layer(L)%hmean= 0  ! layer mean tree height (weighted mean over patches)
             POP%pop_grid(g)%patch(k)%Layer(L)%hmax= 0   ! layer max tree height
             DO c = 1,NCOHORT_MAX
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%id = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%age = 0 ! cohort age
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%biomass = 0 ! cohort biomass
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%density = 0 ! landscape tree density (weighted mean over patches)
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%frac_resource_uptake = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%height = 0
                POP%pop_grid(g)%patch(k)%Layer(L)%cohort(c)%diameter = 0
             ENDDO
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE ZeroPOP
  !*******************************************************************************

  SUBROUTINE InitPOP1D_Poisson(POP, mean_disturbance_interval)
    ! Initialises vector of patches with contrasting prescribed disturbance intervals
    ! Strategy: patch maximum ages should follow a Poisson distribution centred on the mean disturbance interval
    ! starting with 4 years maximum age and with four replicates for each patch maximum age.
    ! Starting year: 1a,2a,3a,4a where a is maximum age.
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) ::  mean_disturbance_interval(:,:)
    INTEGER(i4b) :: j, k, g, ipatch, idist, p, c, n, i
    INTEGER(i4b) :: disturbance_interval
    INTEGER(i4b):: patch_disturbance_interval1(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: patch_first_disturbance_year1(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: patch_disturbance_interval2(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: patch_first_disturbance_year2(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: Poisson_age(1000),Poisson_freq(1000)
    REAL(dp):: Poisson_weight(1000), CumPoisson_weight(1000)
    INTEGER(i4b):: disturbances_per_timebase, timebase
    INTEGER:: i_min, i_max, age_sample(2,NAGE_MAX), tmp(NAGE_MAX)
    INTEGER:: age_tmp, tmp_unique(NAGE_MAX), n_age, np
    REAL(dp):: disturbance_freq, tmp1 ! (BEN)

    np = SIZE(pop%pop_grid)

    DO g=1,np
       idist = 1
       disturbance_freq=1.0/REAL(mean_disturbance_interval(g,idist))
       DO p = 1,1000
          Poisson_age(p) = p
          Poisson_weight(p) = Exponential(disturbance_freq,p)
          CumPoisson_weight(p) = CumExponential(disturbance_freq,REAL(p,dp))
       ENDDO
       ! sample ages with equally spaced cumulative probabilities
       DO k =1,NAGE_MAX
          tmp1 = REAL(k)/REAL(NAGE_MAX)*0.95
          IF (tmp1.GT.CumPoisson_weight(1)) THEN
             i_max = MAXLOC(Poisson_age,1,CumPoisson_weight.LE.tmp1)
             tmp(k) = Poisson_age(i_max)
          ELSE
             tmp(k) = Poisson_age(1)
          ENDIF
       ENDDO

       ! subset to unique ages
       n=0
       age_tmp = -1
       DO i = 1, NAGE_MAX
          IF (tmp(i).NE.age_tmp) n = n+1
          tmp_unique(n) = tmp(i)
          age_tmp = tmp(i)
          n_age = n
       ENDDO
       IF (n_age.LT.NAGE_MAX) THEN
          DO i=n_age+1,NAGE_MAX
             IF ((i.EQ.n_age+1).AND.((tmp_unique(2)-tmp_unique(1)).GT.1)) THEN
                tmp_unique(i) =tmp_unique(1) +1
             ELSEIF ((i-2).GT.0) THEN
                tmp_unique(i) = tmp_unique(i-1) + tmp_unique(i-1)-tmp_unique(i-2)
             ELSE
                tmp_unique(i) = tmp_unique(i-1) + 1
             ENDIF
          ENDDO
          tmp = tmp_unique
       ENDIF

       age_sample(idist,:) = tmp

       k = 0
       DO j=1,NAGE_MAX
          disturbance_interval = age_sample(idist,j)
          DO c = 1,PATCH_REPS
             k = k+1

             POP%pop_grid(g)%patch(k)%disturbance_interval(1) = disturbance_interval
             POP%pop_grid(g)%patch(k)%first_disturbance_year(1) = (disturbance_interval)*c/PATCH_REPS
             POP%pop_grid(g)%patch(k)%age = 0
             POP%pop_grid(g)%patch(k)%id = k

          ENDDO
       ENDDO
       POP%pop_grid(g)%npatch_active = NPATCH

    ENDDO

  END SUBROUTINE InitPOP1D_Poisson
  !*******************************************************************************
  SUBROUTINE InitPOP2D_Poisson(POP, mean_disturbance_interval)
    ! Initialises vector of patches with contrasting prescribed disturbance intervals
    ! Strategy: patch maximum ages should follow a Poisson distribution centred on the mean disturbance interval
    ! starting with 4 years maximum age and with four replicates for each patch maximum age.
    ! Starting year: 1a,2a,3a,4a where a is maximum age.
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) ::  mean_disturbance_interval(:,:)
    INTEGER(i4b) :: j, k, g, ipatch, idist, p, c, n, i
    INTEGER(i4b) :: disturbance_interval
    INTEGER(i4b):: patch_disturbance_interval1(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: patch_first_disturbance_year1(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: patch_disturbance_interval2(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: patch_first_disturbance_year2(NPATCH1D,NPATCH1D)
    INTEGER(i4b):: Poisson_age(1000),Poisson_freq(1000)
    REAL(dp):: Poisson_weight(1000), CumPoisson_weight(1000)
    INTEGER(i4b):: disturbances_per_timebase, timebase
    INTEGER:: i_min, i_max, age_sample(2,NAGE_MAX), tmp(NAGE_MAX)
    INTEGER:: age_tmp, tmp_unique(NAGE_MAX), n_age, np
    REAL(dp):: disturbance_freq, tmp1
    INTEGER:: tmp2

    np = SIZE(POP%pop_grid)

    DO g=1,np

       ! calculate Poisson weights for each of the 2 mean disturbance intervals
       DO idist=1,2
          disturbance_freq=1.0/REAL(mean_disturbance_interval(g,idist))
          DO p = 1,1000
             Poisson_age(p) = p
             Poisson_weight(p) = Exponential(disturbance_freq,p)
             CumPoisson_weight(p) = CumExponential(disturbance_freq,REAL(p,dp))
          ENDDO
          ! sample ages with equally spaced cumulative probabilities
          DO k =1,NAGE_MAX
             tmp1 = REAL(k)/REAL(NAGE_MAX)*0.95
             IF (tmp1.GT.CumPoisson_weight(1)) THEN
                i_max = MAXLOC(Poisson_age,1,CumPoisson_weight.LE.tmp1)
                tmp(k) = Poisson_age(i_max)
             ELSE
                tmp(k) = Poisson_age(1)
             ENDIF
          ENDDO

          ! subset to unique ages
          n=0
          age_tmp = -1
          DO i = 1, NAGE_MAX
             IF (tmp(i).NE.age_tmp) n = n+1
             tmp_unique(n) = tmp(i)
             age_tmp = tmp(i)
             n_age = n
          ENDDO
          IF (n_age.LT.NAGE_MAX) THEN
             DO i=n_age+1,NAGE_MAX
                IF ((i.EQ.n_age+1).AND.((tmp_unique(2)-tmp_unique(1)).GT.1)) THEN
                   tmp_unique(i) =tmp_unique(1) +1
                ELSEIF ((i-2).GT.0) THEN
                   tmp_unique(i) = tmp_unique(i-1) + tmp_unique(i-1)-tmp_unique(i-2)
                ELSE
                   tmp_unique(i) = tmp_unique(i-1) + 1
                ENDIF
             ENDDO
             tmp = tmp_unique
          ENDIF

          age_sample(idist,:) = tmp

          ! get first disturbance year
          k = 0
          DO j=1,NAGE_MAX
             disturbance_interval = age_sample(idist,j)
             DO c = 1,PATCH_REPS
                k = k+1
                IF (idist.EQ.1) THEN
                   patch_disturbance_interval1(k,k)   = disturbance_interval
                   tmp2 = disturbance_interval*c/PATCH_REPS
                   patch_first_disturbance_year1(k,k) = tmp2
                   !write(*,*) tmp2
                ELSEIF (idist.EQ.2) THEN
                   patch_disturbance_interval2(k,k)   = disturbance_interval
                   patch_first_disturbance_year2(k,k) = disturbance_interval*c/PATCH_REPS
                ENDIF
             ENDDO

          ENDDO

       ENDDO


       DO k=1,NPATCH1D
          patch_disturbance_interval1(k,:)          = patch_disturbance_interval1(k,k)
          patch_first_disturbance_year1(k,:)        = patch_first_disturbance_year1(k,k)
          patch_disturbance_interval2(:,k)          = patch_disturbance_interval2(k,k)
          patch_first_disturbance_year2(:,k)        = patch_first_disturbance_year2(k,k)
          POP%pop_grid(g)%patch(k)%Layer(:)%biomass = 0.
       ENDDO

       ipatch=1
       DO k=1,NPATCH1D
          DO j=1,NPATCH1D
             POP%pop_grid(g)%patch(ipatch)%disturbance_interval(1)   = patch_disturbance_interval1(k,j)
             POP%pop_grid(g)%patch(ipatch)%disturbance_interval(2)   = patch_disturbance_interval2(k,j)
             POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(1) = patch_first_disturbance_year1(k,j)
             POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(2) = patch_first_disturbance_year2(k,j)
             POP%pop_grid(g)%patch(ipatch)%age = 0
             POP%pop_grid(g)%patch(ipatch)%id = ipatch
             ipatch = ipatch + 1
          ENDDO
       ENDDO

       IF (NPATCH2D.GT.NPATCH+1) THEN
          DO k=1,NPATCH1D
             POP%pop_grid(g)%patch(ipatch)%disturbance_interval(1)   = 0
             POP%pop_grid(g)%patch(ipatch)%disturbance_interval(2)   = patch_disturbance_interval2(k,k)
             POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(1) = 0
             POP%pop_grid(g)%patch(ipatch)%first_disturbance_year(2) = patch_first_disturbance_year2(k,k)
             POP%pop_grid(g)%patch(ipatch)%age = 0
             POP%pop_grid(g)%patch(ipatch)%id  = ipatch
             ipatch = ipatch + 1
          ENDDO
       ENDIF

       POP%pop_grid(g)%npatch_active = NPATCH

    ENDDO

  END SUBROUTINE InitPOP2D_Poisson
  !*******************************************************************************

  SUBROUTINE POPStep(POP, StemNPP, disturbance_interval, disturbance_intensity,LAI,frac_intensity1,precip )
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    REAL(dp), INTENT(IN) :: StemNPP(:,:)
    REAL(dp), INTENT(IN) :: disturbance_intensity(:,:)
    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:)
    REAL(dp), INTENT(IN) ::  LAI(:)
    REAL(dp), INTENT(IN), OPTIONAL :: frac_intensity1(:), precip(:)
    INTEGER(i4b) :: idisturb, it

    pop%it_pop = pop%it_pop + 1

    it = pop%it_pop

    IF (PRESENT(precip)) THEN
       CALL PatchAnnualDynamics(POP, StemNPP,disturbance_interval, it, precip)
    ELSE
       CALL PatchAnnualDynamics(POP, StemNPP,disturbance_interval, it)
    ENDIF



    IF (NDISTURB.EQ.1) THEN
       IF (PRESENT(precip)) THEN
          CALL Patch_disturb(POP,it,1,precip)
       ELSE
          CALL Patch_disturb(POP,it,1)
       ENDIF
    ELSEIF (NDISTURB.EQ.2) THEN
       IF (PRESENT(frac_intensity1)) THEN
          IF (PRESENT(precip)) THEN
             CALL Patch_partial_disturb(POP,it,1,disturbance_intensity,precip,frac_intensity1)
          ELSE
             CALL Patch_partial_disturb(POP,it,1,disturbance_intensity,frac_intensity1=frac_intensity1)
          ENDIF
       ELSE
          IF (PRESENT(precip)) THEN
             CALL Patch_partial_disturb(POP,it,1,disturbance_intensity,precip=precip)
          ELSE
             CALL Patch_partial_disturb(POP,it,1,disturbance_intensity)
          ENDIF
       ENDIF
       IF (PRESENT(precip)) THEN
          CALL Patch_disturb(POP,it,2,precip)
       ELSE
          CALL Patch_disturb(POP,it,2)
       ENDIF
    ENDIF


    DO idisturb = 1,NDISTURB
       CALL GetUniqueAgeFrequencies(POP, disturbance_interval, idisturb, it)
    ENDDO


    CALL GetPatchFrequencies(POP)

    IF (PRESENT(precip)) THEN
       CALL GetDiagnostics(pop, LAI, it,precip)
    ELSE
       CALL GetDiagnostics(pop, LAI, it)
    ENDIF




  END SUBROUTINE POPStep
  !*******************************************************************************
  SUBROUTINE PatchAnnualDynamics(pop, StemNPP, disturbance_interval, it, precip)
    IMPLICIT NONE

    TYPE( POP_TYPE ), INTENT(INOUT) :: pop
    REAL(dp), INTENT(IN)            :: StemNPP(:,:)
    INTEGER(i4b), INTENT(IN)        ::  disturbance_interval(:,:)
    REAL(dp), INTENT(IN), OPTIONAL  :: precip(:)
    INTEGER(i4b), INTENT(IN)        :: it

    REAL(dp) :: f, mu, densindiv, cmass
    REAL(dp) :: tmp, cmass_stem_sum,cmass_stem_inc
    INTEGER(i4b) :: j, k,c, ncohort
    INTEGER(i4b) :: ivec(NCOHORT_MAX), nc, tmp1(NPATCH2D), tmp2(NPATCH2D), np, idisturb
    REAL(dp) :: growth_efficiency,cmass_stem
    REAL(dp) :: mort, mort_bg, fire_mort
    REAL(dp) :: s2, cpc, crown_area
    REAL(dp) :: mort_cpc
    !Maximum disturbance interval that makes sense is 1.75*2^((NPATCH)/4) as we don't always have a
    !patch that has been undisturbed longer than this
    !This is also the minimum spinup time

    idisturb = 1
    np = SIZE(POP%POP_grid)

    ! growth
    ! Distributes layer biomass increment among cohorts and increments age
    ! calculate fractional resource uptake by each cohort
    DO j=1,np
       pop%pop_grid(j)%freq_old2 = pop%pop_grid(j)%freq_old
       pop%pop_grid(j)%freq_old = pop%pop_grid(j)%freq

       DO k=1,NPATCH2D
          pop%pop_grid(j)%patch(k)%biomass_old2 = pop%pop_grid(j)%patch(k)%biomass_old
          pop%pop_grid(j)%patch(k)%biomass_old = pop%pop_grid(j)%patch(k)%biomass
          pop%pop_grid(j)%patch(k)%age_old = pop%pop_grid(j)%patch(k)%age
          pop%pop_grid(j)%patch(k)%growth = 0.0
          tmp = 0.0
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             tmp = tmp + (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density
          ENDDO
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake = &
                  (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density/tmp

             ! increment biomass in cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass =  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass +  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake*StemNPP(j,1)
             pop%pop_grid(j)%patch(k)%growth = pop%pop_grid(j)%patch(k)%growth +  &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake*StemNPP(j,1)
             ! increment cohort age
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%age = &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%age + 1

          ENDDO
          ! Layer biomass (summed over cohorts)
          nc = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          pop%pop_grid(j)%patch(k)%Layer(1)%biomass = SUM(pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)%biomass)
       ENDDO
    ENDDO

    ! Mortality
    !Implements resource stress mortality and crowding mortality for all cohorts in layer
    ! first, calculate fractional resource uptake by each cohort
    DO j=1,np
       pop%pop_grid(j)%patch(:)%fire_mortality = 0.0
       DO k=1,NPATCH2D

          tmp = 0.0
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             tmp = tmp + (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density
          ENDDO
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake = &
                  (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass/ &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density)**POWERbiomass * &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density/tmp

          ENDDO
       ENDDO
    ENDDO
    DO j=1,np
       DO k=1,NPATCH2D
          nc = 0
          ivec = 0
          pop%pop_grid(j)%patch(k)%stress_mortality = 0.0
          pop%pop_grid(j)%patch(k)%fire_mortality = 0.0
          pop%pop_grid(j)%patch(k)%crowding_mortality = 0.0
          pop%pop_grid(j)%patch(k)%mortality = 0.0
          crown_area = 0.0
          DO c=1,pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             cmass_stem = pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
             cmass_stem_inc=StemNPP(j,1)*pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%frac_resource_uptake

             growth_efficiency=cmass_stem_inc/(cmass_stem**(POWERGrowthEfficiency))
             mort=MORT_MAX/(1.0+(growth_efficiency/GROWTH_EFFICIENCY_MIN)**Pmort)

             pop%pop_grid(j)%patch(k)%stress_mortality = pop%pop_grid(j)%patch(k)%stress_mortality + mort*cmass_stem

             if (ALLOM_SWITCH.eq.1) then
                crown_area = crown_area + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density* &
                     PI*(pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100.*0.14)**2 !! assumes crown radius (m) = 0.14 * dbh (cm)
             else
                crown_area = crown_area + pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density* &
                     k_allom1 * pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter ** k_rp  !! assumes crown radius (m) = 0.07 * dbh (cm)
             endif

             cpc = 1. - exp(-crown_area)
             if (cpc .gt. 1.e-8) then
                mort_cpc = exp(alpha_cpc * (1. - 1./cpc))
             else
                mort_cpc = 0.
             endif
             mort = mort + (mort_cpc*CrowdingFactor)*(1.-mort)
             pop%pop_grid(j)%patch(k)%crowding_mortality = pop%pop_grid(j)%patch(k)%crowding_mortality + &
                  (mort_cpc*CrowdingFactor)*cmass_stem
             mort = min(mort,cmass_stem_inc/cmass_stem)

             pop%pop_grid(j)%patch(k)%mortality = pop%pop_grid(j)%patch(k)%mortality + mort*cmass_stem

             IF (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0).AND. &
                  (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))).OR. &
                  (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))) THEN

                fire_mort = 0.0
             ELSE
                fire_mort = 0.0
             ENDIF
             pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality +fire_mort*cmass_stem

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = cmass_stem*(1.-mort-fire_mort)

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = &
                  pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density*(1.-mort-fire_mort)
             IF (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.LT.DENSINDIV_MIN) THEN
                ! remove cohort
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                pop%pop_grid(j)%patch(k)%stress_mortality = pop%pop_grid(j)%patch(k)%stress_mortality + &
                     pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0
             ELSE
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                nc = nc+1
                ivec(nc)=c
             ENDIF
          ENDDO
          ! SHUFFLE if necessary to remove zero-density cohorts
          IF (nc.LT.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) THEN
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id      = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0
          ENDIF

       ENDDO
    ENDDO

    ! recruitment
    IF (PRESENT(precip)) THEN
       CALL layer_recruitment(pop, precip)
    ELSE
       CALL layer_recruitment(pop)
    ENDIF

    ! Update time since last patch disturbance
    DO j=1,np
       DO k=1,NPATCH

          if(j.eq.1.and.k.eq.42) write(*,*)  pop%pop_grid(j)%patch(k)%age(1)

          IF (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0).AND. &
               (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))).OR. &
               (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))) THEN

             ! pop%pop_grid(j)%patch(k)%age(idisturb)=0
             pop%pop_grid(j)%patch(k)%age(1) = pop%pop_grid(j)%patch(k)%age(1) + 1
          ELSE
             pop%pop_grid(j)%patch(k)%age(1) = pop%pop_grid(j)%patch(k)%age(1) + 1
          ENDIF

          if(j.eq.1.and.k.eq.42) write(*,*)  pop%pop_grid(j)%patch(k)%age(1)

          IF (NDISTURB.EQ.2) THEN
             pop%pop_grid(j)%patch(k)%age(2) = pop%pop_grid(j)%patch(k)%age(2) + 1
          ENDIF

       ENDDO
       IF (NPATCH2D.GT.NPATCH+1) THEN
          DO k= NPATCH+1,NPATCH2D
             pop%pop_grid(j)%patch(k)%age(2) = pop%pop_grid(j)%patch(k)%age(2) + 1
          ENDDO
       ENDIF
    ENDDO

    !tmp1 = pop%pop_grid(1)%patch(:)%age(1)
    !tmp2 = pop%pop_grid(1)%patch(:)%age(2)

  END SUBROUTINE PatchAnnualDynamics
  !*******************************************************************************
  SUBROUTINE GetUniqueAgeFrequencies(pop, disturbance_interval, idisturb, it)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:), idisturb, it
    INTEGER(i4b) :: g, i,j,k,ct,lastct,agecopy,idcopy
    REAL(dp), ALLOCATABLE :: midpoint(:)
    INTEGER(i4b), ALLOCATABLE :: ranked_age(:), ranked_age_init(:)
    !REAL(dp), ALLOCATABLE, SAVE ::freq_ranked_age_unique(:,:), ranked_age_unique(:,:)
    INTEGER(i4b) ::  tmp_count, tmp_i, age_tmp
    INTEGER(i4b), ALLOCATABLE :: ranked_age_unique_id(:), ranked_age_id(:), counter(:)
    REAL(dp), ALLOCATABLE :: tmp(:), freq_tmp(:), freq_tmp1(:)
    REAL(dp) :: p,cump,lastcump, freq, tmp1
    INTEGER(i4b) :: n_age ! number of unique ages
    INTEGER(i4b) :: npatch_active ! number of active patches
    !REAL(dp) :: timebase ! (BEN)
    REAL(dp):: disturbance_freq ! (BEN)
    INTEGER(i4b) :: i_max, age_max, Poisson_age(1000), np
    REAL(dp):: Poisson_weight(1000), CumPoisson_weight(1000)
    INTEGER(i4b), ALLOCATABLE :: bound(:,:), unique_age(:)

    !Fills array freq with weights (frequencies across landscape) for each unique age
    ! given specified mean disturbance interval

    np = SIZE(POP%POP_grid)
    DO g=1,np

       npatch_active = NPATCH2D
       IF (.NOT.ALLOCATED(midpoint)) ALLOCATE(midpoint(npatch_active))
       IF (.NOT.ALLOCATED(counter)) ALLOCATE(counter(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age)) ALLOCATE(ranked_age(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age_init)) ALLOCATE(ranked_age_init(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age_id)) ALLOCATE(ranked_age_id(npatch_active))
       IF (.NOT.ALLOCATED(ranked_age_unique_id)) ALLOCATE(ranked_age_unique_id(npatch_active))
       IF (.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(npatch_active))
       IF (.NOT.ALLOCATED(freq_tmp)) ALLOCATE(freq_tmp(npatch_active))
       IF (.NOT.ALLOCATED(freq_tmp1)) ALLOCATE(freq_tmp1(npatch_active))


       ! rank patches in order of age
       pop%pop_grid(g)%ranked_age_unique(:, idisturb) = 0
       ranked_age_init = pop%pop_grid(g)%patch%age(idisturb)
       ranked_age = pop%pop_grid(g)%patch%age(idisturb)
       ranked_age_id = pop%pop_grid(g)%patch%id
       ranked_age_unique_id = 0
       freq_tmp = 0.0
       freq = 0
       pop%pop_grid(g)%freq_ranked_age_unique(:, idisturb) = 0.0
       midpoint = 0.0


       DO i = 1, npatch_active -1
          DO j = i+1, npatch_active
             IF (ranked_age(i).GT.ranked_age(j)) THEN
                agecopy          = ranked_age(i)
                idcopy           = ranked_age_id(i)
                ranked_age(i)    = ranked_age(j)
                ranked_age_id(i) = ranked_age_id(j)
                ranked_age(j)    = agecopy
                ranked_age_id(j) = idcopy
             ENDIF
          ENDDO
       ENDDO

       ! subset to unique ages
       k=0
       age_tmp = -1
       DO i = 1, npatch_active
          IF (ranked_age(i).NE.age_tmp) k = k+1
          pop%pop_grid(g)%ranked_age_unique(k, idisturb) = ranked_age(i)
          ranked_age_unique_id(k) = ranked_age_id(i)
          age_tmp = ranked_age(i)
          n_age  = k
       ENDDO

       disturbance_freq=1.0/REAL(disturbance_interval(g,idisturb))
       DO i =1,1000
          Poisson_age(i) = i
          CumPoisson_weight(i) = CumExponential(disturbance_freq,REAL(i,dp))
       ENDDO


       ! construct upper and lower bounds for each unique age: these set the range of ages to be
       ! represented by an unique age
       ALLOCATE(bound(n_age,2))
       ALLOCATE (unique_age(n_age))
       bound = 0
       unique_age = pop%pop_grid(g)%ranked_age_unique(1:n_age,idisturb)
       DO i=1,n_age
          IF (unique_age(i).EQ.0) THEN
             bound(i,1) = 0
             bound(i,2) = 0
          ELSEIF ((i.EQ.1).AND.(unique_age(i).GT.0)) THEN
             bound(i,1) = 0
             bound(i,2) = unique_age(i)
          ELSEIF ((unique_age(i).GT.0).AND.(i.GT.1).AND.(unique_age(i-1).EQ.unique_age(i)-1)) THEN
             bound(i,1) = unique_age(i)
             bound(i,2) = unique_age(i)
          ELSEIF ((unique_age(i).GT.0).AND.(i.GT.1).AND.(unique_age(i-1).NE.unique_age(i)-1)) THEN
             bound(i,1) = bound(i-1,2)+1
             IF (i.LT.n_age) THEN
                bound(i,2) = (unique_age(i)+ unique_age(i+1))/2
             ELSE
                i_max = MAXLOC(Poisson_age, 1, CumPoisson_weight.GE.0.99)
                !bound(i, 2) = Poisson_age(i_max)
                bound(i, 2) = bound(i,1)
             ENDIF
          ENDIF

       ENDDO

       ! calculate weighting for each unique age
       DO i=1,n_age
          DO j = bound(i,1),bound(i,2)
             freq_tmp(i) = freq_tmp(i) + REALExponential(disturbance_freq,REAL(j,dp))
             !freq_tmp(i) = freq_tmp(i) + Exponential(disturbance_freq,j)
          ENDDO
       ENDDO

       pop%pop_grid(g)%freq_ranked_age_unique(1:npatch_active,idisturb) = freq_tmp
       pop%pop_grid(g)%n_age(idisturb) = n_age

       DEALLOCATE (bound)
       DEALLOCATE (unique_age)

    ENDDO

  END SUBROUTINE GetUniqueAgeFrequencies

  !*******************************************************************************
  SUBROUTINE GetPatchFrequencies(pop)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b) :: n1, n2, g, REPCOUNT, tmp1(NPATCH1D), np
    REAL(dp) ::  tmp2(NPATCH1D), tmp3(NPATCH1D)

    np = SIZE(Pop%pop_grid)
    DO g=1,np

       pop%pop_grid(g)%freq = 0.0
       IF (NDISTURB.EQ.1) THEN
          DO n1=1,pop%pop_grid(g)%n_age(1)
             repcount = COUNT(pop%pop_grid(g)%patch(:)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1))
             WHERE (pop%pop_grid(g)%patch(:)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1))
                pop%pop_grid(g)%freq = pop%pop_grid(g)%freq_ranked_age_unique(n1,1) /REAL(repcount,dp)
             ENDWHERE
          ENDDO


       ELSEIF (NDISTURB.EQ.2) THEN
          ! first calculate weights for patches with age(2)>age(1)
          DO n1=1,pop%pop_grid(g)%n_age(1)
             DO n2=1,pop%pop_grid(g)%n_age(2)
                repcount = COUNT((pop%pop_grid(g)%patch(1:NPATCH)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1)).AND. &
                     (pop%pop_grid(g)%patch(1:NPATCH)%age(2).EQ.pop%pop_grid(g)%ranked_age_unique(n2,2)))
                WHERE ((pop%pop_grid(g)%patch(1:NPATCH)%age(1).EQ.pop%pop_grid(g)%ranked_age_unique(n1,1)).AND. &
                     (pop%pop_grid(g)%patch(1:NPATCH)%age(2).EQ.pop%pop_grid(g)%ranked_age_unique(n2,2)))
                   pop%pop_grid(g)%freq(1:NPATCH) = &
                        pop%pop_grid(g)%freq_ranked_age_unique(n1,1)*pop%pop_grid(g)%freq_ranked_age_unique(n2,2) &
                        /REAL(repcount,dp)
                ENDWHERE
             ENDDO
          ENDDO


          IF (NPATCH2D.GT.NPATCH) THEN
             ! second calculate weights for patches with age(1)>age(2)
             !tmp1 =  pop%pop_grid(g)%patch(NPATCH+1:NPATCH2D)%age(2)
             tmp2 = 0.
             tmp3 = 0.
             DO n1=1,pop%pop_grid(g)%n_age(1)
                DO n2=1,pop%pop_grid(g)%n_age(2)
                   IF (n1>n2) THEN
                      repcount = COUNT(tmp1.EQ.pop%pop_grid(g)%ranked_age_unique(n2,2))
                      tmp3 = 0.
                      WHERE ((tmp1.EQ.pop%pop_grid(g)%ranked_age_unique(n2,2)))

                         tmp3 =  pop%pop_grid(g)%freq_ranked_age_unique(n1,1)*pop%pop_grid(g)%freq_ranked_age_unique(n2,2)  &
                              /REAL(repcount,dp)
                      ENDWHERE
                      tmp2 = tmp2+tmp3
                   ENDIF
                ENDDO
             ENDDO
             !pop%pop_grid(g)%freq(NPATCH+1:NPATCH2D) = tmp2
          ENDIF

       ENDIF

       pop%pop_grid(g)%freq = pop%pop_grid(g)%freq/SUM(pop%pop_grid(g)%freq)
    ENDDO

  END SUBROUTINE GetPatchFrequencies

  !*******************************************************************************
  SUBROUTINE GetDiagnostics(pop,LAI, it, precip)
    ! Gets diagnostic data for current landscape structure
    IMPLICIT NONE
    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    REAL(dp), INTENT(IN) ::  LAI(:)
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    INTEGER(i4b), INTENT(IN) :: it
    INTEGER(i4b) :: P, g,i,j,ct, ct_highres
    REAL(dp) :: limits(HEIGHT_BINS+1)
    REAL(dp) :: ht, htmax, cmass_stem,densindiv, freq, freq_old
    CHARACTER(len=12) :: string1, string2
    CHARACTER(len=9) :: fmt
    INTEGER(i4b) :: npatch_active  ! number of active patches
    INTEGER(i4b) :: np, i_height
    REAL(dp) :: diam,basal, cump
    REAL(dp) :: patch_crown_area(NPATCH2D), patch_crown_cover(NPATCH2D)
    REAL(dp), ALLOCATABLE :: height_list(:), height_list_weight(:)
    REAL(dp) :: height_copy, weight_copy, Pwc, FAVD
    INTEGER(i4b), PARAMETER :: HEIGHT_BINS_highres=100 ! bins for assessing height_max
    REAL(dp), ALLOCATABLE :: limits_highres(:), DENSINDIV_HIGHRES(:)


    fmt = '(f5.1)'
    limits(1) = 0.
    IF(.NOT.ALLOCATED(limits_highres)) ALLOCATE(limits_highres(HEIGHT_BINS_highres+1))
    IF(.NOT.ALLOCATED(DENSINDIV_HIGHRES)) ALLOCATE(DENSINDIV_HIGHRES(HEIGHT_BINS_highres))
    limits_highres(1) = 0.
    np = SIZE(Pop%pop_grid)

    DO g=1,np
       npatch_active = NPATCH2D
       IF (MAX_HEIGHT_SWITCH.EQ.1) THEN
          ALLOCATE(height_list(NPATCH2D*NCOHORT_MAX))
          ALLOCATE(height_list_weight(NPATCH2D*NCOHORT_MAX))
       ENDIF


       DO i=1,HEIGHT_BINS
          limits(i+1) = BIN_POWER**REAL(i)
          WRITE(string1,fmt) (limits(i))
          WRITE(string2,fmt) (limits(i+1))
          pop%pop_grid(g)%bin_labels(i) = 'Height_'//TRIM(ADJUSTL(string1))//'-'//TRIM(ADJUSTL(string2))//'m'
          pop%pop_grid(g)%cmass_stem_bin(i) = 0.0
          pop%pop_grid(g)%densindiv_bin(i) = 0.0
          pop%pop_grid(g)%cmass_stem_bin(i) = 0.0
          pop%pop_grid(g)%height_bin(i) = REAL(limits(i)+limits(i+1))/2.
          pop%pop_grid(g)%diameter_bin(i) = ((REAL(limits(i))/Kbiometric)**(3/2)+(REAL(limits(i+1))/Kbiometric)**(3/2))/2.
       ENDDO

       DO i=1,HEIGHT_BINS_highres
          limits_highres(i+1) = REAL(i)
       ENDDO

       height_list = 0.0
       height_list_weight = 0.0
       i_height = 0
       pop%pop_grid(g)%cmass_sum = 0.0
       pop%pop_grid(g)%height_mean = 0.0
       pop%pop_grid(g)%fire_mortality = 0.0
       pop%pop_grid(g)%stress_mortality = 0.0
       pop%pop_grid(g)%growth = 0.0
       pop%pop_grid(g)%basal_area = 0.0
       pop%pop_grid(g)%densindiv = 0.0
       pop%pop_grid(g)%height_max = 0.0
       pop%pop_grid(g)%crown_cover = 0.0
       pop%pop_grid(g)%crown_area = 0.0
       pop%pop_grid(g)%crown_volume = 0.0
       densindiv_highres = 0.0
       ! loop through patches
       DO P = 1, npatch_active
          pop%pop_grid(g)%patch(p)%biomass = 0.0
          pop%pop_grid(g)%patch(p)%layer(1)%biomass = 0.0
          pop%pop_grid(g)%patch(p)%layer(1)%density = 0.0
          patch_crown_area(p) = 0.0
          patch_crown_cover(p) = 0.0

          ! loop through cohorts
          DO i = 1, pop%pop_grid(g)%patch(p)%layer(1)%ncohort
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density
             !ht=(Kbiometric**(3.0/5.0))*((2.0/PI)**(2.0/5.0))*((cmass_stem/densindiv/WD)**(2.0/5.0))
             !

             IF (ALLOM_SWITCH.EQ.1.AND.PRESENT(precip)) THEN
                ht=GetHeight(precip(g),cmass_stem,densindiv)
                CALL Allometry(ht,cmass_stem,densindiv,diam,basal)
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%height = ht
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter = diam

             ELSEIF (ALLOM_SWITCH.EQ.0) THEN

                ht=(Kbiometric**(3.0/4.0))*(4.*cmass_stem/(densindiv*WD*PI))**(1.0/4.0)
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%height = ht
                diam = (ht/Kbiometric)**(3/2)
                pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter = diam
                basal=PI*(pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter/2.0)* &
                     (pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%diameter/2.0)*densindiv*1.e4
             ENDIF

             ! get bin
             DO j=1,HEIGHT_BINS
                IF (ht.GT.limits(j)) ct = j
             ENDDO ! bins

             ! get high res bin
             DO j=1,HEIGHT_BINS_highres
                IF (ht.GT.limits_highres(j)) ct_highres = j
             ENDDO ! bins

             pop%pop_grid(g)%patch(p)%layer(1)%biomass = pop%pop_grid(g)%patch(p)%layer(1)%biomass + cmass_stem
             pop%pop_grid(g)%patch(p)%layer(1)%density = pop%pop_grid(g)%patch(p)%layer(1)%density + densindiv

             freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
             freq_old = pop%pop_grid(g)%freq_old(pop%pop_grid(g)%patch(p)%id)

             IF (diam*100..GT.5.) THEN
                patch_crown_area(p) = patch_crown_area(p) + densindiv*PI*(diam*100.*0.1492)**2 ! uses GC relationship
                ! pop%pop_grid(g)%crown_area = pop%pop_grid(g)%crown_area + freq*densindiv*PI*(diam*100.*0.1492)**2
                pop%pop_grid(g)%crown_volume = pop%pop_grid(g)%crown_volume + &
                     freq*densindiv*(4./3.)*PI*(diam*100.*0.1492)**2*(1.5*(diam*100.*0.1492))
                ! assumes vertical radius = 1.5 * horizontal radius
             ENDIF

             ! crown radius (m) = 0.1492 *DBH (cm)


             pop%pop_grid(g)%patch(p)%biomass = pop%pop_grid(g)%patch(p)%biomass + cmass_stem
             pop%pop_grid(g)%cmass_stem_bin(ct) = pop%pop_grid(g)%cmass_stem_bin(ct) + freq*cmass_stem
             pop%pop_grid(g)%densindiv_bin(ct) = pop%pop_grid(g)%densindiv_bin(ct) + freq*densindiv
             pop%pop_grid(g)%cmass_sum = pop%pop_grid(g)%cmass_sum + freq*cmass_stem
             pop%pop_grid(g)%densindiv = pop%pop_grid(g)%densindiv + freq*densindiv
             pop%pop_grid(g)%height_mean = pop%pop_grid(g)%height_mean + ht*freq*densindiv
             pop%pop_grid(g)%basal_area = pop%pop_grid(g)%basal_area +basal*freq
             densindiv_highres(ct_highres) = densindiv_highres(ct_highres) + freq*densindiv
             IF (MAX_HEIGHT_SWITCH.EQ.1) THEN
                i_height = i_height+1
                height_list(i_height) = ht
                height_list_weight(i_height) = densindiv*freq
             ENDIF
          ENDDO ! cohorts

          patch_crown_cover(p) = 1. -EXP(-patch_crown_area(p))
          !
          !pop%pop_grid(g)%crown_cover = pop%pop_grid(g)%crown_cover  + freq*patch_crown_cover(p)
          pop%pop_grid(g)%stress_mortality = pop%pop_grid(g)%stress_mortality + freq*pop%pop_grid(g)%patch(p)%stress_mortality
          pop%pop_grid(g)%fire_mortality = pop%pop_grid(g)%fire_mortality + freq*pop%pop_grid(g)%patch(p)%fire_mortality + &
               pop%pop_grid(g)%patch(p)%biomass_old*(freq_old-freq)
          !WRITE(333,"(6i5, 100e16.6)") it,  p, pop%pop_grid(g)%patch(p)%disturbance_interval,pop%pop_grid(g)%patch(p)%first_disturbance_year, &
          !     pop%pop_grid(g)%patch(p)%age_old(1), &
          !     pop%pop_grid(g)%patch(p)%age(1), &
          !     freq_old, freq, pop%pop_grid(g)%patch(p)%biomass_old, pop%pop_grid(g)%patch(p)%biomass



          pop%pop_grid(g)%growth =  pop%pop_grid(g)%growth + freq*pop%pop_grid(g)%patch(p)%growth

       ENDDO ! patches

       FAVD = LAI(g)/pop%pop_grid(g)%crown_volume ! foliage area volume density
       DO P = 1, npatch_active
          freq =  pop%pop_grid(g)%freq(pop%pop_grid(g)%patch(p)%id)
          ! loop through cohorts
          DO i = 1, pop%pop_grid(g)%patch(p)%layer(1)%ncohort
             cmass_stem = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%biomass
             densindiv = pop%pop_grid(g)%patch(p)%layer(1)%cohort(i)%density
             IF (ALLOM_SWITCH.EQ.1.AND.PRESENT(precip)) THEN
                ht=GetHeight(precip(g),cmass_stem,densindiv)
                CALL Allometry(ht,cmass_stem,densindiv,diam,basal)
             ELSEIF (ALLOM_SWITCH.EQ.0) THEN
                ht=(Kbiometric**(3.0/4.0))*(4.*cmass_stem/(densindiv*WD*PI))**(1.0/4.0)
                diam = (ht/Kbiometric)**(3/2)
             ENDIF

             IF (diam*100..GT.5.) THEN
                Pwc  = EXP(-0.5*(4./3.)*(diam*100.*0.1492)*1.5*FAVD) ! Pwc = exp(-G*FAVD*s)
                pop%pop_grid(g)%crown_area = pop%pop_grid(g)%crown_area + freq*densindiv*PI*(diam*100.*0.1492)**2*(1.-Pwc)
             ENDIF

          ENDDO ! cohorts


       ENDDO ! patches
       pop%pop_grid(g)%crown_cover = 1.-EXP(-pop%pop_grid(g)%crown_area)


       pop%pop_grid(g)%height_mean = pop%pop_grid(g)%height_mean/pop%pop_grid(g)%densindiv

       IF (MAX_HEIGHT_SWITCH.EQ.0) THEN
          ! Set landscape maximum height to centre of bin with <5% of trees in a bin of higher size classes
          cump = 0.
          j = 1
          DO WHILE (cump.LT.0.99)
             cump = cump + pop%pop_grid(g)%densindiv_bin(j)/pop%pop_grid(g)%densindiv
             pop%pop_grid(g)%height_max = pop%pop_grid(g)%height_bin(j)
             j = j+1
          ENDDO
       ELSEIF (MAX_HEIGHT_SWITCH.EQ.1) THEN

          ! sort height list
          DO i = 1, i_height -1
             DO j = i+1, i_height
                IF (height_list(i).GT.height_list(j)) THEN
                   height_copy = height_list(i)
                   weight_copy =  height_list_weight(i)
                   height_list(i) = height_list(j)
                   height_list_weight(i) = height_list_weight(j)
                   height_list(j) = height_copy
                   height_list_weight(j) = weight_copy
                ENDIF
             ENDDO
          ENDDO ! end sort height list

          ! normailse height list weights
          height_list_weight=height_list_weight/SUM(height_list_weight(1:i_height))
          cump = 0.
          j = 1
          DO WHILE (cump.LT.0.99)
             cump = cump + height_list_weight(j)
             pop%pop_grid(g)%height_max = height_list(j)
             j = j+1
          ENDDO
          DEALLOCATE(height_list)
          DEALLOCATE(height_list_weight)

       ELSEIF (MAX_HEIGHT_SWITCH.EQ.2) THEN
          cump = 0.
          j = 1
          densindiv_highres= densindiv_highres/SUM(densindiv_highres)
          DO WHILE ((cump.LT.0.95).AND.(j.LE.HEIGHT_BINS_highres))
             cump = cump + densindiv_highres(j)
             pop%pop_grid(g)%height_max = (limits_highres(j+1) + limits_highres(j))/2.
             j = j+1
          ENDDO
       ENDIF


    ENDDO ! end loop over grid cells

  END SUBROUTINE GetDiagnostics
  !*******************************************************************************
  SUBROUTINE Patch_partial_disturb(pop,it,idisturb,intensity,precip,frac_intensity1)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    INTEGER(i4b), INTENT(IN) :: it, idisturb
    REAL(dp), INTENT(IN) :: intensity(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:), frac_intensity1(:)
    INTEGER(i4b) :: j, k, i, g, c, nc, np
    INTEGER(i4b) ::  ivec(NCOHORT_MAX)
    REAL(dp) :: ht, diam
    REAL(dp) :: Psurvival_l, Psurvival_s, Psurvival, char_height

    np = SIZE(Pop%pop_grid)

    ! Kills a fraction of biomass in patch when prescribed disturbance interval is reached

    DO j=1,np
       DO k=1,NPATCH
          pop%pop_grid(j)%patch(k)%fire_mortality = 0.0

          IF (((pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0).AND. &
               (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))).OR. &
               (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb))) THEN


             ! loop through cohorts
             ivec = 0
             nc = 0
             DO c = 1, pop%pop_grid(j)%patch(k)%layer(1)%ncohort
                ! kill fraction of each cohort
                char_height = 3.7*(1.-EXP(-0.19*Intensity(j,1)))
                ht = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%height
                diam = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%diameter*100. ! diameter in cm
                IF ((ht.GT.8.5).AND.(ht.GT.char_height)) THEN
                   Psurvival_s  =(-0.0011*Intensity(j,1) -0.00002)*ht &
                        +(0.0075*Intensity(j,1)+1.)
                   !Psurvival_s = 1+0.029*Intensity(j,1) - 0.01*Intensity(j,1)*log(diam)
                   !                                Psurvival_s =1.0
                ELSEIF ((ht.LE.8.5).AND.(ht.GT.char_height)) THEN
                   Psurvival_s =(0.0178*Intensity(j,1) + 0.0144)*ht &
                        + (-0.1174*Intensity(j,1)+0.9158)

                ELSE
                   Psurvival_s = 0.0
                ENDIF
                Psurvival_s = MIN(Psurvival_s,1._dp)
                Psurvival_s = MAX(Psurvival_s,1e-3_dp)
                Psurvival = Psurvival_s

                IF (PRESENT(frac_intensity1)) THEN
                   char_height = 3.7*(1.-EXP(-0.19*Intensity(j,2)))
                   ht = pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%height
                   IF ((ht.GT.8.5).AND.(ht.GT.char_height)) THEN
                      Psurvival_s  =(-0.0011*Intensity(j,2) -0.00002)*ht &
                           +(0.0075*Intensity(j,2)+1.)
                      !Psurvival_s = 1+0.029*Intensity(j,2) - 0.01*Intensity(j,2)*log(diam)
                   ELSEIF ((ht.LE.8.5).AND.(ht.GT.char_height)) THEN
                      Psurvival_s =(0.0178*Intensity(j,2) + 0.0144)*ht &
                           + (-0.1174*Intensity(j,2)+0.9158)
                   ELSE
                      Psurvival_s = 0.0
                   ENDIF
                   Psurvival_s = MIN(Psurvival_s,1._dp)
                   Psurvival_s = MAX(Psurvival_s,1e-3_dp)
                   Psurvival = Psurvival_s*(1.-frac_intensity1(j)) + Psurvival*frac_intensity1(j)
                ENDIF

                !Psurvival = 0.0 ! test vh

                IF (it.EQ.10) THEN
                   !        write(333,*) ht, Psurvival
                ENDIF



                pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality + &
                     (1.-Psurvival)*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass = Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density = Psurvival*pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%density
                IF (pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density.LT.DENSINDIV_MIN) THEN
                   ! remove cohort
                   pop%pop_grid(j)%patch(k)%fire_mortality = pop%pop_grid(j)%patch(k)%fire_mortality + &
                        pop%pop_grid(j)%patch(k)%layer(1)%cohort(c)%biomass
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%density = 0.0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 0
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%biomass = 0.0
                ELSE
                   pop%pop_grid(j)%patch(k)%Layer(1)%cohort(c)%id = 1
                   nc = nc+1
                   ivec(nc)=c
                ENDIF
             ENDDO
             ! SHUFFLE if necessary to remove zero-density cohorts
             IF (nc.LT.pop%pop_grid(j)%patch(k)%Layer(1)%ncohort) THEN
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:nc)=pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ivec(1:nc))
                pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = nc

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(nc+1:NCOHORT_MAX)%biomass = 0.0
             ENDIF

             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

             ! understorey recruitment
             IF (PRESENT(precip)) THEN
                CALL layer_recruitment_single_patch(pop,k,j,precip)
             ELSE
                CALL layer_recruitment_single_patch(pop,k,j)
             ENDIF


          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE Patch_partial_disturb
  !*******************************************************************************
  SUBROUTINE Patch_disturb(pop, it,idisturb,precip)
    IMPLICIT NONE

    TYPE(POP_TYPE), INTENT(INOUT)  :: POP
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    INTEGER(i4b), INTENT(IN) :: it,idisturb
    INTEGER(i4b) :: j, k, np, nc

    np = SIZE(Pop%pop_grid)
    ! Kills all biomass in patch when prescribed disturbance interval is reached
    ! Should be called after accounting for this year
    DO j=1,np
       DO k=1,NPATCH2D
          IF (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).NE.0) THEN
             IF (pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb)) THEN
                ! kill entire layer
                nc = pop%pop_grid(j)%patch(k)%layer(1)%ncohort
                pop%pop_grid(j)%patch(k)%fire_mortality = SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%ncohort = 0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%density = 0.0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmean = 0.0
                pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmax  = 0.0
                pop%pop_grid(j)%patch(k)%age(idisturb) = 0
                pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%density = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%id = 0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass = 0.0
                pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%age = 0.0
                ! understorey recruitment
                IF (PRESENT(precip)) THEN
                   CALL layer_recruitment_single_patch(pop,k,j,precip)
                ELSE
                   CALL layer_recruitment_single_patch(pop,k,j)
                ENDIF
             ENDIF
          ELSEIF (pop%pop_grid(j)%patch(k)%disturbance_interval(idisturb).EQ.pop%pop_grid(j)%patch(k)%age(idisturb)) THEN
             ! kill entire layer
             nc = pop%pop_grid(j)%patch(k)%layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%fire_mortality = SUM(pop%pop_grid(j)%patch(k)%layer(1)%cohort(1:nc)%biomass)
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%ncohort = 0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%biomass = 0.0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%density = 0.0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmean = 0.0
             pop%pop_grid(j)%patch(k)%layer(1:NLayer)%hmax  = 0.0
             pop%pop_grid(j)%patch(k)%age(idisturb) = 0
             pop%pop_grid(j)%patch(k)%first_disturbance_year(idisturb) = 0

             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%density = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%id = 0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%biomass = 0.0
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(1:NCOHORT_MAX)%age = 0.0
             ! understorey recruitment
             IF (PRESENT(precip)) THEN
                CALL layer_recruitment_single_patch(pop,k,j,precip)
             ENDIF
          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE Patch_disturb
  !*******************************************************************************
  SUBROUTINE  layer_recruitment(pop,precip)
    IMPLICIT NONE
    TYPE(POP_TYPE), INTENT(INOUT)  :: POP
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    REAL(dp) :: f, mu, densindiv, cmass, ht
    REAL(dp) :: tmp, cmass_stem_sum,cmass_stem_inc
    INTEGER(i4b) :: j, k,c, ncohort, np
    REAL(dp) :: diam,basal

    np = SIZE(Pop%pop_grid)

    DO j=1,np
       !f = exp(-0.6*((pop%pop_grid(j)%cmass_sum)**(0.6667)))
       DO k=1,NPATCH2D
          IF(k.EQ.33) THEN
             ! write(*,*)
          ENDIF
          !write(*,*) k,  pop%pop_grid(j)%patch(k)%Layer(1)%biomass
          pop%pop_grid(j)%patch(k)%factor_recruit = EXP(-0.6*((pop%pop_grid(j)%patch(k)%Layer(1)%biomass)**(0.6667)))
          f = pop%pop_grid(j)%patch(k)%factor_recruit
          mu=EXP(FULTON_ALPHA*(1.0-2*THETA_recruit/(f+1-SQRT((f+1)*(f+1)-4*THETA_recruit*f))));
          densindiv=DENSINDIV_MAX*mu;
          cmass=CMASS_STEM_INIT*densindiv/DENSINDIV_MAX;
          IF (cmass>EPS*10..AND.densindiv>DENSINDIV_MIN.AND.(pop%pop_grid(j)%patch(k)%Layer(1)%ncohort+1).LT.NCOHORT_MAX) THEN
             ! create a new cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort + 1
             ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%biomass = cmass
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%density = densindiv

             ! ht=(Kbiometric**(3.0/5.0))*((2.0/PI)**(2.0/5.0))*((pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass/ &
             !pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density/WD)**(2.0/5.0))
             IF ((ALLOM_SWITCH.EQ.1).AND.(PRESENT(Precip))) THEN

                ht=GetHeight(precip(j),pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass, &
                     pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density)
                CALL Allometry(ht,pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass, &
                     pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density,diam,basal)


                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height = ht
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = (ht/Kbiometric)**(3/2)

             ELSEIF (ALLOM_SWITCH.EQ.0) THEN

                ht=(Kbiometric**(3.0/4.0))*(4.*pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass/ &
                     (pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density*WD*PI))**(1.0/4.0)

                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height = ht
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = (ht/Kbiometric)**(3/2)

             ELSE
                WRITE(*,*) 'Error: top-end allometry requires precip'
                STOP
             ENDIF
          ENDIF

          !if (k==1) then
          ! write(*,*) it, pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          !0 endif
       ENDDO
    ENDDO

  END SUBROUTINE layer_recruitment

  !*******************************************************************************
  SUBROUTINE  layer_recruitment_single_patch(pop, index, grid_index,precip)
    IMPLICIT NONE
    TYPE(POP_TYPE), INTENT(INOUT)  :: POP
    REAL(dp), INTENT(IN), OPTIONAL :: precip(:)
    INTEGER(i4b), INTENT(IN) :: index, grid_index
    REAL(dp) :: f, mu, densindiv, cmass, ht
    REAL(dp) :: tmp, cmass_stem_sum,cmass_stem_inc
    INTEGER(i4b) :: j, k,c, ncohort, np
    REAL(dp) :: diam,basal

    np = SIZE(Pop%pop_grid)
    DO j=grid_index,grid_index
       DO k=index,index
          pop%pop_grid(j)%patch(k)%factor_recruit = EXP(-0.6*((pop%pop_grid(j)%patch(k)%Layer(1)%biomass)**(0.6667)))
          f = pop%pop_grid(j)%patch(k)%factor_recruit
          mu=EXP(FULTON_ALPHA*(1.0-2*THETA_recruit/(f+1-SQRT((f+1)*(f+1)-4*THETA_recruit*f))));
          densindiv=DENSINDIV_MAX*mu;
          cmass=CMASS_STEM_INIT*densindiv/DENSINDIV_MAX;
          IF (cmass>EPS*10..AND.densindiv>DENSINDIV_MIN.AND.(pop%pop_grid(j)%patch(k)%Layer(1)%ncohort+1).LT.NCOHORT_MAX) THEN
             ! create a new cohort
             pop%pop_grid(j)%patch(k)%Layer(1)%ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort + 1
             ncohort = pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%biomass = cmass
             pop%pop_grid(j)%patch(k)%Layer(1)%cohort(ncohort)%density = densindiv

             !ht=(Kbiometric**(3.0/5.0))*((2.0/PI)**(2.0/5.0))*((pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass/ &
             ! pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density/WD)**(2.0/5.0))

             !ht=(Kbiometric**(3.0/4.0))*((4.*pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass/ &
             !pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density/WD/PI)**(1.0/4.0))
             IF ((ALLOM_SWITCH.EQ.1).AND.(PRESENT(Precip))) THEN

                ht=GetHeight(precip(j),pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass, &
                     pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density)
                CALL Allometry(ht,pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass, &
                     pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density,diam,basal)

                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height = ht
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = (ht/Kbiometric)**(3/2)

             ELSEIF (ALLOM_SWITCH.EQ.0) THEN
                ht=(Kbiometric**(3.0/4.0))*(4.*pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%biomass/ &
                     (pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%density*WD*PI))**(1.0/4.0)

                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%height = ht
                pop%pop_grid(j)%patch(k)%layer(1)%cohort(ncohort)%diameter = (ht/Kbiometric)**(3/2)
             ELSE
                WRITE(*,*) 'Error: top-end allometry requires precip'
                STOP
             ENDIF
          ENDIF

          !if (k==1) then
          ! write(*,*) it, pop%pop_grid(j)%patch(k)%Layer(1)%ncohort
          !0 endif
       ENDDO
    ENDDO

  END SUBROUTINE layer_recruitment_single_patch
  !*******************************************************************************
  REAL(dp) FUNCTION SuperPoisson(lambda, x)
    IMPLICIT NONE
    INTEGER(i4b), INTENT(IN) :: x
    REAL(dp), INTENT(IN) ::  lambda
    INTEGER(i4b), PARAMETER :: SWITCH_THRESHOLD=120 ! Switch to normal distribution for x above this value
    REAL(dp) :: p, sd
    REAL(dp) :: tmp

    IF (x.GT.SWITCH_THRESHOLD) THEN ! use normal distribution
       !SuperPoisson=1.0/(lambda**(0.5))/((PI*2.0)**0.5)*exp(-(real(x)-lambda)*(real(x)-lambda)/2.0/lambda)
       SuperPoisson=1.0/(lambda*(PI*2.0)**0.5)*EXP(-(REAL(x)-lambda)*(REAL(x)-lambda)/2.0/lambda)
    ELSE  ! use Poisson distribution
       tmp = factorial(x)
       SuperPoisson= (lambda**x)/(tmp)*EXP(-lambda)
    ENDIF

  END FUNCTION SuperPoisson
  !*******************************************************************************
  ! Exponential distribution (BEN)
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  REAL(dp) FUNCTION Exponential(lambda, x)
    IMPLICIT NONE
    INTEGER(i4b), INTENT(IN) :: x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       Exponential=0.0
    ELSE
       Exponential=lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION Exponential
  !*******************************************************************************
  ! Exponential distribution (BEN)
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  REAL(dp) FUNCTION REALExponential(lambda, x)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) ::  x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       REALExponential=0.0
    ELSE
       REALExponential=lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION REALExponential



  !*******************************************************************************
  REAL(dp) FUNCTION CumExponential(lambda, x)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       CumExponential=0.0
    ELSE
       CumExponential=1.-EXP(-lambda*x)
    ENDIF

  END FUNCTION CumExponential


  !*******************************************************************************

  REAL(dp)  FUNCTION Factorial(n)
    IMPLICIT NONE
    INTEGER , INTENT(IN) :: n
    INTEGER  :: i
    REAL(dp) :: Ans

    Ans = 1
    DO i = 1, n
       Ans = Ans *REAL(i,dp)
    END DO

    Factorial = Ans

  END FUNCTION Factorial
  !*******************************************************************************
  ! TOP-END ALLOMETRY STARTS HERE
  !*******************************************************************************
  ! Tree height based on precipitation and Gary Cook Top-End allometry (BEN)
  ! Bisection solution for tree height (m) based on modified height-DBH relationship
  ! from Garry Cook (pers. comm. 15/4/2013)
  ! "I have been using H=0.054xExp(0.0014xRF)xD + 4.05*exp(-0.00032*rf) with Rf in mm, D in cm and height in m."
  ! Since the above expression does not go to zero at diameter=0 it is linked to a simple linear equation
  ! for initial height growth (H=50*D with D in m) using a non-rectangular hyperbola to smooth between the two
  ! Mathematical derivation: see POP documentation
  ! Arguments:
  ! precip  = annual precipitation (mm)
  ! biomass = tree stem C biomass across patch (kgC/m2)
  ! density = tree density (indiv/m2)

  REAL(dp) FUNCTION GetHeight(precip,biomass,density)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: precip
    REAL(dp), INTENT(IN) :: biomass
    REAL(dp), INTENT(IN) :: density

    REAL(dp),PARAMETER:: THETA=0.99 ! Shape parameter, should be slightly <1
    REAL(dp),PARAMETER:: HMIN=0.001 ! min bound for tree height
    REAL(dp),PARAMETER:: HMAX=100 ! max bound for tree height
    REAL(dp),PARAMETER:: EPS=0.01 ! precision of the root
    INTEGER(i4b), PARAMETER :: MAXTRIES=25

    REAL(dp) :: alpha,beta,delta,rh,st,x1,x2,rtbis,dx,fmid,xmid,lhs,rhs
    INTEGER(i4b) :: b

    alpha=4.05*EXP(-0.00032*precip)
    beta=5.4*EXP(0.0014*precip)
    delta=2.0*SQRT(biomass/density/WD/PI)

    x1=HMIN
    x2=HMAX
    rtbis=x1
    dx=x2-x1
    b=0
    fmid=EPS+1.0

    DO WHILE (ABS(dx).GT.EPS.AND.b.LE.MAXTRIES)
       b=b+1
       dx=dx*0.5
       xmid=rtbis+dx

       ! Evaluate LHS-RHS at height=xmid
       ! LHS-RHS should increase with increasing height

       lhs=xmid
       rh=1.0/SQRT(xmid)
       st=alpha+beta*delta*rh+100*delta*rh
       rhs=1.0/2.0/THETA* &
            (st-SQRT(st*st-400*THETA*alpha*delta*rh- &
            400*THETA*beta*delta*delta/xmid))
       fmid=lhs-rhs

       IF (fmid.LT.0.0) rtbis=xmid

    ENDDO

    GetHeight=xmid

  END FUNCTION GetHeight
  !*******************************************************************************
  ! Top-End Allometry (BEN)
  ! Computes tree stem diameter (m) and basal area (m2/ha)
  ! given height (m), stem biomass (kgC/m2) and tree population density (indiv/m2)

  SUBROUTINE Allometry(height,biomass,density,diam,basal)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: height
    REAL(dp), INTENT(IN) :: biomass
    REAL(dp), INTENT(IN) :: density
    REAL(dp), INTENT(OUT) :: diam
    REAL(dp), INTENT(OUT) :: basal

    REAL(dp) :: delta,rh

    delta=2.0*SQRT(biomass/density/WD/PI)
    rh=1.0/SQRT(height)

    diam=delta*rh
    basal=PI*(diam/2.0)*(diam/2.0)*density*1e4

  END SUBROUTINE Allometry

  !*******************************************************************************

  SUBROUTINE POP_init(POP, disturbance_interval)
    USE POP_types, ONLY: POP_TYPE
    USE TypeDef, ONLY: i4b
    IMPLICIT NONE

    INTEGER(i4b), INTENT(IN) ::  disturbance_interval(:,:)
    TYPE( POP_TYPE )        , INTENT(INOUT) :: POP

    POP%it_pop = 0
    CALL ZeroPOP(pop)
    CALL InitPOP1D_Poisson(pop,INT(disturbance_interval))

  END SUBROUTINE POP_init

  !*******************************************************************************

!!!SUBROUTINE WRITE_POP_OUTPUT_NC ( POP, casamet, YEAR ) !,RST
!!!  USE netcdf
!!!  USE CASAVARIABLE
!!!  USE CABLE_COMMON_MODULE
!!!  IMPLICIT NONE
!!!
!!!  TYPE(POP_TYPE), INTENT(IN) :: POP
!!!  TYPE (casa_met),INTENT(IN) :: casamet
!!!  INTEGER       , INTENT(IN) :: YEAR
!!!!  LOGICAL       , INTENT(IN) :: RST
!!!
!!!  INTEGER                    :: STATUS,i,m,p,l,land_ID,patch_ID,ndis_ID
!!!  INTEGER                    :: ndis1_ID,nlay_ID,hgtb_ID,ncoh_ID,t_ID
!!!  CHARACTER(len=3)           :: typ = 'out'
!!!  CHARACTER                  :: dum*9,fname*100
!!!
!!!  ! 1 dim arrays (mp)
!!!  CHARACTER(len=40),DIMENSION( 2), PARAMETER :: AR0 = (/'latitude','longitude'/)
!!!
!!!  ! LANDSCAPE STRUCTURE
!!!  ! 2 dim arrays (mp,t)
!!!  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI1 = (/ 'npatch_active' /)
!!!  CHARACTER(len=40),DIMENSION(11), PARAMETER :: AR1 = (/ 'cmass_sum',           &
!!!       'densindiv','height_mean','height_max','basal_area','stress_mortality',  &
!!!       'fire_mortality','growth','crown_cover','crown_area','crown_volume' /)
!!!  ! 3 dim arrays (mp,nlayer,t)
!!!  CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR2 = (/ 'biomass','density',   &
!!!       'hmean','hmax' /)
!!!  ! 3 dim arrays (mp,height_bins,t)
!!!  CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR3 = (/ 'cmass_stem_bin',      &
!!!       'densindiv_bin','height_bin','diameter_bin' /)
!!!  ! 3 dim arrays (mp,ndisturb,t)
!!!  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI4 = (/ 'n_age' /)
!!!
!!!  ! PATCH STRUCTURE
!!!  ! 3 dim arrays (mp,npatch2d,t)
!!!  CHARACTER(len=40),DIMENSION(10), PARAMETER :: AR5 = (/ 'patch_freq',          &
!!!       'patch_freq_old','patch_freq_old2','patch_factor_recruit',               &
!!!       'patch_biomass','patch_biomass_old','patch_biomass_old2',                &
!!!       'patch_stress_mortality','patch_fire_mortality','patch_growth' /)
!!!  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI5 = (/ 'patch_id' /)
!!!  ! 4 dim arrays (mp,npatch2d,ndisturb+1,t)
!!!  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI6 =                           &
!!!       (/ 'patch_ranked_age_unique' /)
!!!  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AR6 =                           &
!!!       (/ 'patch_freq_ranked_age_unique' /)
!!!  ! 4 dim arrays (mp,npatch2d,ndisturb,t)
!!!  CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AI7 =                           &
!!!       (/ 'patch_disturbance_interval','patch_first_disturbance_year',          &
!!!       'patch_age','patch_age_old' /)
!!!
!!!  ! LAYER STRUCTURE
!!!  ! 4 dim arrays (mp,npatch2d,nlayer,t)
!!!  CHARACTER(len=40),DIMENSION( 4), PARAMETER :: AR8 = (/ 'layer_biomass',       &
!!!       'layer_density','layer_hmean','layer_hmax' /)
!!!  CHARACTER(len=40),DIMENSION( 1), PARAMETER :: AI8 = (/ 'layer_ncohort' /)
!!!
!!!  ! COHORT STRUCTURE
!!!  ! 5 dim arrays (mp,npatch2d,nlayer,ncohort_max,t)
!!!  CHARACTER(len=40),DIMENSION( 5), PARAMETER :: AR9 = (/ 'cohort_biomass',      &
!!!       'cohort_density','cohort_frac_resource_uptake','cohort_height',          &
!!!       'cohort_diameter' /)
!!!  CHARACTER(len=40),DIMENSION( 2), PARAMETER :: AI9 = (/'cohort_age','cohort_id'/)
!!!
!!!  LOGICAL, SAVE :: CALL1 = .TRUE.
!!!  INTEGER, SAVE :: VIDtime, VIDR0(SIZE(AR0)),VIDR1(SIZE(AR1)),VIDI1(SIZE(AI1)), &
!!!       VIDR2(SIZE(AR2)),VIDR3(SIZE(AR3)),VIDI4(SIZE(AI4)),VIDR5(SIZE(AR5)),     &
!!!       VIDI5(SIZE(AI5)),VIDR6(SIZE(AR6)),VIDI6(SIZE(AI6)),VIDI7(SIZE(AI7)),     &
!!!       VIDR8(SIZE(AR8)),VIDI8(SIZE(AI8)),VIDR9(SIZE(AR9)),VIDI9(SIZE(AI9))
!!!  INTEGER, SAVE :: FILE_ID, CNT = 0
!!!
!!!  CNT = CNT + 1
!!!
!!!  IF ( CALL1 ) THEN
!!!     ! Get File-Name
!!!     WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
!!!     fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
!!!          TRIM(dum)//'_pop_'//typ//'.nc'
!!!
!!!     ! Create NetCDF file:
!!!     STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!
!!!     ! Put the file in define mode:
!!!     STATUS = NF90_redef(FILE_ID)
!!!
!!!     ! GLOBAL ATTRIBUTES
!!!     STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle" , icycle             )
!!!     STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Year"   , YEAR               )
!!!     STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden", CABLE_USER%RunIden )
!!!
!!!     ! Define dimensions:
!!!     STATUS = NF90_def_dim(FILE_ID, 'land'       , mp         , land_ID )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'NPATCH2D'   , NPATCH2D   , patch_ID)
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'NDISTURB'   , NDISTURB   , ndis_ID )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'NDISTURB+1' , NDISTURB+1 , ndis1_ID)
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'NLAYER'     , NLAYER     , nlay_ID )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'HEIGHT_BINS', HEIGHT_BINS, hgtb_ID )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'NCOHORT_MAX', NCOHORT_MAX, ncoh_ID )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     STATUS = NF90_def_dim(FILE_ID, 'time'   ,  NF90_UNLIMITED, t_ID    )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!
!!!     ! Define variables
!!!     STATUS = NF90_def_var(FILE_ID,'Time' ,NF90_INT,(/t_ID/),VIDtime )
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!
!!!     DO i = 1, SIZE(AR0)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR0(i)), NF90_FLOAT,(/land_ID/),VIDR0(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI1)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)), NF90_INT  ,(/land_ID,t_ID/),VIDI1(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR1)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR1(i)), NF90_FLOAT,(/land_ID,t_ID/),VIDR1(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR2)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR2(i)), NF90_FLOAT,(/land_ID,nlay_ID,t_ID/),VIDR2(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR3)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR3(i)), NF90_FLOAT,(/land_ID,hgtb_ID,t_ID/),VIDR3(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI4)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI4(i)), NF90_INT  ,(/land_ID,ndis_ID,t_ID/),VIDI4(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI5)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI5(i)), NF90_INT  ,(/land_ID,patch_ID,t_ID/),VIDI5(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR5)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR5(i)), NF90_FLOAT,(/land_ID,patch_ID,t_ID/),VIDR5(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI6)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI6(i)), NF90_INT  , &
!!!             (/land_ID,patch_ID,ndis1_ID,t_ID/),VIDI6(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR6)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR6(i)), NF90_FLOAT, &
!!!             (/land_ID,patch_ID,ndis1_ID,t_ID/),VIDR6(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI7)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI7(i)), NF90_INT,   &
!!!             (/land_ID,patch_ID,ndis_ID,t_ID/),VIDI7(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI8)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI8(i)), NF90_INT  , &
!!!             (/land_ID,patch_ID,nlay_ID,t_ID/),VIDI8(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR8)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR8(i)), NF90_FLOAT, &
!!!             (/land_ID,patch_ID,nlay_ID,t_ID/),VIDR8(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AI9)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AI9(i)), NF90_INT  , &
!!!             (/land_ID,patch_ID,nlay_ID,ncoh_ID,t_ID/),VIDI9(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!     DO i = 1, SIZE(AR9)
!!!        STATUS = NF90_def_var(FILE_ID,TRIM(AR9(i)), NF90_FLOAT, &
!!!             (/land_ID,patch_ID,nlay_ID,ncoh_ID,t_ID/),VIDR9(i))
!!!        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!     END DO
!!!
!!!     ! End define mode:
!!!     STATUS = NF90_enddef(FILE_ID)
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!
!!!     ! PUT LAT / LON ( mp )
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR0(1), REAL(casamet%lat),      &
!!!          start=(/ 1 /), count=(/ mp /)  )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR0(2), REAL(casamet%lon),      &
!!!          start=(/ 1 /), count=(/ mp /)  )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!     CALL1 = .FALSE.
!!!
!!!  END IF
!!!
!!!  ! WRITE CURRENT STATE
!!!  ! TIME  ( t )
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, YEAR, start=(/ CNT /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!  ! PUT 2D VARS ( mp, t )
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1( 1), POP%pop_grid(:)%npatch_active,      &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 1), POP%pop_grid(:)%cmass_sum,          &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 2), POP%pop_grid(:)%densindiv,          &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 3), POP%pop_grid(:)%height_mean,        &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 4), POP%pop_grid(:)%height_max,         &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 5), POP%pop_grid(:)%basal_area,         &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 6), POP%pop_grid(:)%stress_mortality,   &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 7), POP%pop_grid(:)%fire_mortality,     &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 8), POP%pop_grid(:)%growth,             &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1( 9), POP%pop_grid(:)%crown_cover,        &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(10), POP%pop_grid(:)%crown_area,         &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!  STATUS = NF90_PUT_VAR(FILE_ID, VIDR1(11), POP%pop_grid(:)%crown_volume,       &
!!!       start=(/ 1, CNT /), count=(/ mp, 1 /) )
!!!  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!  ! PUT 3D VARS ( mp,nlayer, t )
!!!  MPS:DO m = 1, mp
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 1), POP%pop_grid(m)%biomass,         &
!!!          start=(/ m, 1, CNT /), count=(/ 1, nlayer, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 2), POP%pop_grid(m)%density,         &
!!!          start=(/ m, 1, CNT /), count=(/ 1, nlayer, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 3), POP%pop_grid(m)%hmean,           &
!!!          start=(/ m, 1, CNT /), count=(/ 1, nlayer, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR2( 4), POP%pop_grid(m)%hmax,            &
!!!          start=(/ m, 1, CNT /), count=(/ 1, nlayer, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!     ! PUT 3D VARS ( mp,height_bins, t )
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 1), POP%pop_grid(m)%cmass_stem_bin,  &
!!!          start=(/ m, 1, CNT /), count=(/ 1, HEIGHT_BINS, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 2), POP%pop_grid(m)%densindiv_bin,   &
!!!          start=(/ m, 1, CNT /), count=(/ 1, HEIGHT_BINS, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 3), POP%pop_grid(m)%height_bin,      &
!!!          start=(/ m, 1, CNT /), count=(/ 1, HEIGHT_BINS, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR3( 4), POP%pop_grid(m)%diameter_bin,    &
!!!          start=(/ m, 1, CNT /), count=(/ 1, HEIGHT_BINS, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!     ! PUT 3D VARS ( mp,ndisturb, t )
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDI4( 1), POP%pop_grid(m)%n_age,           &
!!!          start=(/ m, 1, CNT /), count=(/ 1, ndisturb, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!     ! PATCH STRUCTURE
!!!     ! PUT 3D VARS ( mp,npatch2d, t )
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDI5( 1), POP%pop_grid(m)%patch(:)%id,     &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 1), POP%pop_grid(m)%freq,            &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 2), POP%pop_grid(m)%freq_old,        &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 3), POP%pop_grid(m)%freq_old2,       &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 4), POP%pop_grid(m)%patch(:)%factor_recruit,  &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 5), POP%pop_grid(m)%patch(:)%biomass,         &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 6), POP%pop_grid(m)%patch(:)%biomass_old,     &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 7), POP%pop_grid(m)%patch(:)%biomass_old2,    &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 8), POP%pop_grid(m)%patch(:)%stress_mortality,&
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5( 9), POP%pop_grid(m)%patch(:)%fire_mortality,  &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!     STATUS = NF90_PUT_VAR(FILE_ID, VIDR5(10), POP%pop_grid(m)%patch(:)%growth,          &
!!!          start=(/ m, 1, CNT /), count=(/ 1, npatch2d, 1 /) )
!!!     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!     PAT:DO p = 1, npatch2d
!!!        ! PUT 4D VARS ( mp,npatch2d, ndisturb+1,t )
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDI6( 1), POP%pop_grid(m)%ranked_age_unique(p,:),&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB+1, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDR6( 1), POP%pop_grid(m)%freq_ranked_age_unique(p,:),&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB+1, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!        ! PUT 4D VARS ( mp,npatch2d, ndisturb,t )
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 1), POP%pop_grid(m)%patch(p)%disturbance_interval,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 2), POP%pop_grid(m)%patch(p)%first_disturbance_year,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 3), POP%pop_grid(m)%patch(p)%age,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDI7( 4), POP%pop_grid(m)%patch(p)%age_old,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, NDISTURB, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!        ! LAYER STRUCTURE
!!!        ! PUT 4D VARS ( mp,npatch2d, nlayer,t )
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDI8( 1), POP%pop_grid(m)%patch(p)%layer(:)%ncohort,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, nlayer, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 1), POP%pop_grid(m)%patch(p)%layer(:)%biomass,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, nlayer, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 2), POP%pop_grid(m)%patch(p)%layer(:)%density,&
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, nlayer, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 3), POP%pop_grid(m)%patch(p)%layer(:)%hmean,  &
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, nlayer, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!        STATUS = NF90_PUT_VAR(FILE_ID, VIDR8( 4), POP%pop_grid(m)%patch(p)%layer(:)%hmax,   &
!!!             start=(/ m, p, 1, CNT /), count=(/ 1, 1, nlayer, 1 /) )
!!!        IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!        LAY:DO l = 1, nlayer
!!!           ! COHORT STRUCTURE
!!!           ! PUT 5D VARS ( mp,npatch2d, nlayer,ncohort_max,t )
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( 1), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%age,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDI9( 2), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%id,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 1), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%biomass,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 2), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%density,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 3), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%frac_resource_uptake,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 4), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%height,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!           STATUS = NF90_PUT_VAR(FILE_ID, VIDR9( 5), POP%pop_grid(m)%patch(p)%layer(l)%cohort(:)%diameter,&
!!!                start=(/ m, p, l, 1, CNT /), count=(/ 1, 1, 1, ncohort_max, 1 /) )
!!!           IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!!
!!!        END DO LAY
!!!     END DO PAT
!!!  END DO MPS
!!!
!!!  IF ( YEAR .EQ. cable_user%YEAREND ) THEN
!!!     ! Close NetCDF file:
!!!     STATUS = NF90_close(FILE_ID)
!!!     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!!  ENDIF
!!!
!!!END SUBROUTINE WRITE_POP_OUTPUT_NC
!!!
!!!
!!!
!!!
!!!
!!!
!!!

  !*******************************************************************************

  !CLNMODULE POPvariable
  !CLN
  !CLN    IMPLICIT NONE
  !CLN
  !CLNCONTAINS
  !CLN
  !*******************************************************************************

  SUBROUTINE alloc_POP(POP, arraysize)

    USE TypeDef,   Only: i4b, dp
    USE POP_Types, Only: POP_TYPE
    !    USE POPModule, Only: InitPOP2D_Poisson,InitPOP1D_Poisson, ZeroPOP

    TYPE( POP_TYPE ),INTENT(INOUT) :: POP
    INTEGER,            INTENT(IN) :: arraysize

    IF (.NOT.ALLOCATED(POP%POP_Grid)) ALLOCATE (POP%POP_Grid(arraysize))

  END SUBROUTINE alloc_POP

  !CLNEnd MODULE POPvariable

  !*******************************************************************************
END MODULE POPModule
!*******************************************************************************
