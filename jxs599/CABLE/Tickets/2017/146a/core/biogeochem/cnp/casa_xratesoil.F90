SUBROUTINE casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
!  to account for cold and drought stress on death rate of leaf: xleafcold,xleafdry
!  to account for effects of T and W on litter decomposition: xk, xksurf
!  inputs:
!     ivt(mp)  :       biome type
!     tsoilavg(mp):    soil temperature in K
!     moistavg(mp):    volumetric soil moisture
!
!  outputs
!     xk(mp):          modifier of soil litter decomposition rate (dimensionless)
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xklitter,xksoil
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome

  ! local variables
  INTEGER nland,np
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)
! Kirschbaum function parameters
  REAL(r_2), parameter :: xkalpha=-3.764   ! Kirschbaum (1995, SBB)
  REAL(r_2), parameter :: xkbeta=0.204
  REAL(r_2), parameter :: xktoptc=36.9
  REAL(r_2), DIMENSION(mp)       :: xkwater,xktemp
  REAL(r_2), DIMENSION(mp)       :: fwps,tsavg
  ! Custom soil respiration - see Ticket #42
  REAL(r_2), DIMENSION(mp)       :: smrf,strf,slopt,wlt,tsoil,fcap,sopt
!,tsurfavg  !!, msurfavg
  INTEGER :: npt

  xklitter(:) = 1.0
  xksoil(:)   = 1.0
  fwps(:)     =  casamet%moistavg(:)/soil%ssat(:)
  tsavg(:)    =  casamet%tsoilavg(:)

  ! Custom soil respiration - see Ticket #42
  tsoil(:)    =  tsavg(:)-TKzeroC !tsoil in C
  strf(:)     = 1.0
  smrf(:)     = 1.0
  slopt(:)    = 1.0
  sopt(:)     = 1.0


  ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
    ! Kirschbaum function
    xktemp(npt)  = casabiome%q10soil(veg%iveg(npt))**(0.1*(tsavg(npt)-TKzeroC-35.0))
    xkwater(npt) = ((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
               * ((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
    IF (veg%iveg(npt) == cropland .OR. veg%iveg(npt) == croplnd2) &
               xkwater(npt)=1.0
    xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) * xktemp(npt) * xkwater(npt)

    IF( .NOT. cable_user%SRF) THEN
        ! Use original function, ELSE Ticket #42
       xksoil(npt)   = casabiome%xkoptsoil(veg%iveg(npt))   * xktemp(npt) * xkwater(npt)
    ELSE
    ! Custom soil respiration - see Ticket #42
    ! Implementing alternative parameterizations
      IF(trim(cable_user%SMRF_NAME)=='CASA-CNP') THEN
         smrf(npt)=xkwater(npt)
      ELSE IF (trim(cable_user%SMRF_NAME)=='SOILN') then
         sopt(npt)=0.92
         slopt(npt)=wlt(npt)+0.1          !SLOPT is the lower optimum
         IF (fwps(npt)>sopt(npt)) THEN
           smrf(npt)=0.2+0.8*(1.0-fwps(npt))/(1.0-sopt(npt))
         ELSE IF(slopt(npt)<=fwps(npt) .AND. fwps(npt)<=sopt(npt)) THEN
           smrf(npt) = 1.0
         ELSE IF (wlt(npt)<=fwps(npt) .AND. fwps(npt) <slopt(npt)) THEN
           smrf(npt)=0.01+0.99*(fwps(npt)-wlt(npt))/(slopt(npt)-wlt(npt))
         ELSE IF (fwps(npt)<wlt(npt)) THEN
           smrf(npt) = 0.01
         END IF
      ELSE IF (trim(cable_user%SMRF_NAME)=='TRIFFID') THEN
         sopt(npt) = 0.5 * (1+wlt(npt))
         IF (fwps(npt) > sopt(npt)) THEN
           smrf(npt) =1.0-0.8*(fwps(npt)-sopt(npt))
         ELSE IF (wlt(npt)<fwps(npt) .AND. fwps(npt)<=sopt(npt)) THEN
           smrf(npt)=0.01+0.8*((fwps(npt)-wlt(npt))/(sopt(npt)-wlt(npt)))
         ELSE IF (fwps(npt)<wlt(npt)) THEN
           smrf(npt) = 0.2
         END IF
      END IF

      IF(trim(cable_user%STRF_NAME)=='CASA-CNP') THEN
        strf(npt)=xktemp(npt)
      ELSE if (trim(cable_user%STRF_NAME)=='K1995') THEN
      !Kirschbaum from Kirschbaum 1995, eq (4) in SBB, .66 is to collapse smrf
      !to same area
        strf(npt)=exp(-3.764+0.204*tsoil(npt)*(1-0.5*tsoil(npt)/36.9))/.66
      ELSE IF (trim(cable_user%STRF_NAME)=='PnET-CN') THEN
        strf(npt)=0.68*exp(0.1*(tsoil(npt)-7.1))/12.64
      END IF
      xksoil(npt) = casabiome%xkoptsoil(veg%iveg(npt))*strf(npt)*smrf(npt)
    END IF
  END IF
  END DO

END SUBROUTINE casa_xratesoil

