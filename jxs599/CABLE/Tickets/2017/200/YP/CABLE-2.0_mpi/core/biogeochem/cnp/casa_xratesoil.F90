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
!,tsurfavg  !!, msurfavg
  INTEGER :: npt

!  print *, 'within casa_xratesoil'

  xklitter(:) = 1.0
  xksoil(:)   = 1.0
  fwps(:)     =  casamet%moistavg(:)/soil%ssat(:)
  tsavg(:)    =  casamet%tsoilavg(:) 

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
    xksoil(npt)   = casabiome%xkoptsoil(veg%iveg(npt))   * xktemp(npt) * xkwater(npt)
  END IF
  END DO
!  WHERE(casamet%iveg2/=icewater)  
!!    ! Kirschbaum function
!!    xktemp(:) = exp(xkalpha + xkbeta*(tsavg(:)-TKzeroC) &
!!              * (1.0-0.5*(tsavg(:)-TKzeroc)/xktoptc))
!    ! add by ypwang on 3/april/2009
!!    xktemp(:) = xkoptcoeff(veg%iveg(:))*exp(xkbeta*(tsavg(:)-TKzeroC-xktoptc))
!    xktemp(:)  = q10soil**(0.1*(tsavg(:)-TKzeroC-35.0))
!    xkwater(:) = ((fwps(:)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
!               * ((fwps(:)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
!    WHERE(veg%iveg==12)
!      xkwater(:)=1.0
!    ENDWHERE
!    xklitter(:) = xkoptlitter(veg%iveg(:)) * xktemp(:) * xkwater(:)
!    xksoil(:)   = xkoptsoil(veg%iveg(:))   * xktemp(:) * xkwater(:)
!  ENDWHERE
!   npt =26493
!  print *, 'xratesoil', npt, casamet%moistavg(npt),soil%ssat(npt),casamet%tsoilavg(npt),xklitter(npt),tsavg(npt),xksoil(npt)

END SUBROUTINE casa_xratesoil


