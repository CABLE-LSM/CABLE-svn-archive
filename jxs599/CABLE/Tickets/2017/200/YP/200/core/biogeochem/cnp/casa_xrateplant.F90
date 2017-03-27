SUBROUTINE casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                           casamet,phen)
! use xleafcold and xleafdry to account for
! cold and drought stress on death rate of leaf
! inputs:
!     ivt(mp)  :       biome type
!     phase(mp):       leaf growth stage
!     tairk(mp)    :   air temperature in K
! outputs
!     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
!     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
!     xkleaf(mp):      set to 0.0 during maximal leaf growth phase

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleafcold
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleafdry
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleaf
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mp)              :: xcoldleaf
  INTEGER :: npt

  xkleafcold(:) = 0.0
  xkleafdry(:)  = 0.0
  xkleaf(:)     = 1.0

  ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
  !    following the formulation of Arora (2005) on the
  !    effect of cold or drought stress on leaf litter fall
  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
    IF(casamet%tairk(npt)>=phen%TKshed(veg%iveg(npt))) THEN
      xcoldleaf(npt) = 1.0
    ELSE
      IF(casamet%tairk(npt)<=(phen%TKshed(veg%iveg(npt))-5.0)) THEN
        xcoldleaf(npt)=0.0
      ELSE
        xcoldleaf(npt) = (casamet%tairk(npt)-phen%TKshed(veg%iveg(npt))-5.0)/5.0
      ENDIF
    ENDIF
    xcoldleaf(npt) = min(1.0,max(0.0,xcoldleaf(npt)))
    xkleafcold(npt) = casabiome%xkleafcoldmax(veg%iveg(npt)) &
                    * (1.0-xcoldleaf(npt)) &
                    ** casabiome%xkleafcoldexp(veg%iveg(npt))
    xkleafdry(npt)  = casabiome%xkleafdrymax(veg%iveg(npt)) &
                    * (1.0-casamet%btran(npt))&
                    ** casabiome%xkleafdryexp(veg%iveg(npt))
    IF (phen%phase(npt)==1) xkleaf(npt)= 0.0
    !Ticket200
    ! vh: account for high rate of leaf loss during senescence
    if (trim(cable_user%PHENOLOGY_SWITCH)=='climate') then
       IF (phen%phase(npt)==3.or.phen%phase(npt)==0) xkleaf(npt)= 100.0
    endif
  END IF
  END DO

!  WHERE(casamet%iveg2/=icewater)
!  !    following the formulation of Arora (2005) on the
!  !    effect of cold or drought stress on leaf litter fall
!  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
!    WHERE(casamet%tairk(:)>=phen%TKshed(veg%iveg(:)))
!      xcoldleaf(:) = 1.0
!    ELSEWHERE
!      WHERE(casamet%tairk(:)<=(phen%TKshed(veg%iveg(:))-5.0))
!        xcoldleaf(:)=0.0
!      ELSEWHERE
!        xcoldleaf(:) = (casamet%tairk(:)-phen%TKshed(veg%iveg(:))-5.0)/5.0
!      ENDWHERE
!    ENDWHERE
!    xcoldleaf(:) = min(1.0,max(0.0,xcoldleaf(:)))
!    xkleafcold(:) = casabiome%xkleafcoldmax(veg%iveg(:)) * (1.0-xcoldleaf(:)) &
!                 ** casabiome%xkleafcoldexp(veg%iveg(:))
!    xkleafdry(:) = casabiome%xkleafdrymax(veg%iveg(:))*(1.0-casamet%btran(:))&
!                 ** casabiome%xkleafdryexp(veg%iveg(:))
!    WHERE(phen%phase(:)==1) xkleaf(:)= 0.0
!  ENDWHERE

END SUBROUTINE casa_xrateplant


