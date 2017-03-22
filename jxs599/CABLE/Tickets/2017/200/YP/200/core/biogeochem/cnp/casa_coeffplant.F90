MODULE casa__mod

USE cable_def_types_mod
USE casadimension
USE casaparm
USE casavariable
USE phenvariable
USE cable_common_module, only: cable_user ! Custom soil respiration: Ticket #42

IMPLICIT NONE
  REAL(r_2), PARAMETER :: zero = 0.0_r_2
  REAL(r_2), PARAMETER :: one  = 1.0_r_2

CONTAINS

SUBROUTINE casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                           casaflux,casamet,phen)
! calculate the plant litter fall rate, litter fall and sOM decomposition rate (1/day)
! and the transfer coefficients between different pools
!
! inputs:
!     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
!     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
!     xkleaf(mp):      set to 0.0 during maximal leaf growth phase
!
! outputs:
!     kplant(mp,mplant):        senescence rate of plant pool (1/day)
!     fromPtoL(mp,mlitter,mplant): fraction of senesced plant biomass to litter pool (fraction)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(IN)    :: xkleafcold,xkleafdry,xkleaf
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (phen_variable),       INTENT(IN) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mp)  :: xk
  REAL(r_2), DIMENSION(mp,mplant)         :: ratioLignintoN
  INTEGER npt

  casaflux%fromPtoL(:,:,:)      = 0.0
  casaflux%kplant(:,:)          = 0.0   ! (BPjun2010)

  WHERE(casamet%iveg2/=icewater)
  ! using max function to avoid dividing by zero, ypw 14/may/2008
    ratioLignintoN(:,leaf) = (casapool%Cplant(:,leaf) &
                             /(max(1.0e-10,casapool%Nplant(:,leaf)) *casabiome%ftransNPtoL(veg%iveg(:),leaf))) &
                             * casabiome%fracLigninplant(veg%iveg(:),leaf)
    ratioLignintoN(:,froot)= (casapool%Cplant(:,froot)&
                             /(max(1.0e-10,casapool%Nplant(:,froot))*casabiome%ftransNPtoL(veg%iveg(:),froot))) &
                             * casabiome%fracLigninplant(veg%iveg(:),froot)

    casaflux%fromPtoL(:,metb,leaf)    = max(0.001, 0.85 - 0.018 *ratioLignintoN(:,leaf))
    casaflux%fromPtoL(:,metb,froot)   = max(0.001, 0.85 - 0.018 *ratioLignintoN(:,froot))
    casaflux%fromPtoL(:,str,leaf)    = 1.0 - casaflux%fromPtoL(:,metb,leaf)
    casaflux%fromPtoL(:,str,froot)   = 1.0 - casaflux%fromPtoL(:,metb,froot)
    casaflux%fromPtoL(:,cwd,wood)    = 1.0

    casaflux%kplant(:,leaf)        = casabiome%plantrate(veg%iveg(:),leaf)*xkleaf(:) &
                                   + xkleafcold(:) + xkleafdry(:)

    casaflux%kplant(:,wood)        = casabiome%plantrate(veg%iveg(:),wood)
    casaflux%kplant(:,froot)       = casabiome%plantrate(veg%iveg(:),froot)
  ENDWHERE


  ! When glai<glaimin,leaf biomass will not decrease anymore. (Q.Zhang 10/03/2011)
  DO npt = 1,mp
    if(casamet%glai(npt).le.casabiome%glaimin(veg%iveg(npt))) casaflux%kplant(npt,leaf) = 0.0
  ENDDO
  ! end change

END SUBROUTINE casa_coeffplant


END MODULE casa__mod
