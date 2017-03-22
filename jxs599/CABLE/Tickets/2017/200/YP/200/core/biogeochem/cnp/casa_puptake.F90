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

SUBROUTINE casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
! (1) compute  P uptake by plants;
! (2) allocation of uptaken P to plants
!
  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  REAL(r_2), DIMENSION(mp),     INTENT(IN)    :: xkNlimiting


  ! local variables
  INTEGER                   nland,np,ip
  REAL(r_2), DIMENSION(mp,mplant) :: Preqmax,Preqmin,PtransPtoP,xPuptake
  REAL(r_2), DIMENSION(mp)        :: totPreqmax,totPreqmin
  REAL(r_2), DIMENSION(mp)        :: xpCnpp

  Preqmin(:,:)             = 0.0
  Preqmax(:,:)             = 0.0
  PtransPtoP(:,:)          = 0.0
  casaflux%Plabuptake(:)   = 0.0
  casaflux%fracPalloc(:,:) = 0.0
  totPreqmax               = 0.0
  totPreqmin               = 0.0

  xpCnpp = max(0.0_r_2,casaflux%cnpp)
  call casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
  WHERE(casamet%iveg2/=icewater)
    totPreqmax(:) = Preqmax(:,leaf)+Preqmax(:,wood)+Preqmax(:,froot)
    totPreqmin(:) = Preqmin(:,leaf)+Preqmin(:,wood)+Preqmin(:,froot)

    xpuptake(:,leaf) = Preqmin(:,leaf) + xkNlimiting(:)* (Preqmax(:,leaf)-Preqmin(:,leaf))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))
    xpuptake(:,wood) = Preqmin(:,wood) + xkNlimiting(:)* (Preqmax(:,wood)-Preqmin(:,wood))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))
    xpuptake(:,froot) = Preqmin(:,froot) + xkNlimiting(:)* (Preqmax(:,froot)-Preqmin(:,froot))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))

    casaflux%Plabuptake(:) = xpuptake(:,leaf) + xpuptake(:,wood) + xpuptake(:,froot)+1.0e-10
    casaflux%fracPalloc(:,leaf)  = xpuptake(:,leaf)/casaflux%Plabuptake(:)
    casaflux%fracPalloc(:,wood)  = xpuptake(:,wood)/casaflux%Plabuptake(:)
    casaflux%fracPalloc(:,froot) = xpuptake(:,froot)/casaflux%Plabuptake(:)

  ENDWHERE

  casaflux%Pupland = casaflux%Plabuptake

!  np=1
!  write(*,911) casapool%Psoillab(np),casaflux%Plabuptake(np), &
!               xpuptake(np,leaf), xpuptake(np,wood), xpuptake(np,froot), &
!               casaflux%fracPalloc(np,leaf),casaflux%fracPalloc(np,wood), &
!               casaflux%fracPalloc(np,froot)
!911 format('P uptake:',100(f8.3,2x))

!  ! only used in spinning up the model
!  DO  np=1,mp
!    casaflux%Plabuptake(np) = TotPreqmax(np)
!    casaflux%Pupland(np)    = TotPreqmax(np)
!    casaflux%Pwea(np)       = TotPreqmax(np)
!  ENDDO

END SUBROUTINE casa_puptake


END MODULE casa__mod
