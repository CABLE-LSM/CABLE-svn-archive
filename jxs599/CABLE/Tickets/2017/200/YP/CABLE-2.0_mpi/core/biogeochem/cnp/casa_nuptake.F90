SUBROUTINE casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
! (1) compute (1)N uptake by plants; 
! (2) allocation of uptaken N to plants 
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
  REAL(r_2), DIMENSION(mp,mplant)      :: Nreqmax, Nreqmin, NtransPtoP, xnuptake
  REAL(r_2), DIMENSION(mp)             :: totNreqmax,totNreqmin
  REAL(r_2), DIMENSION(mp)             :: xnCnpp

  Nreqmin(:,:)       = 0.0
  Nreqmax(:,:)       = 0.0
  NtransPtoP(:,:)    = 0.0
  totNreqmax = 0.0
  totNreqmin = 0.0

  casaflux%Nminuptake(:)     = 0.0
  casaflux%fracNalloc(:,:)   = 0.0
  xnCnpp = max(0.0_r_2,casaflux%Cnpp)
  call casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
  
  DO np=1,mp
  IF(casamet%iveg2(np)/=icewater) THEN 
    totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
    totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)

    xnuptake(np,leaf) = Nreqmin(np,leaf) + xkNlimiting(np)* (Nreqmax(np,leaf)-Nreqmin(np,leaf))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))
    xnuptake(np,wood) = Nreqmin(np,wood) + xkNlimiting(np)* (Nreqmax(np,wood)-Nreqmin(np,wood))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))
    xnuptake(np,froot) = Nreqmin(np,froot) + xkNlimiting(np)* (Nreqmax(np,froot)-Nreqmin(np,froot))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))

    casaflux%Nminuptake(np) = xnuptake(np,leaf) + xnuptake(np,wood) + xnuptake(np,froot)+1.0e-10
    casaflux%fracNalloc(np,leaf)  = xnuptake(np,leaf)/casaflux%Nminuptake(np)
    casaflux%fracNalloc(np,wood)  = xnuptake(np,wood)/casaflux%Nminuptake(np)
    casaflux%fracNalloc(np,froot) = xnuptake(np,froot)/casaflux%Nminuptake(np)
  ENDIF
  ENDDO

!  np=1
!  write(*,911) casapool%nsoilmin(np),casaflux%Nminuptake(np),xnuptake(np,leaf), xnuptake(np,wood), xnuptake(np,froot), &
!               casaflux%fracNalloc(np,leaf),casaflux%fracNalloc(np,wood),casaflux%fracNalloc(np,froot)

  casaflux%Nupland = casaflux%Nminuptake

911 format('N uptake:',100(f8.3,2x))
END SUBROUTINE casa_nuptake


