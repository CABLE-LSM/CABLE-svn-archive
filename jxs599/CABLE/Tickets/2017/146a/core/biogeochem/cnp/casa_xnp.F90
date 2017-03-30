SUBROUTINE casa_xNP(xNPlimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xnplimit
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xNPuptake
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  INTEGER :: np
  REAL(r_2), DIMENSION(mp)        :: xnlimit,xplimit
  REAL(r_2), DIMENSION(mp)        :: xncleaf,xpcleaf
  REAL(r_2), DIMENSION(mp)        :: xnCnpp,xpCnpp
  REAL(r_2), DIMENSION(mp,mplant) :: Nreqmax, Nreqmin, NtransPtoP
  REAL(r_2), DIMENSION(mp)        :: totNreqmax,totNreqmin
  REAL(r_2), DIMENSION(mp)        :: xNuptake,xPuptake
  REAL(r_2), DIMENSION(mp,mplant) :: Preqmax, Preqmin, PtransPtoP
  REAL(r_2), DIMENSION(mp)        :: totPreqmax,totPreqmin

  xnlimit  = 1.0
  xplimit  = 1.0
  xnplimit = 1.0
  casaflux%fracClabile(:) = 0.0

  SELECT CASE(icycle)
  CASE(2)
    WHERE(casamet%iveg2/=icewater)
      xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
      !xnlimit(:) = xncleaf(:)/(xncleaf(:)+casabiome%KminN(veg%iveg(:)))
      xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01)
      xplimit(:) = 1.0
      xnplimit(:) =min(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
    ENDWHERE
  CASE(3)
    WHERE(casamet%iveg2/=icewater)
      xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
      xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01)
      xpcleaf(:) = casapool%pplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
      !xplimit(:) = xpcleaf(:)/(xpcleaf(:)+casabiome%Kuplabp(veg%iveg(:)))
      xplimit(:) = xpcleaf(:)/(xpcleaf(:)+0.0006)
      !xnplimit(:) = min(1.0,casabiome%Kuptake(veg%iveg(:))*min(xnlimit(:),xplimit(:)))
      xnplimit(:) =min(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
    ENDWHERE
  END SELECT

  ! now check if soil nutrient supply can meet the plant uptake,
  ! otherwise reduce NPP
  xNuptake = 1.0
  xPuptake = 1.0

  IF(icycle >1) THEN
    Nreqmin(:,:)    = 0.0
    Nreqmax(:,:)    = 0.0
    NtransPtoP(:,:) = 0.0
    totNreqmax = 0.0
    totNreqmin = 0.0
    xNuptake   = 1.0

    xnCnpp = max(0.0,casaflux%Cnpp)
    call casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
    DO np=1,mp
      IF(casamet%iveg2(np)/=icewater) THEN
        totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
        totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)
        xNuptake(np)   = MAX(0.0,MIN(1.0,casapool%Nsoilmin(np) &
                                         /(totNreqmin(np)*deltpool+1.0e-10)))
      ENDIF
    ENDDO
  ENDIF
  IF(icycle >2) THEN
    Preqmin(:,:)       = 0.0
    Preqmax(:,:)       = 0.0
    PtransPtoP(:,:)    = 0.0
    totPreqmax = 0.0
    totPreqmin = 0.0
    xPuptake   = 1.0
    xpCnpp = max(0.0,casaflux%Cnpp)
    call casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                       casabiome,casapool,casaflux,casamet)
    DO np=1,mp
      IF(casamet%iveg2(np)/=icewater) THEN
        totPreqmax(np) = Preqmax(np,leaf)+Preqmax(np,wood)+Preqmax(np,froot)
        totPreqmin(np) = Preqmin(np,leaf)+Preqmin(np,wood)+Preqmin(np,froot)
        xPuptake(np)   = MAX(0.0,MIN(1.0,casapool%psoillab(np) &
                                         /(totPreqmin(np)*deltpool+1.0e-10)))
      ENDIF
    ENDDO
  ENDIF

  xnplimit(:)  = 1.0
  xNPuptake(:)     = min(xnuptake(:), xpuptake(:))
  do np =1, mp
     if(casamet%iveg2(np)/=icewater.and.casaflux%cnpp(np) > 0.0.and.xNPuptake(np) < 1.0) then
        casaflux%fracClabile(np) =min(1.0,max(0.0,(1.0- xNPuptake(np)))) * max(0.0,casaflux%cnpp(np)) &
                                 /(casaflux%cgpp(np) +1.0e-10)
        casaflux%cnpp(np)    = casaflux%cnpp(np) - casaflux%fracClabile(np) * casaflux%cgpp(np)
     endif

   enddo

!  casaflux%cnpp(:) = xNPuptake(:) * xnplimit(:) * casaflux%cnpp(:)


END SUBROUTINE casa_xnp

