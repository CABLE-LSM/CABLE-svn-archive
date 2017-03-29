SUBROUTINE casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                         casabiome,casapool,casaflux,casamet)
!
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xnCnpp
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Nreqmax, Nreqmin
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: NtransPtoP
  TYPE (veg_parameter_type),             INTENT(INOUT) :: veg
  TYPE (casa_biome),                     INTENT(INOUT) :: casabiome
  TYPE (casa_pool),                      INTENT(INOUT) :: casapool
  TYPE (casa_flux),                      INTENT(INOUT) :: casaflux
  TYPE (casa_met),                       INTENT(INOUT) :: casamet

  ! local variable
  INTEGER :: np
  REAL(r_2), DIMENSION(mp,mplant)     :: ncplantmax

  Nreqmin(:,:)    = 0.0
  Nreqmax(:,:)    = 0.0
  NtransPtoP(:,:) = 0.0

  DO np=1,mp
  IF(casamet%iveg2(np)/=icewater) THEN
    if(casapool%Nsoilmin(np)<2.0) then
       ncplantmax(np,leaf) =casabiome%ratioNCplantmin(veg%iveg(np),leaf)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),leaf)-casabiome%ratioNCplantmin(veg%iveg(np),leaf)) &
                             * min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
       ncplantmax(np,wood) =casabiome%ratioNCplantmin(veg%iveg(np),wood)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),wood)-casabiome%ratioNCplantmin(veg%iveg(np),wood)) &
                             * min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
       ncplantmax(np,froot) =casabiome%ratioNCplantmin(veg%iveg(np),froot)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),froot)-casabiome%ratioNCplantmin(veg%iveg(np),froot)) &
                             * min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
    else
      ncplantmax(np,leaf)  = casabiome%ratioNCplantmax(veg%iveg(np),leaf)
      ncplantmax(np,wood)  = casabiome%ratioNCplantmax(veg%iveg(np),wood)
      ncplantmax(np,froot) = casabiome%ratioNCplantmax(veg%iveg(np),froot)
    endif

    Nreqmax(np,leaf)  = xnCnpp(np)* casaflux%fracCalloc(np,leaf) *ncplantmax(np,leaf)
    Nreqmax(np,wood)  = xnCnpp(np)* casaflux%fracCalloc(np,wood) *ncplantmax(np,wood)
    Nreqmax(np,froot) = xnCnpp(np)* casaflux%fracCalloc(np,froot)*ncplantmax(np,froot)

    Nreqmin(np,leaf) =  xnCnpp(np)* casaflux%fracCalloc(np,leaf) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),leaf)
    Nreqmin(np,wood) =  xnCnpp(np)* casaflux%fracCalloc(np,wood) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),wood)
    Nreqmin(np,froot) =  xnCnpp(np)* casaflux%fracCalloc(np,froot) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),froot)

    NtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Nplant(np,leaf) &
                       * (1.0-casabiome%ftransNPtoL(veg%iveg(np),leaf))
    NtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Nplant(np,wood) &
                       * (1.0-casabiome%ftransNPtoL(veg%iveg(np),wood))
    NtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Nplant(np,froot) &
                       * (1.0-casabiome%ftransNPtoL(veg%iveg(np),froot))

    Nreqmax(np,leaf)  = max(0.0,Nreqmax(np,leaf) - NtransPtoP(np,leaf))
    Nreqmax(np,wood)  = max(0.0,Nreqmax(np,wood) - NtransPtoP(np,wood))
    Nreqmax(np,froot) = max(0.0,Nreqmax(np,froot) - NtransPtoP(np,froot))
    Nreqmin(np,leaf)  = max(0.0,Nreqmin(np,leaf) - NtransPtoP(np,leaf))
    Nreqmin(np,wood)  = max(0.0,Nreqmin(np,wood) - NtransPtoP(np,wood))
    Nreqmin(np,froot) = max(0.0,Nreqmin(np,froot) - NtransPtoP(np,froot))

    if(casapool%nplant(np,leaf)/(casapool%cplant(np,leaf)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),leaf)) then
       Nreqmax(np,leaf) = 0.0
       Nreqmin(np,leaf) =0.0
    endif
    if(casapool%nplant(np,wood)/(casapool%cplant(np,wood)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),wood)) then
       Nreqmax(np,wood) = 0.0
       Nreqmin(np,wood) =0.0
    endif
    if(casapool%nplant(np,froot)/(casapool%cplant(np,froot)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),froot)) then
       Nreqmax(np,froot) = 0.0
       Nreqmin(np,froot) =0.0
    endif

  ENDIF
  ENDDO

END SUBROUTINE casa_Nrequire

