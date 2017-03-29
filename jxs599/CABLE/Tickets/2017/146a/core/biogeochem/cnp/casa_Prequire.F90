SUBROUTINE casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                         casabiome,casapool,casaflux,casamet)
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xpCnpp
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Preqmax, Preqmin
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: PtransPtoP
  TYPE (veg_parameter_type),             INTENT(INOUT) :: veg
  TYPE (casa_biome),                     INTENT(INOUT) :: casabiome
  TYPE (casa_pool),                      INTENT(INOUT) :: casapool
  TYPE (casa_flux),                      INTENT(INOUT) :: casaflux
  TYPE (casa_met),                       INTENT(INOUT) :: casamet

  ! local variables
  INTEGER :: nland,np,ip

  Preqmin(:,:)       = 0.0
  Preqmax(:,:)       = 0.0
  PtransPtoP(:,:)    = 0.0
  do np=1,mp
  IF(casamet%iveg2(np)/=icewater) then
    Preqmax(np,leaf) = xpCnpp(np)* casaflux%fracCalloc(np,leaf) &
                    * (casapool%Nplant(np,leaf)/(casapool%Cplant(np,leaf)+1.0e-10))/casabiome%ratioNPplantmin(veg%iveg(np),leaf)
    Preqmax(np,wood) = xpCnpp(np)* casaflux%fracCalloc(np,wood) &
                    * (casapool%Nplant(np,wood)/(casapool%Cplant(np,wood)+1.0e-10))/casabiome%ratioNPplantmin(veg%iveg(np),wood)
    Preqmax(np,froot) = xpCnpp(np)* casaflux%fracCalloc(np,froot) &
                    * (casapool%Nplant(np,froot)/(casapool%Cplant(np,froot)+1.0e-10))/casabiome%ratioNPplantmin(veg%iveg(np),froot)

    Preqmin(np,leaf) = xpCnpp(np) * casaflux%fracCalloc(np,leaf) &
                    * (casapool%Nplant(np,leaf)/(casapool%Cplant(np,leaf)+1.0e-10))/casabiome%ratioNPplantmax(veg%iveg(np),leaf)
    Preqmin(np,wood) = xpCnpp(np) * casaflux%fracCalloc(np,wood) &
                    * (casapool%Nplant(np,wood)/(casapool%Cplant(np,wood)+1.0e-10))/casabiome%ratioNPplantmax(veg%iveg(np),wood)
    Preqmin(np,froot) = xpCnpp(np) * casaflux%fracCalloc(np,froot) &
                    * (casapool%Nplant(np,froot)/(casapool%Cplant(np,froot)+1.0e-10))/casabiome%ratioNPplantmax(veg%iveg(np),froot)

    PtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Pplant(np,leaf) &
                       * (1.0-casabiome%ftransPPtoL(veg%iveg(np),leaf))
    PtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Pplant(np,wood) &
                       * (1.0-casabiome%ftransPPtoL(veg%iveg(np),wood))
    PtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Pplant(np,froot) &
                       * (1.0-casabiome%ftransPPtoL(veg%iveg(np),froot))

    Preqmax(np,leaf)    = max(0.0,Preqmax(np,leaf) - PtransPtoP(np,leaf))
    Preqmax(np,wood)    = max(0.0,Preqmax(np,wood) - PtransPtoP(np,wood))
    Preqmax(np,froot)    = max(0.0,Preqmax(np,froot) - PtransPtoP(np,froot))

    Preqmin(np,leaf)    = max(0.0,Preqmin(np,leaf) - PtransPtoP(np,leaf))
    Preqmin(np,wood)    = max(0.0,Preqmin(np,wood) - PtransPtoP(np,wood))
    Preqmin(np,froot)    = max(0.0,Preqmin(np,froot) - PtransPtoP(np,froot))


    if(casapool%pplant(np,leaf)/(casapool%nplant(np,leaf)+1.0e-10)> 1.0/casabiome%ratioNPplantmin(veg%iveg(np),leaf)) then
       Preqmax(np,leaf) = 0.0
       Preqmin(np,leaf) =0.0
    endif
    if(casapool%pplant(np,wood)/(casapool%nplant(np,wood)+1.0e-10)> 1.0/casabiome%ratioNPplantmin(veg%iveg(np),wood)) then
       Preqmax(np,wood) = 0.0
       Preqmin(np,wood) =0.0
    endif
    if(casapool%pplant(np,froot)/(casapool%nplant(np,froot)+1.0e-10)> 1.0/casabiome%ratioNPplantmin(veg%iveg(np),froot)) then
       Preqmax(np,froot) = 0.0
       Preqmin(np,froot) =0.0
    endif

  endif
  ENDDO

END SUBROUTINE casa_Prequire

