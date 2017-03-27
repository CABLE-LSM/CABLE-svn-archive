! changed by yp wang following Chris Lu 5/nov/2012
SUBROUTINE biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                      casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
                      cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                      nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                      pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
  USE cable_def_types_mod
  USE casadimension
  USE casa_cnp_module
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ktau
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: idoy
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables added by ypwang following Chris Lu 5/nov/2012

  real, dimension(mp), INTENT(OUT)   :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                                        nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                                        pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

  ! local variables
  REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
  REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
  REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
  INTEGER  npt,j

!  npt =33748

  xKNlimiting = 1.0
  call phenology(idoy,veg,phen)
  call avgsoil(veg,soil,casamet)
  call casa_rplant(veg,casabiome,casapool,casaflux,casamet)

!  print *, 'biogeochem1', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)
  call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen)

  call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                       casamet,phen)

!  print *, 'biogeochem2', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)

  call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                       casaflux,casamet)

!  print *, 'biogeochem3', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)

  call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

!  print *, 'biogeochem4', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)
!  print *, 'calling casa_xratesoil ???'

  call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)

!  print *, 'biogeochem5', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)

  call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casapool,casaflux,casamet)

!  print *, 'biogeochem6', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)

  IF (icycle>1) THEN
    call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)

 !   print *, 'biogeochem7', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)

    DO j=1,mlitter
      casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
    ENDDO
 !   print *, 'biogeochem8', casaflux%ksoil(npt,:),casapool%Nsoil(npt,:)
    call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
    IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
                                     casapool,casaflux,casamet)
  ENDIF 

  ! changed by ypwang following Chris Lu on 5/nov/2012
!   write(*,900) ktau,idoy,npt,casapool%cplant(npt,:),casapool%nplant(npt,:), casapool%pplant(npt,:), &
!               casaflux%cgpp(npt),casaflux%Cnpp(npt),casaflux%crmplant(npt,:),casaflux%Crgplant(npt), &
!               casaflux%nupland(npt),casaflux%pupland(npt),xkNlimiting(npt),xnplimit(npt),xNPuptake(npt)

  call casa_delplant(veg,casabiome,casapool,casaflux,casamet,                &
                         cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

  !  call casa_delplant(veg,casabiome,casapool,casaflux,casamet)

  call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)

  call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)
  ! modified by ypwang following Chris Lu on 5/nov/2012

  IF (icycle<3) then
      call casa_pdummy(casapool)
      IF (icycle<2) call casa_ndummy(casapool)
  ENDIF

  call casa_cnpbal(casapool,casaflux,casabal)
  call casa_cnpflux(casaflux,casapool,casabal)

!  ! for spinning up only
!  casapool%Nsoilmin = max(casapool%Nsoilmin,0.5)
  casapool%Psoillab = max(casapool%Psoillab,0.01)

!    write(*,901) ktau,idoy,npt,phen%phase(npt),casapool%cplant(npt,:),casapool%nplant(npt,:), casapool%pplant(npt,:), &
!               casaflux%cgpp(npt),casaflux%Cnpp(npt),casaflux%crmplant(npt,:),casaflux%Crgplant(npt), &
!               casaflux%fracCalloc(npt,:),casaflux%fracClabile(npt),               &
!               casapool%Nsoilmin(npt), casaflux%Nupland(npt),                       &
!               casapool%psoillab(npt), casaflux%Pupland(npt),                       &
!               casamet%glai(npt),casabiome%glaimin(veg%iveg(npt)),casabiome%glaimax(veg%iveg(npt))

901 format('after delplant: ',4(i6,1x),100(f9.2,1x))
900 format('before delplant: ',3(i6,2x),100(f8.4,2x))
END SUBROUTINE biogeochem


