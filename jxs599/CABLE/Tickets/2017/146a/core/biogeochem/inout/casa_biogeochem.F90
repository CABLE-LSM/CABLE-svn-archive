SUBROUTINE biogeochem(ktau,dels,idoY,LALLOC,veg,soil,casabiome,casapool,casaflux, &
     casamet,casabal,phen,POP,climate,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
     cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                      pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
!Ticket146:YP's vn below - consolidate arg list w call from _driver
!SUBROUTINE biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
!                      casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
!                      cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
!                      nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
!                      pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)                      
  USE cable_def_types_mod
  USE casadimension
  USE casa_cnp_module
  USE POP_TYPES,            ONLY: POP_TYPE
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ktau
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: idoy
  INTEGER, INTENT(IN)    :: LALLOC
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen
  TYPE(POP_TYPE),             INTENT(IN) :: POP
  TYPE(climate_TYPE),             INTENT(IN) :: climate

  ! local variables added by ypwang following Chris Lu 5/nov/2012

  real, dimension(mp), INTENT(OUT)   :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

  ! local variables
  REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
  REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
  REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
  INTEGER  npt,j
  REAL, ALLOCATABLE :: tmp(:)
  logical :: Ticket146 = .false.
  
  xKNlimiting = 1.0

  !Ticket146:NB YP does not call
 ! zero annual sums
  if (idoy==1) CALL casa_cnpflux(casaflux,casapool,casabal,.true.)

  !Ticket146:NB YP does not have condititonal switch
  IF (cable_user%PHENOLOGY_SWITCH.eq.'MODIS') THEN
     call phenology(idoy,veg,phen)
  ENDIF
  
  call avgsoil(veg,soil,casamet)
  
  !Ticket146:NB YP does not pass climate%
  if(.NOT. Ticket146 ) then
    call casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)
  else
    call casa_rplant(veg,casabiome,casapool,casaflux,casamet)
  endif
  
  !Ticket146:NB YP does not call POP
   IF (.NOT.cable_user%CALL_POP) THEN
      call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
   ENDIF

   call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
        casamet,phen)

   !Ticket146:NB YP does not pass phen%
   call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
        casaflux,casamet,phen)

   call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

   IF (cable_user%CALL_POP) THEN 

      call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)

      WHERE (pop%pop_grid(:)%cmass_sum_old.gt.0.001 .and. pop%pop_grid(:)%cmass_sum.gt.0.001 )
         
         casaflux%frac_sapwood(POP%Iwood) = POP%pop_grid(:)%csapwood_sum/ POP%pop_grid(:)%cmass_sum
         casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6)
         veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
          
         WHERE (pop%pop_grid(:)%LU ==2)

            casaflux%kplant(POP%Iwood,2) =  1.0 -  &
              (1.0-  max( min((POP%pop_grid(:)%stress_mortality + &
              POP%pop_grid(:)%crowding_mortality+ &
              + POP%pop_grid(:)%fire_mortality ) &
              /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth) + &
              1.0/veg%disturbance_interval(POP%Iwood,1), 0.99), 0.0))**(1.0/365.0)

         ELSEWHERE !(pop%pop_grid(:)%LU ==2)
         
            casaflux%kplant(POP%Iwood,2) =  1.0 -  &
              (1.0-  max( min((POP%pop_grid(:)%stress_mortality + &
              POP%pop_grid(:)%crowding_mortality+ &
              + POP%pop_grid(:)%fire_mortality+POP%pop_grid(:)%cat_mortality  ) &
              /(POP%pop_grid(:)%cmass_sum+POP%pop_grid(:)%growth), 0.99), 0.0))**(1.0/365.0)

         ENDWHERE !(pop%pop_grid(:)%LU ==2)

         veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
      
      ELSEWHERE
      
         casaflux%frac_sapwood(POP%Iwood) = 1.0
         casaflux%sapwood_area(POP%Iwood) = max(POP%pop_grid(:)%sapwood_area/10000., 1e-6)
         casaflux%kplant(POP%Iwood,2) = 0.0
         veg%hc(POP%Iwood) = POP%pop_grid(:)%height_max
      
      ENDWHERE
    
   ENDIF ! (cable_user%CALL_POP) 

  call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
  
  if (.NOT. Ticket146) then
    call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
  else
    call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casapool,casaflux,casamet)
  endif

  IF (icycle>1) THEN
    call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)

    DO j=1,mlitter
      casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
    ENDDO

    call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
    IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
                                     casapool,casaflux,casamet)
  ENDIF

  call casa_delplant(veg,casabiome,casapool,casaflux,casamet,                &
       cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
       nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
       pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

  if (.NOT. Ticket146) then
    casaflux%Cplant_turnover_disturbance = 0
    casaflux%Cplant_turnover_crowding = 0
    casaflux%Cplant_turnover_resource_limitation = 0
  endif
  
  if (cable_user%CALL_POP) THEN
     if (.not.allocated(tmp)) allocate(tmp(size(POP%pop_grid)))
     tmp = (POP%pop_grid(:)%stress_mortality + POP%pop_grid(:)%crowding_mortality &
          +POP%pop_grid(:)%cat_mortality &
          + POP%pop_grid(:)%fire_mortality  )
     where (tmp.gt. 1.e-12)
        casaflux%Cplant_turnover_disturbance(POP%Iwood) =  &
             casaflux%Cplant_turnover(POP%Iwood,2)*(POP%pop_grid(:)%cat_mortality &
             + POP%pop_grid(:)%fire_mortality  )/tmp
        casaflux%Cplant_turnover_crowding(POP%Iwood) =  &
             casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%crowding_mortality/tmp
        casaflux%Cplant_turnover_resource_limitation(POP%Iwood) = &
             casaflux%Cplant_turnover(POP%Iwood,2)*POP%pop_grid(:)%stress_mortality/tmp
     endwhere
  endif

  call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)

  if (.NOT. Ticket146) then
    call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet, LALLOC)
  else
    call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)
  endif
  
  !Ticket146  
  !ndummy must be before pdummy
  IF (icycle<3) then
      IF (icycle<2) call casa_ndummy(casapool)
      call casa_pdummy(casapool)
  ENDIF

  call casa_cnpbal(casapool,casaflux,casabal)

  call casa_cnpflux(casaflux,casapool,casabal,.false.)

  if (Ticket146) &
    casapool%Psoillab = max(casapool%Psoillab,0.01)

END SUBROUTINE biogeochem

