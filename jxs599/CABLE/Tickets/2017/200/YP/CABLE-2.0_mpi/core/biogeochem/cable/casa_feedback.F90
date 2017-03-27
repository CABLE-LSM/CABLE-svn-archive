 SUBROUTINE casa_feedback(ktau,veg,casabiome,casapool,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE casa_cnp_module, ONLY: vcmax_np
  USE cable_common_module,  ONLY:  CABLE_USER
  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(IN) :: casapool
  TYPE (casa_met),            INTENT(IN) :: casamet

  integer np,ivt
  real, dimension(mp)  :: ncleafx,npleafx, pleafx, nleafx ! local variables
  real, dimension(17)                   ::  xnslope
  data xnslope/0.80,1.00,2.00,1.00,1.00,1.00,0.50,1.00,0.34,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/

  ! first initialize
  ncleafx(:) = casabiome%ratioNCplantmax(veg%iveg(:),leaf)
  if(Ticket200) then
    npleafx(:) = casabiome%ratioNPplantmin(veg%iveg(:),leaf)
  else 
    npleafx = 14.2
  endif
  
  DO np=1,mp
    ivt=veg%iveg(np)
    IF (casamet%iveg2(np)/=icewater &
        .AND. casamet%glai(np)>casabiome%glaimin(ivt)  &
        .AND. casapool%cplant(np,leaf)>0.0) THEN
      ncleafx(np) = MIN(casabiome%ratioNCplantmax(ivt,leaf), &
                        MAX(casabiome%ratioNCplantmin(ivt,leaf), &
                            casapool%nplant(np,leaf)/casapool%cplant(np,leaf)))
      IF (icycle>2 .AND. casapool%pplant(np,leaf)>0.0) THEN
        npleafx(np) = MIN(30.0,MAX(8.0,real(casapool%nplant(np,leaf) &
                /casapool%pplant(np,leaf))))
      ENDIF
    ENDIF

    IF (TRIM(cable_user%vcmax).eq.'standard') then
    
       IF (casamet%glai(np) > casabiome%glaimin(ivt)) THEN
          IF (ivt/=2) THEN
            veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                        + casabiome%nslope(ivt)     &
                        * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
      ELSE
        IF (casapool%nplant(np,leaf)>0.0.AND.casapool%pplant(np,leaf)>0.0) THEN
          veg%vcmax(np) = ( casabiome%nintercept(ivt)  &
                          + casabiome%nslope(ivt)*(0.4+9.0/npleafx(np)) &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ELSE
          veg%vcmax(np) = ( casabiome%nintercept(ivt) &
                          + casabiome%nslope(ivt)     &
                          * ncleafx(np)/casabiome%sla(ivt) ) * 1.0e-6
        ENDIF
      ENDIF
    ENDIF

!    veg%vcmax(np) = ( nintercept(ivt)  &
!                  + nslope(ivt)*(0.4+8.5/npleafx(np)) &
!                  * ncleafx(np)/casabiome%sla(ivt))*(1.0e-6)
!    veg%vcmax(np) =veg%vcmax(np)* xnslope(ivt)

!write(*,991) np, ivt,veg%vlai(np),veg%vcmax(np)*1.0e6
!write(*,891) np,ivt,casapool%cplant(np,leaf),casapool%nplant(np,leaf),casapool%pplant(np,leaf)
!891 format(2(i6),3(f9.3,2x))
  ENDDO

  veg%ejmax = 2.0 * veg%vcmax
!991 format(i6,2x,i4,2x,2(f9.3,2x))
 END SUBROUTINE casa_feedback


