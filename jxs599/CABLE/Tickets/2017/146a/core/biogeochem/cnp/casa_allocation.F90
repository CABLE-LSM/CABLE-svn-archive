SUBROUTINE casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen,LALLOC)
! compute fraction of net photosynthate allocated to leaf, wood and froot
!
! inputs
!   moistavg(mp)           as an argument (volume fraction)
!   tsoilavg(mp)           as an argument (K)
!   btran(mp)              as an argument (dimensionless)
! outputs:
!   fracCalloc(mp,mplant1)
!
! modified Piere's alocation scheme
! input: leaf stage
!        leaf area

  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (phen_variable),       INTENT(INOUT) :: phen
  INTEGER , INTENT(IN) :: LALLOC
  ! local variables
  INTEGER :: npt,ns,is,iv
  REAL(r_2), DIMENSION(mp,mplant) :: fracCallocx
  REAL(r_2), DIMENSION(mp,mplant) :: delc
  REAL(r_2), DIMENSION(mp)        :: ctotal
  REAL(r_2), DIMENSION(mp)        :: xLalloc,xwsalloc,xTalloc
  REAL(r_2), DIMENSION(mp)        :: xWorNalloc,xNalloc,xWalloc
  REAL(r_2), DIMENSION(mp)        :: totfracCalloc
  REAL(r_2), DIMENSION(mp)        :: newLAI
  logical :: Ticket146 = .false.
  ! initlization
  casaflux%fracCalloc  = 0.0
  !Ticket146:VH comments out init
  !casaflux%fracClabile = 0.0
  fracCallocx = 0.0
  newLAI = 0.0
  SELECT CASE (LALLOC)

  CASE(2)   !
    ! calculate the allocation coefficients
    call casa_wolf(veg,casabiome,casaflux,casapool,casamet)

  CASE(1)   ! dynamic allocation
    WHERE(casamet%iveg2/=icewater)
      xLalloc(:) = min(1.0,max(0.0,exp(-0.5*casamet%glai(:))))   ! L limiting
      ! Pseudo-nutrient limitation calculation
      WHERE(casamet%tsoilavg > 0.0)
        xwsalloc(:) = min( max(casamet%moistavg(:)-soil%swilt(:),0.0) &
                         /(soil%sfc(:)-soil%swilt(:)), 1.0 )
      ELSE WHERE
        xwsalloc(:) = 0.01
      END WHERE
      xTalloc(:)    = min(1.0,max(0.0,Q10alloc** &
                      ((casamet%tsoilavg(:)-TkzeroC-30.0)/10.0) )) !T limiting
      xNalloc(:)    = min(1.0,max(0.0,xwsalloc(:)*xTalloc(:)))     !N limiting
      xWalloc(:)    = min(1.0,max(0.0,casamet%btran(:)))           !W limiting
      xWorNalloc(:) = min(xWalloc(:),xNalloc(:))
      WHERE(casamet%lnonwood==0)
        casaflux%fracCalloc(:,FROOT) = R0 * 3.0 * xLalloc(:) &
                                     / (xLalloc(:)+ 2.0*xWorNalloc(:))
        casaflux%fracCalloc(:,WOOD)  = S0 * 3.0 * xWorNalloc(:) &
                                     / (2.0*xLalloc(:)+ xWorNalloc(:))
        casaflux%fracCalloc(:,LEAF)  = 1.0 - casaflux%fracCalloc(:,FROOT) &
                                     - casaflux%fracCalloc(:,WOOD)
      ELSE WHERE
        casaflux%fracCalloc(:,FROOT) = R0 * 3.0 * xLalloc(:) &
                                     / (xLalloc(:)+2.0*xWorNalloc(:))
        casaflux%fracCalloc(:,WOOD)  = 0.0
        casaflux%fracCalloc(:,LEAF)  = 1.0 - casaflux%fracCalloc(:,FROOT)
      END WHERE
    END WHERE
  CASE (0)   ! fixed allocation
    casaflux%fracCalloc(:,:) = casabiome%fracnpptop(veg%iveg(:),:)
!Ticket146: VH implements case(3)
  CASE (3) ! leaf:wood allocation set to maintain LA:SA ratio
     ! below target value of 4000, where phen%phase = 1 or 2 
     !(requires casaflux%sapwood_area, which is inherited from the 
     ! POP tree demography module. (Ticket #61)
    WHERE(casamet%lnonwood==0)
        casaflux%fracCalloc(:,FROOT) =  casabiome%fracnpptop(veg%iveg(:),FROOT)
        casaflux%fracCalloc(:,WOOD) = 0.01
        casaflux%fracCalloc(:,LEAF) = 1.0 - casaflux%fracCalloc(:,FROOT) - &
             casaflux%fracCalloc(:,WOOD)
        newLAI =casamet%glai + (casaflux%fracCalloc(:,LEAF) *casaflux%cnpp- &
             casaflux%kplant(:,leaf) *casapool%cplant(:,LEAF) )*casabiome%sla(veg%iveg(:))
        where (casaflux%sapwood_area.gt.1.e-6 .and. newLAI.gt.(4000.*casaflux%sapwood_area) &
             .and. casaflux%cnpp.gt.0.0)

           casaflux%fracCalloc(:,LEAF) = ((4000.*casaflux%sapwood_area - casamet%glai)/ &
                casabiome%sla(veg%iveg(:)) &
             + casaflux%kplant(:,leaf) *casapool%cplant(:,LEAF)  )/casaflux%cnpp

           casaflux%fracCalloc(:,LEAF) = max(0.0,  casaflux%fracCalloc(:,LEAF) )
           casaflux%fracCalloc(:,LEAF) = min(1.0 - casaflux%fracCalloc(:,FROOT) - &
                casaflux%fracCalloc(:,WOOD) ,&
             casaflux%fracCalloc(:,LEAF) )

           casaflux%fracCalloc(:,WOOD) = 1.0 -  casaflux%fracCalloc(:,FROOT) - &
                casaflux%fracCalloc(:,LEAF)
        end where


     ELSEWHERE

        casaflux%fracCalloc(:,FROOT) =  casabiome%fracnpptop(veg%iveg(:),FROOT)
        casaflux%fracCalloc(:,WOOD) = 0.0
        casaflux%fracCalloc(:,LEAF) =  casabiome%fracnpptop(veg%iveg(:),LEAF)

     ENDWHERE

  END SELECT

  ! during leaf growth phase 0 or 3, no carbon is allocated to leaf,
  ! during maximal leaf growth phase, all C is allocated to leaf
  ! during steady growth period, C allocation is estimated in such
  ! a way that approach the allometric relationship
  ! the relationships are:(all pools in g C/m2)
  ! for forests
  !   fineroot/totalC C=0.3192-0.0485*(totalC)^0.1755, see mokany et al. (2003)
  !   fineroot = ratiofrootleaf*cleaf
  ! for grassland
  !   root=ratiofinerootleaf*cleaf

if( Ticket46 ) then

  WHERE(casamet%iveg2/=icewater) 
    WHERE(phen%phase==0) 
      casaflux%fracCalloc(:,leaf)  = 0.0
      casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                                   /(casaflux%fracCalloc(:,froot) &
                                     +casaflux%fracCalloc(:,wood))
      casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
    END WHERE 

    WHERE(phen%phase==1)          
      casaflux%fracCalloc(:,leaf)  = 0.8
      WHERE(casamet%lnonwood==0)  !woodland or forest
        casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
        casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
      ELSEWHERE !grassland
        casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
      ENDWHERE
    END WHERE

    WHERE(phen%phase==3) 
!      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
      casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,wood) 
      casaflux%fracCalloc(:,leaf)  = 0.0
    ENDWHERE

  ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
  ! Q.Zhang 17/03/2011
    WHERE(casamet%glai(:)>=0.95*casabiome%glaimax(veg%iveg(:)))
      casaflux%fracCalloc(:,leaf)  = 0.0
      casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                                   /(casaflux%fracCalloc(:,froot) &
                                     +casaflux%fracCalloc(:,wood))
      casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
    ENDWHERE

    WHERE(casamet%glai(:)<1.05*casabiome%glaimin(veg%iveg(:)))
      casaflux%fracCalloc(:,leaf)  = 0.8
      WHERE(casamet%lnonwood==0)  !woodland or forest
        casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
        casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
      ELSEWHERE !grassland
        casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
      ENDWHERE
    ENDWHERE

    !! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
    WHERE(casaflux%Cnpp<0.0.and.sum(casapool%cplant,2)>0.0)
 !      casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/sum(casaflux%Crmplant,2)
 !      casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)/sum(casaflux%Crmplant,2)
 !      casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/sum(casaflux%Crmplant,2)
       casaflux%fracCalloc(:,leaf)  = sum(casaflux%Crmplant,2) * casapool%cplant(:,leaf)/sum(casapool%cplant,2)
       casaflux%fracCalloc(:,wood)  = sum(casaflux%Crmplant,2) * casapool%cplant(:,wood)/sum(casapool%cplant,2)
       casaflux%fracCalloc(:,froot) = sum(casaflux%Crmplant,2) * casapool%cplant(:,froot)/sum(casapool%cplant,2)
    ENDWHERE

  ENDWHERE

else (Ticket46)

! vh edit to avoid overwriting CASE(3) for woody veg
!! vh_js !!
  IF (LALLOC.ne.(3)) THEN

     WHERE(casamet%iveg2/=icewater)
        WHERE(phen%phase==0)
           casaflux%fracCalloc(:,leaf)  = 0.0
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
        END WHERE

        WHERE(phen%phase==1)
           casaflux%fracCalloc(:,leaf)  = 0.8
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
              casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
           ELSEWHERE !grassland
              casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
           ENDWHERE
        END WHERE

        WHERE(phen%phase==3)
           !      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
           casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,wood)
           casaflux%fracCalloc(:,leaf)    = 0.0
        ENDWHERE


        ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
        ! Q.Zhang 17/03/2011
        WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
           casaflux%fracCalloc(:,leaf)  = 0.0
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
        ENDWHERE

        ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
        WHERE(casaflux%Cnpp<0.0)
           casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/sum(casaflux%Crmplant,2)
           casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)/sum(casaflux%Crmplant,2)
           casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/sum(casaflux%Crmplant,2)
        ENDWHERE

        !! vh_js !!
        !! as long as biomass is positive, adjust allocation to be
        !! proportional to stock when NPP -ve   (Ticket#108)
        WHERE(casaflux%Cnpp<0.0 .and. sum(casapool%Cplant,2)>0  )
           casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)/sum(casapool%Cplant,2)
           casaflux%fracCalloc(:,wood)  = casapool%Cplant(:,wood)/sum(casapool%Cplant,2)
           casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot)/sum(casapool%Cplant,2)
        ENDWHERE
     ENDWHERE

  ELSE
     WHERE(casamet%iveg2/=icewater)
        WHERE(phen%phase==0)
           casaflux%fracCalloc(:,leaf)  = 0.0
           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                /(casaflux%fracCalloc(:,froot) &
                +casaflux%fracCalloc(:,wood))
           WHERE (casamet%lnonwood==0)
              casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
           ELSEWHERE
              casaflux%fracCalloc(:,wood) = 0.0
           ENDWHERE
        END WHERE

        WHERE(phen%phase==1.and.casamet%lnonwood==1)

           casaflux%fracCalloc(:,leaf)  = 0.8
           casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
           casaflux%fracCalloc(:,wood) = 0.0
        ENDWHERE

        WHERE(phen%phase==3)
           !      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
           casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,wood)
           casaflux%fracCalloc(:,leaf)  = 0.0
        ENDWHERE

!! vh !! don't require this fix for LALLOC = 3 (POP allocation scheme)
!! Thiss fix can lead to over-allocation to roots, in turn bumping up N-uptake
!! , leading to decline in mineral nitrogen availability and spikes in fracCalloc,
!! causing spikes in tree mortality and lack of model convergence in productive
!! regions where LAI is hitting LAImax.
!!$        ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
!!$        ! Q.Zhang 17/03/2011
!!$        WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
!!$           casaflux%fracCalloc(:,leaf)  = 0.0
!!$           casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
!!$                /(casaflux%fracCalloc(:,froot) &
!!$                +casaflux%fracCalloc(:,wood))
!!$           WHERE (casamet%lnonwood==0)
!!$              casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
!!$           ELSEWHERE
!!$              casaflux%fracCalloc(:,wood) = 0.0
!!$           ENDWHERE
!!$        ENDWHERE

        WHERE(casamet%glai(:)<casabiome%glaimin(veg%iveg(:)))
           casaflux%fracCalloc(:,leaf)  = 0.8
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
              casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
           ELSEWHERE !grassland
              casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
              casaflux%fracCalloc(:,wood) = 0.0
           ENDWHERE
        ENDWHERE
        ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
        WHERE(casaflux%Cnpp<0.0)
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/sum(casaflux%Crmplant,2)
              casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)/sum(casaflux%Crmplant,2)
              casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/sum(casaflux%Crmplant,2)
           ELSEWHERE
              casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/sum(casaflux%Crmplant,2)
              casaflux%fracCalloc(:,wood)  = 0.0
              casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/sum(casaflux%Crmplant,2)
           ENDWHERE
        ENDWHERE
        
        !! vh_js !!  Ticket#108 
        WHERE(casaflux%Cnpp<0.0 .and. sum(casapool%Cplant,2)>0  )
           WHERE(casamet%lnonwood==0)  !woodland or forest
              casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)/sum(casapool%Cplant,2)
              casaflux%fracCalloc(:,wood)  = casapool%Cplant(:,wood)/sum(casapool%Cplant,2)
              casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot)/sum(casapool%Cplant,2)
           ELSEWHERE
              casaflux%fracCalloc(:,leaf)  = casapool%Cplant(:,leaf)/sum(casapool%Cplant,2)
              casaflux%fracCalloc(:,wood)  = 0.0
              casaflux%fracCalloc(:,froot) = casapool%Cplant(:,froot)/sum(casapool%Cplant,2)
           ENDWHERE
        ENDWHERE
        
        
     ENDWHERE
  !   write(*,*) 'alloc2',  casaflux%fracCalloc(1,2), casaflux%Cnpp(1), casapool%Cplant(1,:), &
  !        casamet%lnonwood(1)
!if (ANY(casapool%Cplant(1,:).NE.casapool%Cplant(1,:))) then
!write(*,*) 'cplant', casapool%Cplant(1,:)
!stop
!endif
  ENDIF ! LALLOC=3
ENDIF ! TIcket46 
  



  ! normalization the allocation fraction to ensure they sum up to 1
  totfracCalloc(:) = sum(casaflux%fracCalloc(:,:),2)
  casaflux%fracCalloc(:,leaf) = casaflux%fracCalloc(:,leaf)/totfracCalloc(:)
  casaflux%fracCalloc(:,wood) = casaflux%fracCalloc(:,wood)/totfracCalloc(:)
  casaflux%fracCalloc(:,froot) = casaflux%fracCalloc(:,froot)/totfracCalloc(:)

END SUBROUTINE casa_allocation

