SUBROUTINE casa_rplant(veg,casabiome,casapool,casaflux,casamet)
! maintenance respiration of woody tisse and fineroots 
! see Sitch et al. (2003), GCB, reqn (23)

  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (climate_type),            INTENT(IN) :: climate
  INTEGER :: npt, ivt

  real(r_2), dimension(mp)        :: Ygrow        ! growth efficiency Q.Zhang 22/02/2011
  real(r_2), dimension(mp,mplant) :: ratioPNplant ! Q.Zhang 22/02/2011
  real(r_2), dimension(mp)        :: delcrmleaf, delcrmwood,delcrmfroot    ! reduction in wood and root respiration when NPP <0.0   
  real(r_2), dimension(mp)        :: resp_coeff_root, resp_coeff_sapwood, resp_coeff
  real,  dimension(mp)        :: nleaf, pleaf, vcmaxmax

  resp_coeff = 1
  resp_coeff_root = 1
  resp_coeff_sapwood = 1
  ratioPNplant = 0.0
  Ygrow        = 0.0

  WHERE(casapool%Nplant>0.0)
    ratioPNplant = casapool%Pplant/(casapool%Nplant+ 1.0e-10)
  ENDWHERE

  Ygrow(:) = 0.65+0.2*ratioPNplant(:,leaf)/(ratioPNplant(:,leaf)+1.0/15.0)

  casaflux%crmplant(:,wood) = 0.0
  casaflux%crmplant(:,froot) = 0.0
  delcrmleaf   = 0.0
  delcrmwood   = 0.0
  delcrmfroot  = 0.0
  casaflux%crgplant = 0.0
  casaflux%clabloss = 0.0

  WHERE(casamet%iveg2/=icewater) 
    WHERE(casamet%tairk >250.0) 
      WHERE(casapool%cplant(:,wood)>1.0e-6)
        casaflux%crmplant(:,wood)  = 
             casabiome%rmplant(veg%iveg(:),wood) &
                                   * casapool%nplant(:,wood)             &
                                   * exp(308.56*(1.0/56.02-1.0           &
                                   / (casamet%tairk(:)+46.02-tkzeroc)))
      ENDWHERE
      casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                            * max(0.0,casapool%Clabile(:))      &
                            * exp(308.56*(1.0/56.02-1.0         &
                            / (casamet%tairk(:)+46.02-tkzeroc)))
    ENDWHERE
    WHERE(casamet%tsoilavg >250.0.and.casapool%cplant(:,froot)>1.0e-6) 
      casaflux%crmplant(:,froot) = casabiome%rmplant(veg%iveg(:),froot) &
                                 * casapool%nplant(:,froot)             &
                                 * exp(308.56*(1.0/56.02-1.0            &
                                 / (casamet%tsoilavg(:)+46.02-tkzeroc)))
    ENDWHERE
!    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + casaflux%clabloss(:)

    WHERE((casaflux%Cgpp-SUM(casaflux%crmplant,2))>0.0)
    !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
    ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
      casaflux%crgplant(:)  = (1.0-Ygrow(:))* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))

    ELSEWHERE
      casaflux%crgplant(:) = 0.0
    ENDWHERE

!    casaflux%Cnpp(:) = MAX(0.0,(casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2)&
!                     - casaflux%crgplant(:))) 
    ! changes made by yp wang 5 april 2013
    ! to reduce wood and root resp when NPP < 0.0
    casaflux%Cnpp(:) = casaflux%Cgpp(:) - SUM(casaflux%crmplant(:,:),2) &
                     - casaflux%crgplant(:)
    WHERE(casaflux%Cnpp < 0.0)
! change made here by ypw on 11-7-2016 to include leaf maintenance respiration    
      delcrmleaf(:)  = casaflux%Cnpp(:) * casaflux%crmplant(:,leaf) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))
      delcrmwood(:)  = casaflux%Cnpp(:) * casaflux%crmplant(:,wood) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))
      delcrmfroot(:) = casaflux%Cnpp(:) * casaflux%crmplant(:,froot) &
                     / max(0.01,(casaflux%crmplant(:,leaf)+casaflux%crmplant(:,wood) &
                               + casaflux%crmplant(:,froot)))

      casaflux%crmplant(:,leaf)  = casaflux%crmplant(:,leaf)  + delcrmleaf(:)
      casaflux%crmplant(:,wood)  = casaflux%crmplant(:,wood)  + delcrmwood(:)
      casaflux%crmplant(:,froot) = casaflux%crmplant(:,froot) + delcrmfroot(:)
  !    casaflux%Cnpp(:) = casaflux%Cnpp(:) -delcrmwood(:)-delcrmfroot(:)
      casaflux%crgplant(:) = 0.0
    ENDWHERE
  ENDWHERE
  casaflux%Cnpp(:) = casaflux%Cgpp(:) - SUM(casaflux%crmplant(:,:),2) &
                   - casaflux%crgplant(:)
!  npt = 19360 
!  write(67,671) npt,veg%iveg(npt),casamet%tairk(npt),casamet%lat(npt),casamet%lon(npt), &
!                casapool%cplant(npt,:),casapool%nplant(npt,:), &
!                casaflux%Cgpp(npt), casaflux%Cnpp(npt), &
!                casaflux%crmplant(npt,:),casaflux%crgplant(npt),Ygrow(npt),&
!                ratioPNplant(npt,leaf),casapool%pplant(npt,leaf)

671 format('calling r plant',i6,1x,i3,1x,100(f10.4,2x))
END SUBROUTINE casa_rplant


