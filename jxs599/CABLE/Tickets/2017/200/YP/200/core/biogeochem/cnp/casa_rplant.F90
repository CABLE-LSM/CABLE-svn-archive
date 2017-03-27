SUBROUTINE casa_rplant(veg,casabiome,casapool,casaflux,casamet,climate)
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
  
  if (cable_user%CALL_climate) then
  ! coefficients required to implement T-acclimation of autotrophic respiration (Ticket # 110)
  ! adapted from Atkin et al., New Phyt., 2015)
     DO npt = 1, mp
        ivt=veg%iveg(npt)
        ! max leaf N in g N m-2 leaf
        nleaf(npt) =  casabiome%ratioNCplantmax(ivt,leaf)/casabiome%sla(ivt) 
        ! max leaf P in g P m-2 leaf
        pleaf(npt) = casabiome%ratioPcplantmax(ivt,leaf)/casabiome%sla(ivt)  
        if (ivt .EQ. 7) then
           ! special for C4 grass: set here to value from  parameter file
           vcmaxmax(npt) = 1.0e-5 
        else
           vcmaxmax(npt) = vcmax_np(nleaf(npt), pleaf(npt))
        endif
        if (veg%iveg(npt).eq.2 .or. veg%iveg(npt).eq. 4  ) then 
           ! broadleaf forest

           resp_coeff_root(npt) = (1.2818 * 1.e-6 *casapool%nplant(npt,froot)/ &
                vcmaxmax(npt)/0.0116   + &
                casapool%nplant(npt,froot)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )  

           resp_coeff_sapwood(npt) = (1.2818 * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      ) 


        elseif (veg%iveg(npt).eq.1 .or. veg%iveg(npt).eq. 3  ) then 
           ! needleleaf forest

           resp_coeff_root(npt) = (1.2877 * 1.e-6 *casapool%nplant(npt,froot) &
                /vcmaxmax(npt)/0.0116   + &
                casapool%nplant(npt,froot)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )  

           resp_coeff_sapwood(npt) = (1.2877 * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      ) 



        elseif (veg%iveg(npt).eq.6 .or. veg%iveg(npt).eq.8 .or. veg%iveg(npt).eq. 9  ) then 
           ! C3 grass, tundra, crop

           resp_coeff_root(npt) = (1.6737 * 1.e-6 *casapool%nplant(npt,froot)/ &
                vcmaxmax(npt)/0.0116   + &
                casapool%nplant(npt,froot)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )  

           resp_coeff_sapwood(npt) = (1.6737 * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      ) 
        else 
           ! shrubs and other (C4 grass and crop)
           resp_coeff_root(npt) = (1.5758 * 1.e-6 *casapool%nplant(npt,froot)/ &
                vcmaxmax(npt)/0.0116   + &
                casapool%nplant(npt,froot)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 * &
                casapool%nplant(npt,froot)/vcmaxmax(npt)/0.0116      )  

           resp_coeff_sapwood(npt) = (1.5758 * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116  + &
                casapool%nplant(npt,wood) * casaflux%frac_sapwood(npt)  - &
                0.0334* climate%qtemp_max_last_year(npt) * 1.e-6 *casapool%nplant(npt,wood) * &
                casaflux%frac_sapwood(npt)/vcmaxmax(npt)/0.0116      ) 
        endif
     ENDDO
     resp_coeff = 0.50
  ENDIF  ! end coefficients for acclimation of autotrophic respiration Ticket #110

!Ticket200 
IF (cable_user%CALL_climate) then 
   !  acclimation of autotrophic respiration Ticket #110
     WHERE(casamet%iveg2/=icewater)
        WHERE(casamet%tairk >250.0)
           WHERE(casapool%cplant(:,wood)>1.0e-6)
              casaflux%crmplant(:,wood)  =  resp_coeff  * &
              resp_coeff_sapwood * &
                   casabiome%rmplant(veg%iveg(:),wood) &
                   * exp(308.56*(1.0/56.02-1.0           &
                   / (casamet%tairk(:)+46.02-tkzeroc)))

           ENDWHERE
           !vh! prevent floating underflow with this mask
           WHERE (casapool%Clabile(:).gt.1.e-8) &
              casaflux%clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                   * max(0.0,casapool%Clabile(:))      &
                   * exp(308.56*(1.0/56.02-1.0         &
                   / (casamet%tairk(:)+46.02-tkzeroc)))
           

        ENDWHERE

        WHERE(casamet%tsoilavg >250.0.and.casapool%cplant(:,froot)>1.0e-6)

           casaflux%crmplant(:,froot) =  resp_coeff * resp_coeff_root * &
                casabiome%rmplant(veg%iveg(:),froot) &
                * exp(308.56*(1.0/56.02-1.0            &
                / (casamet%tsoilavg(:)+46.02-tkzeroc)))

        ENDWHERE

        WHERE((casaflux%Cgpp-SUM(casaflux%crmplant,2))>0.0)
           !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
           ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
           casaflux%crgplant(:)  = (1.0-Ygrow(:))* &
                max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))

        ELSEWHERE
           casaflux%crgplant(:) = 0.0
        ENDWHERE
     ENDWHERE

     Casaflux%cnpp(:) = casaflux%Cgpp(:)-Sum(casaflux%crmplant(:,:),2) - casaflux%crgplant(:)

  ELSE


     WHERE(casamet%iveg2/=icewater)
        WHERE(casamet%tairk >250.0)
           WHERE(casapool%cplant(:,wood)>1.0e-6)
             !Ticket200
              casaflux%crmplant(:,wood)  =  resp_coeff * 
                   casaflux%frac_sapwood(:) * &
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

           casaflux%crmplant(:,froot) =  resp_coeff * casabiome%rmplant(veg%iveg(:),froot) &
                * casapool%nplant(:,froot)             &
                * exp(308.56*(1.0/56.02-1.0            &
                / (casamet%tsoilavg(:)+46.02-tkzeroc)))

        ENDWHERE
        !    casaflux%crmplant(:,leaf) = casaflux%crmplant(:,leaf) + casaflux%clabloss(:)

        WHERE((casaflux%Cgpp-SUM(casaflux%crmplant,2))>0.0)
           !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))
           ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
           casaflux%crgplant(:)  = (1.0-Ygrow(:))* &
                max(0.0,casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2))

        ELSEWHERE
           casaflux%crgplant(:) = 0.0
        ENDWHERE


        !casaflux%Cnpp(:) = MAX(0.0,(casaflux%Cgpp(:)-SUM(casaflux%crmplant(:,:),2) &
        !                 - casaflux%crgplant(:)))
        ! changes made by yp wang 5 april 2013
        Casaflux%cnpp(:) = casaflux%Cgpp(:)-Sum(casaflux%crmplant(:,:),2) - casaflux%crgplant(:)


     ENDWHERE

  ENDIF

END SUBROUTINE casa_rplant


