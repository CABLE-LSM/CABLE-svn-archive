SUBROUTINE casa_wolf(veg,casabiome,casaflux,casapool,casamet)
   ! carbon allocation based on WFB2011
   ! Wolf,Field and Berry, 2011. Ecological Applications, p1546-1556
   ! Wolf et al. 2011. Global Biogeochemical Cycles, 25, GB3015,
   ! doi:10.1019/2010GB003917
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(IN) :: casabiome
  TYPE (casa_met),            INTENT(IN) :: casamet
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux

   real, parameter :: wolf_alpha1=6.22
   real, parameter :: wolf_beta=-1.33
   real, parameter :: wolf_c1=-wolf_alpha1/(1+wolf_beta)
   real, parameter :: wolf_c2=1.0/(1.0+wolf_beta)

   real  totleaf,totwood,totcroot,totfroot,totnpp
   real  fracleaf,fracwood,fraccroot,fracfroot
   !
   ! local variables
   integer   npt
   real(r_2), dimension(mp)  ::  totbmdm,ntree,nppdm
   real(r_2), dimension(mp)  ::  gleaf,gwood,gcroot,gfroot,gtot
   !
   ! input
   !  totleaf, totwood, totcroot, totfroot :    g C m-2
   !  totnpp:                                   g C m-2 d-1
   ! output
   !  fracleaf,fracwood, fraccroot, fracfroot:  fractions
   !

    do npt=1,mp
       IF(casamet%iveg2(npt)==3.and.casaflux%cnpp(npt)>0.0001) THEN  !forest types
          totbmdm(npt) = sum(casapool%cplant(npt,:)) *10.0 / fracCbiomass      !10.0 for convert gc/m2 to kg/ha
          totbmdm(npt) = max(30000.0, totbmdm(npt))
          ! calculate tree stocking density
           ntree(npt) = 10**(wolf_c1+wolf_c2*log10(totbmdm(npt)))   ! tree ha-1, based on eqn (4) of WFB2011
           ntree(npt) = min(200000.0,ntree(npt))
           ! changed by ypw 23/april/2012 to avoid negative npp
           nppdm(npt)  = (abs(casaflux%cnpp(npt)) *365.0*0.001/fracCbiomass)/(0.0001*ntree(npt))  ! in kg dm tree-1 yr-1

           gleaf(npt)  = 0.156*(nppdm(npt)**1.106)     ! Figure 2a of WFB2011
           gwood(npt)  = 0.232*(nppdm(npt)**1.165)     ! Figure 2b of WFB2011
           gcroot(npt) = 0.0348*(nppdm(npt)**1.310)    ! Figure 2d of WFB2011
           gfroot(npt) = 0.247*(nppdm(npt)**0.987)     ! Figure 2c of WFB2011
           gtot(npt)   = gleaf(npt) + gwood(npt) + gcroot(npt) + gfroot(npt)

           casaflux%fracCalloc(npt,leaf)  = gleaf(npt)/gtot(npt)
           casaflux%fracCalloc(npt,wood)  = gwood(npt)/gtot(npt)
           casaflux%fracCalloc(npt,froot) = (gcroot(npt)+gfroot(npt))/gtot(npt)

!        write(87,*) 'allocation = ',npt,casamet%iveg2(npt), totbmdm(npt),ntree(npt),nppdm(npt),casaflux%fracCalloc(npt,:)

        ELSE                ! other types
           casaflux%fracCalloc(npt,:) = casabiome%fracnpptop(veg%iveg(npt),:)
        ENDIF
    enddo

END SUBROUTINE casa_wolf


