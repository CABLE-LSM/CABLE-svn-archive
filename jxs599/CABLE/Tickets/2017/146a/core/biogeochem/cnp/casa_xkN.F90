SUBROUTINE casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
! computing the reduction in litter and SOM decomposition
! when decomposition rate is N-limiting
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(INOUT) :: xkNlimiting
  TYPE (casa_pool),         INTENT(INOUT) :: casapool
  TYPE (casa_flux),         INTENT(INOUT) :: casaflux
  TYPE (casa_met),          INTENT(INOUT) :: casamet
  TYPE (casa_biome),        INTENT(INOUT) :: casabiome
!
  TYPE (veg_parameter_type),   INTENT(IN) :: veg  ! vegetation parameters

  ! local variables
  INTEGER i,j,k,kk,iv,thepoint,nland
  REAL(r_2), DIMENSION(mp)         :: xFluxNlittermin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilmin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilimm
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilminnet
! A maximum Clitter set to avoid positive feedback for litter accumulation
! when N mode is activated. (Q.Zhang 23/05/2011)
!  real(r_2), dimension(17)         :: xClitter
!  data xClitter/100.0,100.0,100.0,100.0,50.0,150.0,150.0,100.0,&
!                150.0,150.0,100.0, 20.0,20.0, 20.0, 20.0, 20.0,20.0/

  xkNlimiting  = 1.0
!  set N mineral N fluxes to zero
  xFluxNlittermin(:)  = 0.0
  xFluxNsoilmin(:)    = 0.0
  xFluxNsoilimm(:)    = 0.0  !negative for microbial upatek and postive for release of mineral N
  xFluxNsoilminnet(:) = 0.0
!   PRINT *, 'within casa_xkN'

!  calculate gross mineralisation
  DO nland=1,mp
  IF (casamet%iveg2(nland)/=icewater) THEN

    ! calculate C:N ratio of newly formed SOM as function of soil mineral N pool
    IF (casapool%Nsoilmin(nland) < 2.0) THEN
      casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmin(nland,:)  &
                                       + (casapool%ratioNCsoilmax(nland,:) &
                                         -casapool%ratioNCsoilmin(nland,:)) &
                                       * max(0.0,casapool%Nsoilmin(nland)) / 2.0
    ELSE
      casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmax(nland,:)
    ENDIF

    DO j=1,mlitter
      xFluxNlittermin(nland) = xFluxNlittermin(nland) &
                         + casaflux%klitter(nland,j) * casapool%Nlitter(nland,j)
    ENDDO
    DO k=1,msoil
      xFluxNsoilmin(nland)   = xFluxNsoilmin(nland)   &
                         + casaflux%ksoil(nland,k)   * casapool%Nsoil(nland,k)
    ENDDO

    ! calculate N immobilisation from L to S and S to S
    DO kk=1,msoil
      DO j=1,mlitter    ! immobilisation from litter to soil
        xFluxNsoilimm(nland) = xFluxNsoilimm(nland) &
             - casaflux%fromLtoS(nland,kk,j) * casaflux%klitter(nland,j) &
             * casapool%Clitter(nland,j) * casapool%ratioNCsoilnew(nland,kk)
      ENDDO
      DO k=1,msoil      ! immobilisation from soil to soil
        IF(k.ne.kk) THEN
          xFluxNsoilimm(nland) = xFluxNsoilimm(nland) &
               - casaflux%fromStoS(nland,kk,k) * casaflux%ksoil(nland,k) &
               * casapool%Csoil(nland,k) * casapool%ratioNCsoilnew(nland,kk)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  ENDDO

  ! now check if there is sufficient mineral N
  xFluxNsoilminnet(:) = xFluxNlittermin(:) + xFluxNsoilmin(:) + xFluxNsoilimm(:)
!   PRINT *, 'casamet%iveg2 = ', casamet%iveg2
!   PRINT *, 'deltpool = ',deltpool
!   PRINT *, 'xFluxNsoilminnet = ', xFluxNsoilminnet
! WHERE(casamet%iveg2(:)/=icewater)
!    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0)
!      xkNlimiting(:) =1.0
!    ELSEWHERE
!      xkNlimiting(:) =max(0.0, - (casapool%Nsoilmin(:)-2.0)/(deltpool*xFluxNsoilminnet(:)))
!      xkNlimiting(:) =MIN(1.0,xkNlimiting(:))
!    ENDWHERE
! ENDWHERE

! Q.Zhang 23/05/2011 test code according to YPW
  WHERE(casamet%iveg2(:)/=icewater)
    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0 &
          .OR. xFluxNsoilminnet(:) .ge. 0.0)
      xkNlimiting(:) =1.0
    ELSEWHERE
      xkNlimiting(:) =MAX(0.0, - (casapool%Nsoilmin(:)-0.5) &
                                /(deltpool*xFluxNsoilminnet(:)))
      xkNlimiting(:) =MIN(1.0,xkNlimiting(:))
    ENDWHERE
! Q.Zhang 23/05/2011 test
! If pool size larger than xClitter, turnover rate will not constrained by Nsoilmin.
!    where(casapool%clitter(:,1) > xClitter(veg%iveg(:)))
!     xkNlimiting(:) = 1.0
!    end where
! end (Q.Zhang 23/05/2011)
    where(sum(casapool%clitter,2) > casabiome%maxfinelitter(veg%iveg(:)) + casabiome%maxcwd(veg%iveg(:)))
     xkNlimiting(:) = 1.0
    end where
  ENDWHERE


END SUBROUTINE casa_xkN

