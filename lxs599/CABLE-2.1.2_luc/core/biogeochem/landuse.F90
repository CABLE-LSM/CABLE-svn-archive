SUBROUTINE newplant(cplant_x,frac_x,ifpre_x, &
                    cplant_y,frac_y,ifpre_y,logc)
! Used for LAND USE CHANGE SIMULATION
! Call by casa_reinit
! Re-allcate plant C,N,P pools to new patch array. 
! Q.Zhang @ 29/05/2011
  USE cable_def_types_mod
  USE casadimension
  USE casaparm

  implicit none

  logical,DIMENSION(mvtype),INTENT(in) :: ifpre_x,ifpre_y
  real,DIMENSION(mvtype),INTENT(in) :: frac_x,frac_y
  real(r_2),DIMENSION(mvtype,mplant),INTENT(inout) :: cplant_x,cplant_y
  real(r_2),DIMENSION(mlogmax),INTENT(inout) :: logc
  ! local variable
  integer p 

  !!! TEST !!!!
 
  ! write(6,*)"TEST subroutine newplant called"

  !!! TEST !!!!

  DO p = 1, mvtype
     ! exist in both years    
     IF (ifpre_x(p) .and. ifpre_y(p)) THEN
        IF (abs(frac_x(p)-0.0)<1.e-8 .or. abs(frac_y(p)-0.0)<1.e-8) THEN
print *, 'Lest Veg0', p,ifpre_x(p),ifpre_y(p),frac_x(p),frac_y(p)
          STOP "vegetation fraction .eq. 0"
        END IF
        IF ((frac_x(p)-frac_y(p))>0.) THEN  ! patch weight decrease 
          ! New pools
          cplant_y(p,:) = cplant_x(p,:)
          ! Save wood log
          IF (p<=mlogmax) logc(p) = cplant_x(p,wood) * (frac_x(p) - frac_y(p))
        ELSE ! patch weight incease
          cplant_y(p,:) = cplant_x(p,:)*frac_x(p)/frac_y(p)
        END IF
   ! plant clear in the second year
     ELSEIF (ifpre_x(p) .and. .not.ifpre_y(p)) THEN
        IF (p<=mlogmax) logc(p) = cplant_x(p,wood) * (frac_x(p) - frac_y(p))
    ! does not exist in both years
     ELSE 

     END IF

  END DO  ! end pft loop

END SUBROUTINE newplant


SUBROUTINE newlitter( casabiome,frac_x,ifpre_x,frac_y,ifpre_y, & 
                      cplant_x,nplant_x,pplant_x,cplant_y,nplant_y,pplant_y, &
                      clitter_x,nlitter_x,plitter_x,clitter_y,nlitter_y,plitter_y )
! Used for LAND USE CHANGE SIMULATION
! Call by casa_reinit
! Transfer the deforest C to litter, and re-allocate litter pools. 
! Q.Zhang @ 29/05/2011
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable

  implicit none

  TYPE (casa_biome),        INTENT(IN) :: casabiome
  logical,DIMENSION(mvtype),INTENT(in) :: ifpre_x,ifpre_y
  real,DIMENSION(mvtype),INTENT(in) :: frac_x,frac_y
  real(r_2),DIMENSION(mvtype,mplant),INTENT(in) :: cplant_x,nplant_x,pplant_x
  real(r_2),DIMENSION(mvtype,mlitter),INTENT(in) :: clitter_x,nlitter_x,plitter_x
  real(r_2),DIMENSION(mvtype,mplant),INTENT(inout) :: cplant_y,nplant_y,pplant_y
  real(r_2),DIMENSION(mvtype,mlitter),INTENT(inout) :: clitter_y,nlitter_y,plitter_y 

  ! local variable
  real(r_2),DIMENSION(mvtype,mlitter,mplant) :: fromPtoL
  real(r_2),DIMENSION(mvtype,mplant) :: dcplant,dnplant,dpplant,ratioLignintoN
  real(r_2),DIMENSION(mvtype,mlitter) :: dnlitter, dplitter,clitter_g,nlitter_g,plitter_g
  real(r_2),DIMENSION(mlitter) :: dcY, dnY, dpY
  integer nL, nP, nv
  
  dcplant = 0.
  dnplant = 0.
  dpplant = 0.
  ratioLignintoN = 0.
  fromPtoL = 0.
  dnlitter = 0.
  dplitter = 0.
  dcY = 0.
  dnY = 0.
  dpY = 0.

! I. transfer removed plant to litter
  DO nP =1,mplant
    dcplant(:,nP) = cplant_x(:,nP) * frac_x(:) - cplant_y(:,nP) * frac_y(:)
    IF (icycle > 1) dnplant(:,nP) = nplant_x(:,nP) * frac_x(:) - nplant_y(:,nP) * frac_y(:)
    IF (icycle > 2) dpplant(:,nP) = pplant_x(:,nP) * frac_x(:) - pplant_y(:,nP) * frac_y(:)
  END DO
! NB: logged wood should not be transfered to litter
  dcplant(1:mlogmax,wood) = 0.
  IF (icycle > 1) dnplant(1:mlogmax,wood) = 0.
  IF (icycle > 2) dpplant(1:mlogmax,wood) = 0.

  WHERE(sum(dcplant,2) > 0.)
  ! In land use, all plant nutient is allocated to litter pools without re-asorbsion.Q.Zhang 11/08/2011 
    ratioLignintoN(:,leaf) = cplant_x(:,leaf)/max(1.0e-10,nplant_x(:,leaf)) &
                             * casabiome%fracLigninplant(:,leaf)
    ratioLignintoN(:,froot)= cplant_x(:,froot)/max(1.0e-10,nplant_x(:,froot)) &
                             * casabiome%fracLigninplant(:,froot)

    fromPtoL(:,metb,leaf)   = max(0.001, 0.85 - 0.018 *ratioLignintoN(:,leaf))
    fromPtoL(:,metb,froot)  = max(0.001, 0.85 - 0.018 *ratioLignintoN(:,froot))
    fromPtoL(:,str,leaf)    = 1.0 - fromPtoL(:,metb,leaf)
    fromPtoL(:,str,froot)   = 1.0 - fromPtoL(:,metb,froot)
    fromPtoL(:,cwd,wood)    = 1.0
  ENDWHERE

  DO nv=1, mvtype
  ! average litter pools on gridcell
     clitter_g(nv,:) = clitter_x(nv,:) *  frac_x(nv)
     nlitter_g(nv,:) = nlitter_x(nv,:) *  frac_x(nv)
     plitter_g(nv,:) = plitter_x(nv,:) *  frac_x(nv)
  ! transfer removed C,N,P pools from plant to litter 
    IF(ifpre_x(nv) .and. frac_x(nv)>frac_y(nv))THEN

      DO nL=1,mlitter
        DO nP=1,mplant
           clitter_g(nv,nL) = clitter_g(nv,nL) + fromPtoL(nv,nL,nP) * dcplant(nv,nP)
        ENDDO
      ENDDO
  
      IF(icycle > 1) THEN
         dnlitter(nv,str) = (fromPtoL(nv,str,leaf) * dcplant(nv,leaf) &
                           + fromPtoL(nv,str,froot) * dcplant(nv,froot)) * ratioNCstrfix
         dnlitter(nv,metb) = dnplant(nv,leaf) + dnplant(nv,froot) - dnlitter(nv,str)
         dnlitter(nv,CWD) = dnplant(nv,wood)
      ENDIF !end "icycle >1"
    
      IF(icycle > 2) THEN
         dplitter(nv,str) = (fromPtoL(nv,str,leaf) * dcplant(nv,leaf) &
                           + fromPtoL(nv,str,froot)* dcplant(nv,froot)) * ratioPCstrfix
         dplitter(nv,metb) = dpplant(nv,leaf) + dpplant(nv,froot) -dplitter(nv,str)
         dplitter(nv,CWD) = dpplant(nv,wood)
      ENDIF  !of "icycle >2"
    ENDIF
  END DO

  IF (icycle > 1) nlitter_g = nlitter_g + dnlitter
  IF (icycle > 2) plitter_g = plitter_g + dplitter
   
! II. re-allocate litter pools according to patch weights. 
! average pool variables from gridcell to new patches.
  DO nv=1,mvtype
    IF (ifpre_x(nv) .and. .not.ifpre_y(nv)) THEN
      dcY(:) = dcY(:) + clitter_g(nv,:)
      IF (icycle > 1) dnY(:) = dnY(:) + nlitter_g(nv,:)
      IF (icycle > 2) dpY(:) = dpY(:) + plitter_g(nv,:)
    ENDIF
  END DO

  DO nv=1,mvtype
    IF (ifpre_y(nv)) THEN   ! pft exist in the 2nd year
      clitter_y(nv,:) = clitter_g(nv,:)/frac_y(nv) + dcY(:)
      IF (icycle > 1) nlitter_y(nv,:) = nlitter_g(nv,:)/frac_y(nv) + dnY(:)
      IF (icycle > 2) plitter_y(nv,:) = plitter_g(nv,:)/frac_y(nv) + dpY(:)
    ENDIF
  END DO

END SUBROUTINE newlitter


SUBROUTINE newsoil(nd,csoil_x,frac_x,ifpre_x,csoil_y,frac_y,ifpre_y)
! Used for LAND USE CHANGE SIMULATION
! Call by casa_reinit
! Re-allocate soil C,N and P pools
! Q.Zhang @ 29/05/2011
  USE cable_def_types_mod
  USE casadimension
  USE casaparm

  implicit none

  integer,INTENT(in) :: nd  ! dimension of soil pool
  integer,DIMENSION(mvtype),INTENT(in) :: ifpre_x,ifpre_y
  real,DIMENSION(mvtype),INTENT(in) :: frac_x,frac_y
  real(r_2),DIMENSION(mvtype,nd),INTENT(inout) :: csoil_x,csoil_y
! local variable
  real(r_2),DIMENSION(nd) :: dcY
  integer nv

  dcY = 0. 

! plant clear, pft patch 'disappear'
  DO nv=1,mvtype
    IF (ifpre_x(nv) .and. .not.ifpre_y(nv)) THEN
      dcY(:) = dcY(:) + csoil_x(nv,:) * frac_x(nv)
    ENDIF
  END DO
! Q.Zhang: should not put the bellow IF() into the above DO loop. 
  DO nv=1,mvtype
    IF (ifpre_y(nv)) THEN   ! pft exist in the 2nd year
      csoil_y(nv,:) = csoil_x(nv,:) * frac_x(nv) / frac_y(nv) + dcY(:)
    ENDIF
  END DO

END SUBROUTINE newsoil
