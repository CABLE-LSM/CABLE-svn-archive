#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine LOTKA --------------------------------------------------
!
! Purpose : Updates fractional coverage of each functional type.
!           Based on the Lotka-Volterra equations of interspecies
!           competition.
!
! -------------------------------------------------------------------
SUBROUTINE lotka (land_pts,trif_pts,trif_index                    &
,                 c_veg,forw,frac_vs,frac_agric,GAMMA,lai,pc_s    &
,                 frac,dfrac)

USE nstypes, ONLY :                                               &
!      imported scalars with intent(in)
  npft,ntype

USE pftparm
USE seed
USE sigm
USE trif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Total number of land points.
,trif_pts                                                         &
                            ! IN Number of points on which
!                                 !    TRIFFID may operate
,trif_index(land_pts)                                             &
                            ! IN Indices of land points on
!                                 !    which TRIFFID may operate
,k,l,m,n,t                                                        &
                            ! WORK Loop counters.
,dom(land_pts,npft)         ! WORK Dominance hierachy.

REAL                                                              &
 c_veg(land_pts,npft)                                             &
                            ! IN Carbon content of vegetation
                            !    (kg C/m2).
,forw                                                             &
                            ! IN Forward timestep weighting.
,frac_vs(land_pts)                                                &
                            ! IN Total fractional cover of
!                                 !    vegetation and soil.
,frac_agric(land_pts)                                             &
                            ! IN Fraction of agriculture.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,lai(land_pts,npft)                                               &
                            ! IN Leaf area index.
,pc_s(land_pts,npft)                                              &
                            ! IN Net carbon flux available for
                            !    spreading (kg C/m2/360days).
,frac(land_pts,ntype)                                             &
                            ! INOUT Fractional cover of each
!                                 !       Functional Type.
,dfrac(land_pts,npft)                                             &
                            ! OUT Increment to the areal fraction
!                                 !     during the timestep (/timestep).
,b(land_pts,npft)                                                 &
                            ! WORK Mean rate of change of
!                                 !      vegetation fraction over
!                                 !      the timestep (kg C/m2/360days).
,db_dfrac(land_pts,npft,npft)                                     &
!                                 ! WORK Rate of change of B
!                                 !      with vegetation fraction.
,com(land_pts,npft,npft)                                          &
                            ! WORK Coefficients representing
!                                 !      the influence of one type
!                                 !      (second argument) on another
!                                 !      (first argument).
,diff_sum                                                         &
                            ! WORK Difference divided by sum
!                                 !      for competing canopy heights.
,hc1,hc2,hc3,hc4                                                  &
                            ! WORK Competing canopy heights (m).
,nosoil(land_pts)                                                 &
                            ! WORK Fractional area not available
!                                 !      to vegetation.
,space(land_pts,npft)       ! WORK Space available for invasion.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Define competition coefficients and the dominance hierachy
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('LOTKA',zhook_in,zhook_handle)

DO n=1,npft
  DO m=1,npft
    DO t=1,trif_pts
      l=trif_index(t)
      com(l,n,m) = 1.0
    END DO
  END DO
END DO

DO t=1,trif_pts
  l=trif_index(t)

  hc1 = a_wl(1)/(a_ws(1)*eta_sl(1))*(lai(l,1)**(b_wl(1)-1))
  hc2 = a_wl(2)/(a_ws(2)*eta_sl(2))*(lai(l,2)**(b_wl(2)-1))
  diff_sum = (hc1-hc2)/(hc1+hc2)

  com(l,1,2) = 1.0/(1+EXP(pow*diff_sum))    ! BT vs NT
  com(l,1,3) = 0.0                          ! BT vs C3G
  com(l,1,4) = 0.0                          ! BT vs C4G
  com(l,1,5) = 0.0                          ! BT vs S

  com(l,2,1) = 1.0-com(l,1,2)               ! NT vs BT
  com(l,2,3) = 0.0                          ! NT vs C3G
  com(l,2,4) = 0.0                          ! NT vs C4G
  com(l,2,5) = 0.0                          ! NT vs S

  hc3 = a_wl(3)/(a_ws(3)*eta_sl(3))*(lai(l,3)**(b_wl(3)-1))
  hc4 = a_wl(4)/(a_ws(4)*eta_sl(4))*(lai(l,4)**(b_wl(4)-1))
  diff_sum = (hc3-hc4)/(hc3+hc4)

  com(l,3,4) = 1.0/(1+EXP(pow*diff_sum))    ! C3G vs C4G
  com(l,4,3) = 1.0-com(l,3,4)               ! C4G vs C3G

  com(l,5,3) = 0.0                          ! S vs C3G
  com(l,5,4) = 0.0                          ! S vs C4G

  IF (hc1  >=  hc2) THEN
    dom(l,1) = 1
    dom(l,2) = 2
  ELSE IF (hc1  <   hc2) THEN
    dom(l,1) = 2
    dom(l,2) = 1
  END IF

  dom(l,3) = 5

  IF (hc3  >=  hc4) THEN
    dom(l,4) = 3
    dom(l,5) = 4
  ELSE IF (hc3  <   hc4) THEN
    dom(l,4) = 4
    dom(l,5) = 3
  END IF

END DO

!----------------------------------------------------------------------
! Calculate the space available for the expansion of each FT
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  nosoil(l) = 1 - frac_vs(l)
END DO

!----------------------------------------------------------------------
! Exclude non-crop types from agricultural regions
!----------------------------------------------------------------------
DO k=1,npft
  DO t=1,trif_pts
    l=trif_index(t)
    n=dom(l,k)
    space(l,n)=1.0-nosoil(l)-frac_agric(l)*(1-crop(n))            &
                            -frac_min*(npft-k)
  END DO
END DO

DO n=1,npft
  DO m=1,npft
    DO t=1,trif_pts
      l=trif_index(t)
      space(l,n)=space(l,n)-com(l,n,m)*frac(l,m)
    END DO
  END DO
END DO

!----------------------------------------------------------------------
! Calculate the variables required for the implicit calculation.
! Divide the update equation by FRAC to eliminate the (unstable)
! bare soil solution.
!----------------------------------------------------------------------
DO n=1,npft
  DO t=1,trif_pts
    l=trif_index(t)
    b(l,n) = pc_s(l,n)*space(l,n)/c_veg(l,n)-g_area(n)


    DO m=1,npft
      db_dfrac(l,n,m) = -com(l,n,m)*pc_s(l,n)/c_veg(l,n)
    END DO
  END DO
END DO

!----------------------------------------------------------------------
! Update the areal fractions
!----------------------------------------------------------------------
! DEPENDS ON: compete
CALL compete (dom,land_pts,trif_pts,trif_index                    &
,             b,db_dfrac,forw,GAMMA,nosoil                        &
,             frac,dfrac)

IF (lhook) CALL dr_hook('LOTKA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lotka
#endif
