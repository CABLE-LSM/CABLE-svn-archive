#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine COMPETE ------------------------------------------------
!
! Purpose : Updates fractional coverage of each functional type.
!           Requires a dominance hierachy as input.
!
!--------------------------------------------------------------------
SUBROUTINE compete (dom,land_pts,trif_pts,trif_index              &
,                   b,db_dfrac,forw,GAMMA,nosoil                  &
,                   frac,dfrac)

USE descent
USE seed

USE nstypes, ONLY :                                               &
!      imported scalars with intent(in)
  npft,ntype,soil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Total number of land points.
,trif_pts                                                         &
                            ! IN Number of points on which
!                                 !    TRIFFID may operate
,k,l,m,n,t                  ! WORK Loop counters.

INTEGER                                                           &
 dom(land_pts,npft)                                               &
                            ! IN Dominance hierachy.
,trif_index(land_pts)       ! IN Indices of land points on
!                                 !    which TRIFFID may operate

REAL                                                              &
 b(land_pts,npft)                                                 &
                            ! IN Mean rate of change of
!                                 !    vegetation fraction over
!                                 !    the timestep (kg C/m2/360days).
,db_dfrac(land_pts,npft,npft)                                     &
!                                 ! IN Rate of change of B
!                                 !    with vegetation fraction.
,forw                                                             &
                            ! IN Forward weighting factor.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,nosoil(land_pts)                                                 &
                            ! IN Fractional area not available
!                                 !    to vegetation.
,frac(land_pts,ntype)                                             &
                            ! INOUT Updated areal fraction.
,dfrac(land_pts,npft)                                             &
                            ! OUT Increment to areal fraction.
,denom                                                            &
                            ! WORK Denominator of update
!                                 !      equation.
,fracn,fracm                                                      &
                            ! WORK Fractions used in the spreading
!                                 !      calculation.
,numer                                                            &
                            ! WORK Numerator of the update
!                                 !      equation.
,space(land_pts)                                                  &
                            ! WORK Available space.
,p1,p2,q1,q2,r1,r2          ! WORK Coefficients in simultaneous
!                                 !      equations.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Initialisations. Set increments to zero and define the space
! available to the dominant type leaving space for the seeds of others.
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('COMPETE',zhook_in,zhook_handle)

DO t=1,trif_pts
  l=trif_index(t)
  DO n=1,npft
    dfrac(l,n) = 0.0
  END DO
  space(l) = 1-nosoil(l)-frac_min*(npft-1)
END DO

!----------------------------------------------------------------------
! Calculate the increments to the tree fractions
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  n = dom(l,1)
  m = dom(l,2)

  fracn=frac(l,n)
  fracn=MAX(fracn,frac_seed)

  fracm=frac(l,m)
  fracm=MAX(fracm,frac_seed)

  p1 = GAMMA/fracn-forw*db_dfrac(l,n,n)
  p2 = GAMMA/fracm-forw*db_dfrac(l,m,m)
  q1 = -forw*db_dfrac(l,n,m)
  q2 = -forw*db_dfrac(l,m,n)
  r1 = b(l,n)
  r2 = b(l,m)
  DO k=1,npft
    r1 = r1+forw*(db_dfrac(l,n,k)*dfrac(l,k))
    r2 = r2+forw*(db_dfrac(l,m,k)*dfrac(l,k))
  END DO

  numer = r1-(q1/p2)*r2
  denom = p1-(q1/p2)*q2
  denom = MAX(denom,denom_min)
  dfrac(l,n) = numer/denom
  frac(l,n) = frac(l,n)+dfrac(l,n)

  IF (frac(l,n) <  frac_min) THEN
    dfrac(l,n) = dfrac(l,n)+(frac_min-frac(l,n))
    frac(l,n) = frac_min
  ELSE IF (frac(l,n) >  space(l)) THEN
    dfrac(l,n) = dfrac(l,n)+(space(l)-frac(l,n))
    frac(l,n) = space(l)
  END IF

  space(l) = space(l)-frac(l,n)+frac_min

  numer = r2-q2*dfrac(l,n)
  denom = p2
  denom = MAX(denom,denom_min)
  dfrac(l,m) = numer/denom
  frac(l,m) = frac(l,m)+dfrac(l,m)

  IF (frac(l,m) <  frac_min) THEN
    dfrac(l,m) = dfrac(l,m)+(frac_min-frac(l,m))
    frac(l,m) = frac_min
  ELSE IF (frac(l,m) >  space(l)) THEN
    dfrac(l,m) = dfrac(l,m)+(space(l)-frac(l,m))
    frac(l,m) = space(l)
  END IF

  space(l) = space(l)-frac(l,m)+frac_min

END DO

!----------------------------------------------------------------------
! Calculate the increment to the shrub fraction
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  n = dom(l,3)

  fracn=frac(l,n)
  fracn=MAX(fracn,frac_seed)

  denom = GAMMA/fracn-forw*db_dfrac(l,n,n)
  denom = MAX(denom,denom_min)

  numer = b(l,n)
  DO k=1,npft
    numer = numer+forw*(db_dfrac(l,n,k)*dfrac(l,k))
  END DO

  dfrac(l,n) = numer/denom
  frac(l,n) = frac(l,n)+dfrac(l,n)

  IF (frac(l,n) <  frac_min) THEN
    dfrac(l,n) = dfrac(l,n)+(frac_min-frac(l,n))
    frac(l,n) = frac_min
  ELSE IF (frac(l,n) >  space(l)) THEN
    dfrac(l,n) = dfrac(l,n)+(space(l)-frac(l,n))
    frac(l,n) = space(l)
  END IF

  space(l) = space(l)-frac(l,n)+frac_min
END DO


!----------------------------------------------------------------------
! Calculate the increments to the grass fractions
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  n = dom(l,4)
  m = dom(l,5)

  fracn=frac(l,n)
  fracn=MAX(fracn,frac_seed)

  fracm=frac(l,m)
  fracm=MAX(fracm,frac_seed)

  p1 = GAMMA/fracn-forw*db_dfrac(l,n,n)
  p2 = GAMMA/fracm-forw*db_dfrac(l,m,m)
  q1 = -forw*db_dfrac(l,n,m)
  q2 = -forw*db_dfrac(l,m,n)
  r1 = b(l,n)
  r2 = b(l,m)
  DO k=1,npft
    r1 = r1+forw*(db_dfrac(l,n,k)*dfrac(l,k))
    r2 = r2+forw*(db_dfrac(l,m,k)*dfrac(l,k))
  END DO

  numer = r1-(q1/p2)*r2
  denom = p1-(q1/p2)*q2
  denom = MAX(denom,denom_min)
  dfrac(l,n) = numer/denom
  frac(l,n) = frac(l,n)+dfrac(l,n)

  IF (frac(l,n) <  frac_min) THEN
    dfrac(l,n) = dfrac(l,n)+(frac_min-frac(l,n))
    frac(l,n) = frac_min
  ELSE IF (frac(l,n) >  space(l)) THEN
    dfrac(l,n) = dfrac(l,n)+(space(l)-frac(l,n))
    frac(l,n) = space(l)
  END IF

  space(l) = space(l)-frac(l,n)+frac_min

  numer = r2-q2*dfrac(l,n)
  denom = p2
  denom = MAX(denom,denom_min)
  dfrac(l,m) = numer/denom
  frac(l,m) = frac(l,m)+dfrac(l,m)

  IF (frac(l,m) <  frac_min) THEN
    dfrac(l,m) = dfrac(l,m)+(frac_min-frac(l,m))
    frac(l,m) = frac_min
  ELSE IF (frac(l,m) >  space(l)) THEN
    dfrac(l,m) = dfrac(l,m)+(space(l)-frac(l,m))
    frac(l,m) = space(l)
  END IF

  space(l) = space(l)-frac(l,m)+frac_min

END DO

!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  frac(l,soil) = 1.0-nosoil(l)
  DO n=1,npft
    frac(l,soil) = frac(l,soil)-frac(l,n)
  END DO
END DO

IF (lhook) CALL dr_hook('COMPETE',zhook_out,zhook_handle)
RETURN
END SUBROUTINE compete
#endif
