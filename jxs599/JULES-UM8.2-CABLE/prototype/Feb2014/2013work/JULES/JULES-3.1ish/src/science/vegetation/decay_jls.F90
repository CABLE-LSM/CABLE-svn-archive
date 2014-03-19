#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine DECAY --------------------------------------------------
!
! Purpose : Updates carbon contents of the soil.
!
!--------------------------------------------------------------------
SUBROUTINE decay (land_pts,trif_pts,trif_index                    &
,                 dpc_dcs,forw,GAMMA,pc,cs)

USE descent
USE csmin

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
,l,t,n                      ! WORK Loop counters

REAL                                                              &
 dpc_dcs(land_pts,4)                                              &
                            ! IN Rate of change of PC with
!                                 !    soil carbon (yr).
,forw                                                             &
                            ! IN Forward timestep weighting.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,pc(land_pts,4)                                                   &
                            ! IN Net carbon flux into the
!                                 !    soil (kg C/m2/360days).
,cs(land_pts,4)                                                   &
                            ! INOUT Soil carbon (kg C/m2).
,denom                                                            &
                            ! WORK Denominator of update
!                                 !      equation.
,numer                      ! WORK Numerator of the update
!                                 !      equation.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('DECAY',zhook_in,zhook_handle)

DO n=1,4
  DO t=1,trif_pts
    l=trif_index(t)

    numer = pc(l,n)
    denom = GAMMA+forw*dpc_dcs(l,n)
    denom = MAX(denom,denom_min)

    cs(l,n) = cs(l,n)+numer/denom

    cs(l,n) = MAX(cs_min,cs(l,n))

  END DO
END DO

IF (lhook) CALL dr_hook('DECAY',zhook_out,zhook_handle)
RETURN
END SUBROUTINE decay
#endif
