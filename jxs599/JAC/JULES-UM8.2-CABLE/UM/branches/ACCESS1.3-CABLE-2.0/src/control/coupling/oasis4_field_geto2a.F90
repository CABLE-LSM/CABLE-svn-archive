#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine OASIS4_field_geto2a(                    &
     &  g_p_field,                                       &
#include "argd1.h" 
     &  cmessage)

! DEPENDS ON: OASIS4_atm_data_mod
      USE OASIS4_atm_data_mod 


      Implicit none
!
! Description: This subroutine controls the receiving
!              of incoming coupling data to be coupled via OASIS4
!              from other component models or to be read from files.
!                                   
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.4    Jan 2007  Original code. R. Hill, M. Christoforou 
!================================================================

! Declarations:
#include "parvars.h"   
! Includes an include that defines halo_type
#include "decomptp.h"  
! Contains definition of decomp_standard_atm
#include "decompdb.h"  
! decomp_db_glsize
#include "atm_lsm.h"   
! Contains definitions for atmos_landmask_lo

#include "cmaxsize.h"  
!
#include "typsize.h"   
!
#include "typd1.h"     
! d1
#include "typptra.h"   
! Contains jfrac_land but TODO don't know wh
#include "typcona.h"   
!
#include "typlndm.h"   
!

#include "caoptr.h"    
! Contains ja_solar etc. set in inita2nemo.
#include "c_lheat.h"   
! Contains code that set LC (Latent Heat Con
#include "c_mdi.h"     
! Contains definition and setting of RMDI
#include "cntlatm.h"   
! l_ctile
#include "cntlocn.h"   
!

!     Subroutine arguments
      Integer :: g_p_field      ! Global horiz domain for atmos

!     Local dynamic arrays:
!     Coupling fields on atmosphere grid being sent into NEMO ocean model
!     One dimensional versions of array for setting from d1 array.

      Logical :: amasktp(g_p_field) ! Atmos model land-sea mask for TP
!     Two dimensional versons that we will dynamically allocate 
      Real,dimension(:,:),allocatable :: ocn_sst
      Real,dimension(:,:),allocatable :: ocn_sst_orig
      Real,dimension(:,:),allocatable :: ocn_freeze
      Real,dimension(:,:),allocatable :: ocn_hice
      Real,dimension(:,:),allocatable :: ocn_u
      Real,dimension(:,:),allocatable :: ocn_v

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: put_data

      Real    :: fland_loc(g_p_field)
      Real    :: fland_glo(g_p_field)
      Integer :: nlandpt

      Real :: latentHeatOfCond=lc

      Integer :: i              ! Loop counter
      Integer :: j              ! Loop counter
      Integer :: o4_jmt
      Integer :: o4_jmt_u
      Integer :: o4_jmt_v
      Integer :: o4_imt
      Integer :: o4error,ft
      Integer :: icode          ! Error return code (=0 is OK)
      Integer :: o4info         ! Return code from prism calls
      Integer :: infosst, infou, infov, infofreeze
      Integer :: gather_pe      ! Processor for gathering
      Character*(*) :: cmessage ! OUT - Error return message
 
      Integer :: POINTER  ! Pointer to D1 for incoming data fields
      Integer :: POINTER_HICE  ! Pointer to D1 for incoming ice
                               ! depth - allows more than one D1
                               ! field to be updated in the same loop.

      ! Set up local domain sizes for t, u  and v points.
      o4_imt=lasize(1,fld_type_p,halo_type_no_halo)
      o4_jmt=lasize(2,fld_type_p,halo_type_no_halo)
      o4_jmt_u=lasize(2,fld_type_u,halo_type_no_halo)
      o4_jmt_v=lasize(2,fld_type_v,halo_type_no_halo)

      allocate(ocn_sst(o4_imt,o4_jmt))
      allocate(ocn_sst_orig(o4_imt,o4_jmt))
      allocate(ocn_freeze(o4_imt,o4_jmt))
      allocate(ocn_hice(o4_imt,o4_jmt))
      allocate(ocn_u(o4_imt,o4_jmt_u))   ! u grid points
      allocate(ocn_v(o4_imt,o4_jmt_v))   ! v grid points

      ocn_sst = 0.0
      ocn_u = 0.0
      ocn_v =0.0
      ocn_freeze = 0.0
      ocn_hice = 0.0

      ! Move all existing field values to our incoming array
      ! this saves us having to do an explicit mask because
      ! OASIS4 will just overlay the sea points 
      ! NOTE: for pre 2007 versions this cannot be relied upon
      ! when applying a mask in the atmos SMIOC to 
      ! SST fields we get correct values over sea but zero
      ! over land i.e. it doesn't simply not overwrite
      ! land values but actively writes zeros which is just
      ! plain wrong! It is hoped this will be corrected 
      ! with the conservative coupler.
      ! The UM expects values in K (not Celsius)


      IF (L_CTILE) THEN
         POINTER = JTSTAR_SEA
      ELSE
         POINTER = JTSTAR
      END IF

      DO J = 1, o4_jmt
         DO I = 1, o4_imt
            ocn_sst(I,J)=D1(pointer)
            ocn_sst_orig(I,J) = D1(pointer)
            POINTER = POINTER+1
         END DO
      END DO

      ! Prepare the ice fraction field (0->1)
      POINTER = JA_AICE
      DO J = 1, o4_jmt
         DO I = 1, o4_imt
            ocn_freeze(I,J)=D1(pointer)
            POINTER = POINTER+1
         END DO
      END DO


      ! Incoming currents expected in m/s
      POINTER = JU_SEA
      DO J = 1, o4_jmt_u
         DO I = 1, o4_imt
            ! Don't preset for now ocn_u(I,J) = D1(pointer)
            POINTER = POINTER+1
         END DO
      END DO

      POINTER = JV_SEA
      DO J = 1, o4_jmt_v
         DO I = 1, o4_imt
            ! Don't preset for now ocn_v(I,J) = D1(pointer)
            POINTER = POINTER+1
         END DO
      END DO

      ft = fld_type_p 
      var_ind = vind_ocn_sst
! DEPENDS ON: OASIS4_get64
      CALL OASIS4_get64(ocn_sst,o4_imt,o4_jmt,infosst,o4error,ft   &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)
   
      var_ind = vind_ocn_freeze
! DEPENDS ON: OASIS4_get64
      CALL OASIS4_get64(ocn_freeze,o4_imt,o4_jmt,infofreeze,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

 
      ft = fld_type_u 
      var_ind = vind_ocn_u
! DEPENDS ON: OASIS4_get64
      CALL OASIS4_get64(ocn_u,o4_imt,o4_jmt_u,infou,o4error,ft    &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)


      ft = fld_type_v
      var_ind = vind_ocn_v
! DEPENDS ON: OASIS4_get64
      CALL OASIS4_get64(ocn_v,o4_imt,o4_jmt_v,infov,o4error,ft    &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

      
      ! Having got our new values we need to move them
      ! to D1 where they'll be picked up next time the appropriate
      ! calculation is called. 
      IF (L_CTILE) THEN
         POINTER = JTSTAR_SEA
      ELSE
         POINTER = JTSTAR
      END IF

      ! We should ultimately make the following D1 updates dependent 
      ! on INFOU etc being equal to PRISM_Cpl or similar.
      ! However, we don't have these PRISM values immediately available
      ! yet. 


      ! This is a temporary hack to get sensible
      ! numbers in the polar row -  it is NOT intended to be 
      ! a permanent solution once the conservative coupler is 
      ! available. 
      IF (MYPE.EQ.NPROC-1) THEN
         DO I = 1, o4_imt
            ocn_sst(i,o4_jmt)=ocn_sst(i,o4_jmt-1)
            ocn_freeze(i,o4_jmt)=ocn_freeze(i,o4_jmt-1)
         END DO
      END IF

      ! If this PE deals with the North polar row we must
      ! set a uniform temperature in all grid points along that 
      ! row to prevent crashes in the atmos dynamics.
      ! That's all jolly well in a 1xn decomposition, but a 
      ! nx1 or nxm composition makes this more complicated. 
      ! Any implications for conservation - there shouldn't be really.
      ! The easiest way to do this is with a 
      ! gather-mean-scatter operation. Of course we have to do this
      ! before assigning to D1. 
      ! Note: MEAN_POLAR_ROW must be called for all T grid point 
      ! fields which we receive.
!      IF (INFOSST.EQ.PRISM_Cpl) THEN

! DEPENDS  ON: MEAN_POLAR_ROW
          CALL MEAN_POLAR_ROW(OCN_SST)
 
          ! The incoming SST is expected in K NOT degrees C
          ! Any necessary transforms will be done in the PSMILE
          ! controlled via the SMIOC.
          DO J = 1, o4_jmt
             DO I = 1, o4_imt
                IF (.NOT.UM_TMASK(I,J,1)) THEN
! Special form of masking for T points - if
! value returned from coupling is approx abs zero
! we assume this is an unchanged point masked during
! the coupling process and we reinstate
! our original value.
                   IF (ocn_sst(I,J).LT.1.0) THEN
                      ocn_sst(I,J)=ocn_sst_orig(I,J)
                   END IF
                   D1(pointer)=ocn_sst(I,J)
                END IF 
                POINTER = POINTER+1
             END DO
          END DO
!      END IF

       ! Update ice fraction
!      IF (INFOFREEZE.EQ.PRISM_Cpl) THEN

! DEPENDS ON: MEAN_POLAR_ROW
          CALL MEAN_POLAR_ROW(OCN_FREEZE) 

          POINTER = JA_AICE
          POINTER_HICE = JICE_THICKNESS
          ocn_hice(:,:) = ocn_freeze(:,:)
          DO J = 1, o4_jmt
             DO I = 1, o4_imt
                IF (.NOT.UM_TMASK(I,J,1)) THEN
                   D1(pointer)=ocn_freeze(I,J)
                   D1(pointer_hice)=ocn_hice(I,J)
                END IF
                POINTER = POINTER+1
                POINTER_HICE = POINTER_HICE+1
             END DO
          END DO
!      END IF

       ! Update our surface currents (U)
       POINTER = JU_SEA
!      IF (INFOU.EQ.PRISM_Cpl) THEN
         ! Incoming currents expected in m/s

         ! This is a temporary hack to get vaguely sensible
         ! numbers in the polar row -  it is NOT to be submitted in
         ! any permanent base code. This is a fix
         ! to make up for shortcomings in OASIS4 remapping
         IF (MYPE.EQ.NPROC-1) THEN
            DO I = 1, o4_imt
               ocn_u(i,o4_jmt)=0.0
            END DO
         END IF

         DO J = 1, o4_jmt_u
            DO I = 1, o4_imt
               IF (.NOT.UM_UMASK(I,J,1)) THEN
                  D1(pointer)=ocn_u(I,J)
               END IF
               POINTER = POINTER+1
            END DO
         END DO
!      END IF

      ! Update our surface currents (V)
      POINTER = JV_SEA

!      IF (INFOV.EQ.PRISM_Cpl) THEN
         ! Incoming currents expected in m/s
         DO J = 1, o4_jmt_v
            DO I = 1, o4_imt
               IF (UM_VMASK(I,J,1)) THEN
                  D1(pointer)=ocn_v(I,J)
               END IF
               POINTER = POINTER+1
            END DO
         END DO
!      END IF


      ! We don't need these arrays any more, deallocate them
      ! to avoid memory leaks.
      deallocate(ocn_v)
      deallocate(ocn_u)
      deallocate(ocn_hice)
      deallocate(ocn_freeze)
      deallocate(ocn_sst)
      deallocate(ocn_sst_orig)

      End Subroutine OASIS4_field_geto2a
#endif
