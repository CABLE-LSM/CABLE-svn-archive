#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      subroutine OASIS4_um_init(comm64,ierror,cmessage)

! Description: This subroutine calls the PRISM initialization
!              routines needed when using the OASIS4
!              coupler. 
!
!              PRISM calls refer to the prebuilt OASIS4 library
!              and associated include files. The PRISM modules 
!              and library must be specified in the compile and
!              load options. FCM must be prevented from 
!              performing dependency checks on PRISM.
!
!              Link stage will also require linking to an existing
!              XML library and potentially netcdf library depending
!              on the options used to build the OASIS4 PRISM libs. 
!              (All done, typically, via a compile override file). 
!
! History:
! UM Version    Date    Description
! ----------  --------  ----------------------------------
!    6.4      Jan 2007  Original Code. R. Barnes, R. Hill
!=================================================================

      USE PRISM
! DEPENDS ON: OASIS4_atm_data_mod
      USE OASIS4_atm_data_mod

      implicit none

#include "c_kinds.h"

      integer (kind=integer64), intent(out) :: ierror
      integer (kind=integer32)              :: ierr32
      integer (kind=integer32)              :: comm32

      LOGICAL                               :: pinit

      integer (kind=integer64), intent(out) :: comm64
      character(len=*),intent(out) :: cmessage

!
! Init Component
!
      CHARACTER(LEN=128) :: model_name , comp_name
      INTEGER (kind=integer32) :: comp_id
!

      ! Check whether PRISM has been initialized
      call PRISM_initialized(pinit, ierr32)

      WRITE(6,*) "UM PRISM INITIALIZED",pinit, ierr32


      ! If PRISM not already initialised, do the necessary
      ! things at this stage.
      IF (.NOT.pinit) THEN

         ! Set the model name - hard wired to um for now.
         ! it seems unlikely that this would ever change
         ! but this situation should be kept under review
         model_name="um"
    
         call PRISM_Init (model_name, ierr32 )

         write(6,*) 'PRISM_init',model_name,ierr32

         if (ierr32 /= 0) then
            cmessage = "error in PRISM_init from init_OASIS4"
         else
            cmessage = "PRISM_init called from init_OASIS4"
         end if

         ! Set the component name. This is really something
         ! which probably ought to be user specified
         ! to be more informative (e.g. "hadgam1") rather 
         ! thean using a generic atmos.
         comp_name="atmos"

         call PRISM_Init_comp (comp_id, comp_name, ierr32 )

         write(6,*)'PRISM_init_comp',comp_id,ierr32

         if (ierr32 /= 0) then
           cmessage = "error in PRISM_init_comp from init_OASIS4"
         else
           cmessage = "PRISM_init_comp called from init_OASIS4"
         end if


          ! Find out the MPI communicator which this
          ! component must use in place of MPI_COMM_WORLD 
          call PRISM_get_localcomm(comp_id,comm32,ierr32 )
          comm64=comm32

          WRITE(6,*) "ATMOS LOCALCOMM IS",comm64 
          OASIS4_comp_id = comp_id

      END IF ! PRISM not already initialised

      ierror = ierr32

      end subroutine OASIS4_um_init
#endif
