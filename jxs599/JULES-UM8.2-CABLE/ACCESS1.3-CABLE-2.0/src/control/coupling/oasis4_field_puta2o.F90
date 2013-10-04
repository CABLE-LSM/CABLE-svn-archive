#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine OASIS4_field_puta2o(                                              &  
     &  g_p_field,                                                    &                          
#include "argd1.h" 
     &  cmessage)

! DEPENDS ON: OASIS4_atm_data_mod
      USE OASIS4_atm_data_mod 

      Implicit none
!
! Description: This subroutine controls the setting up and passing
!              of outgoing coupling data to be coupled via OASIS4
!              to other component models or to be output to files.
!                                   
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.4    Jan 2007  Original code. R. Hill, M. Christoforou 
!================================================================

!
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
! d1?
#include "typptra.h"   
! Contains jfrac_land but TODO don't know wh
#include "typcona.h"   
!
#include "typlndm.h"   
!

#include "caoptr.h"    
! Contains ja_solar etc. set in inita2NEMO.
#include "c_lheat.h"   
! Contains code that set LC (Latent Heat Con
#include "c_mdi.h"     
! Contains defination and setting of RMDI
#include "cntlatm.h"   
! l_ctile
#include "cntlocn.h"   
!

!     Subroutine arguments
      Integer :: g_p_field      ! Global horiz domain for atmos


!     Local dynamic arrays:
!     Coupling fields on atmosphere grid being sent to ocean model.
!     One dimensional versions of array for setting from d1 array.

      Logical :: amasktp(g_p_field) ! Atmos model land-sea mask for TP.
!     Two dimensional versions that we will dynamically allocate


      Real,dimension(:,:),allocatable :: taux
      Real,dimension(:,:),allocatable :: tauy

      Real,dimension(:,:),allocatable :: solar2d
      Real,dimension(:,:),allocatable :: blue2d
      Real,dimension(:,:),allocatable :: evap2d
      Real,dimension(:,:),allocatable :: longwave2d
      Real,dimension(:,:),allocatable :: sensible2d
      Real,dimension(:,:),allocatable :: heatflux

      Real,dimension(:,:),allocatable :: rainls
      Real,dimension(:,:),allocatable :: snowls
      Real,dimension(:,:),allocatable :: rainconv
      Real,dimension(:,:),allocatable :: snowconv
      Real,dimension(:,:),allocatable :: pme 
      Real,dimension(:,:),allocatable :: empmr
      Real,dimension(:,:),allocatable :: riverout
      Real,dimension(:,:),allocatable :: wme

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
      Integer :: o4error, ft
      Integer :: icode          ! Error return code (=0 is OK)
      Integer :: info           ! Return code from MPP
      Integer :: gather_pe      ! Processor for gathering
      Character*(*) :: cmessage ! OUT - Error return message

!     Set the arrays from the d1 array. Note that ja_solar etc. are 
!     defined in include file caoptr.

!     Our row lengths and column lengths must be the local domain sizes
!     for each PE. We don't gather and scatter. That's done by OASIS4. 

      o4_imt=lasize(1,fld_type_p,halo_type_no_halo) !L_ROW_LENGTH
      o4_jmt=lasize(2,fld_type_p,halo_type_no_halo) !L_P_ROWS
      o4_jmt_u=lasize(2,fld_type_u,halo_type_no_halo) ! U rows
      o4_jmt_v=lasize(2,fld_type_v,halo_type_no_halo) ! V rows

    
!     Now that we know what the dimensions are going to be, allocate 
!     the 2d arrays.


      allocate(taux(o4_imt,o4_jmt_u))
      allocate(tauy(o4_imt,o4_jmt_v))

      allocate(solar2d(o4_imt,o4_jmt))
      allocate(blue2d(o4_imt,o4_jmt))
      allocate(evap2d(o4_imt,o4_jmt))
      allocate(longwave2d(o4_imt,o4_jmt))
      allocate(sensible2d(o4_imt,o4_jmt))     

      allocate(rainls(o4_imt,o4_jmt))
      allocate(snowls(o4_imt,o4_jmt))
      allocate(rainconv(o4_imt,o4_jmt))
      allocate(snowconv(o4_imt,o4_jmt))
      allocate(pme(o4_imt,o4_jmt))
      allocate(empmr(o4_imt,o4_jmt))

      allocate(riverout(o4_imt,o4_jmt)) 
      allocate(wme(o4_imt,o4_jmt)) 

      ! Prepare fields for putting to the coupler
      ! by copying to our temporary arrays.
      ! In principle, we could probably just pass 
      ! D1 with an appropriate pointer for some fields
      ! to the put, but we still have to combine
      ! many fields for e.g. in getting the PME 
      ! field etc.
      ! Doing it this way is slightly long winded
      ! but much clearer.

      ! Set up fields from the U grid
      Do j=1,o4_jmt_u
         Do i=1,o4_imt
            taux(i,j)=D1(ja_taux+i-1+((j-1)*o4_imt))
         End do
      End do

      ! Set up fields from the V grid
      Do j=1,o4_jmt_v
         Do i=1,o4_imt
            tauy(i,j)=D1(ja_tauy+i-1+((j-1)*o4_imt))
         End do
      End do

      ! Set up fields from the T grid
      Do j=1,o4_jmt
         Do i=1,o4_imt
            ! Copy various heat flux components 
            solar2d(i,j)=D1(ja_solar+i-1+((j-1)*o4_imt))
            blue2d(i,j)=D1(ja_blue+i-1+((j-1)*o4_imt))
            evap2d(i,j)=D1(ja_evap+i-1+((j-1)*o4_imt))
            longwave2d(i,j)=D1(ja_longwave+i-1+((j-1)*o4_imt))
            sensible2d(i,j)=D1(ja_sensible+i-1+((j-1)*o4_imt)) 

            ! PME components
            rainls(i,j) = D1(ja_lsrain+i-1+((j-1)*o4_imt))
            snowls(i,j) = D1(ja_lssnow+i-1+((j-1)*o4_imt))
            rainconv(i,j) = D1(ja_cvrain+i-1+((j-1)*o4_imt))
            snowconv(i,j) = D1(ja_cvsnow+i-1+((j-1)*o4_imt))

            ! River runoff
            riverout(i,j) = D1(ja_riverout+i-1+((j-1)*o4_imt))
            ! Wind mixing energy WME
            wme(i,j) = D1(ja_windmix+i-1+((j-1)*o4_imt))
         End do
      End do


!     NON-PENETRATIVE SURFACE HEAT FLUXES.
!
!     NOTICE THAT THE BLUE END OF THE SOLAR SPECTRUM HAS TO BE
!     SUBTRACTED OUT HERE. IT SHOULD ALSO BE POINTED OUT THAT AT
!     SEA-ICE POINTS (IF THEY EXIST), THE SENSIBLE HEAT FLUX AND
!     THE EVAPORATION WERE ALREADY WEIGHTED BY THE FRACTIONAL AREA
!     OF LEADS WHEN THEY WERE DIAGNOSED, SO NO SPECIAL CODE IS
!     NECESSARY HERE.
!       


      ! Allocate array which actually gets used for putting data.
      allocate(heatflux(o4_imt,o4_jmt))     

      Do j=1,o4_jmt
         Do i=1,o4_imt
            ! We need to mask out land points
            ! in heat flux otherwise we end up using MDIs
            If (UM_TMASK(i,j,1)) Then
               heatflux(i,j)=0.0
            Else
! NEMO on the other hand needs Total Heat rather than non-pen solar
! plus heat flux: Effectively this just means we no longer subtract
! the sw (blue) term.
               heatflux(i,j)=solar2d(i,j)+longwave2d(i,j)         &
     &           -(sensible2d(i,j)+latentHeatOfCond*evap2d(i,j))
            End if

! Calculate total PME. What about Sea ice contribution?
! How do we add that because we dont have ocean ice concs etc.
            pme(i,j)=SNOWLS(I,J) + SNOWCONV(I,J)

! IF SEAICE ACTIVE:
!           pme(i,J)=pme(i,J)*(1.0 - AICE(I,J))
!
!
            pme(i,j)=pme(i,j) + RAINLS(I,J) + RAINCONV(I,J)       &
     &                          - EVAP2d(I,J)

! The thing is that NEMO refers to emp NOT pme so 
! as a minimum we have to reverse the sign. (* by -1.0)
! We could, however do this via "mult_scalar" in the SMIOC.
! so no particular special code needed for that, but !

! This is wrong for NEMO bulk momthly forcing mode
! we don't do anything with runoff - that's treated separately
!           empmr(i,j)=riverout(i,j)-pme(i,j)

            empmr(i,j)=pme(i,j)
         End do
      End do

      ! Now we have the fields which are going to be put
      ! it just remains to issue the puts

      ft = fld_type_p
      var_ind = vind_pen_solar
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_put64(blue2d,o4_imt,o4_jmt,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

       ! Put total heat flux
      var_ind = vind_tot_atm_hflux
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_put64(heatflux,o4_imt,o4_jmt,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

       ! Put total EMP
      var_ind = vind_pme
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_put64(empmr,o4_imt,o4_jmt,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

       ! Put riverrunoff
      var_ind = vind_runoff 
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_put64(riverout,o4_imt,o4_jmt,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

      ft = fld_type_u
      ! Put x windstress
      var_ind = vind_taux
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_put64(taux,o4_imt,o4_jmt_u,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

      ft = fld_type_v

      ! Put y windstress
      var_ind = vind_tauy
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_put64(tauy,o4_imt,o4_jmt_v,o4error,ft  &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

      ft = fld_type_p
      ! Put wme
      var_ind = vind_wme
! DEPENDS ON: OASIS4_put64
      CALL OASIS4_PUT64_(wme,o4_imt,o4_jmt,o4error,ft   &
     &,halo_type_no_halo,gc_all_proc_group,nproc,mype)

      ! Deallocate all the arrays which only need to 
      ! exist in this routine.
      deallocate(wme)
      deallocate(riverout)

      deallocate(empmr)
      deallocate(pme)
      deallocate(snowconv)
      deallocate(rainconv)
      deallocate(snowls)
      deallocate(rainls)
      deallocate(heatflux)
      deallocate(sensible2d)
      deallocate(longwave2d)
      deallocate(evap2d)
      deallocate(blue2d)
      deallocate(solar2d)

      deallocate(taux)
      deallocate(tauy)

      End Subroutine OASIS4_field_puta2o
#endif
