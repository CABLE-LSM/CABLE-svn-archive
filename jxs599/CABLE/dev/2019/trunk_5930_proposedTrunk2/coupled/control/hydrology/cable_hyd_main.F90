module cable_hyd_main_mod
  
contains

SUBROUTINE cable_hyd_main( land_pts, ntiles, lying_snow, SNOW_surft, SURF_ROFF,&
                           SUB_SURF_ROFF, TOT_TFALL )
  
  !subrs called 
  USE cable_hyd_driv_mod, ONLY : cable_hyd_driver
  
  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
  USE cable_data_module, ONLY : cable

!H!  !diag 
!H!  USE cable_fprint_module, ONLY : cable_fprintf
!H!  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
!H!  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
!H!                                 L_cable_Pyfprint, unique_subdir
  USE cable_def_types_mod, ONLY : mp !only need for fprint here
  
!data
USE cbl_masks_mod, ONLY : L_tile_pts
 
  implicit none
 
  !___ re-decl input args

  integer :: land_pts, ntiles
  
  real :: snow_surft(land_pts,ntiles)
  
  real, dimension(land_pts) ::                                                 &
    lying_snow,    & ! OUT Gridbox snowmass (kg/m2)        
    sub_surf_roff,  & ! OUT Sub-surface runoff (kg/m2/s).
    surf_roff,      & ! OUT Surface runoff (kg/m2/s).
    tot_tfall         ! OUT Total throughfall (kg/m2/s).

  !___ local vars
  logical, save :: first_call = .true.
  
  ! std template args 
  character(len=*), parameter :: subr_name = "cable_hyd_main"
 
!H!# include "../../../core/utils/diag/cable_fprint.txt"
  
  !-------- Unique subroutine body -----------
  
  !--- initialize cable_runtime% switches 
  cable_runtime%um =          .TRUE.
  cable_runtime%um_hydrology =.TRUE.
  
  CALL cable_hyd_driver( land_pts, ntiles, L_tile_pts, lying_snow, SNOW_surft, &
                         SURF_ROFF, SUB_SURF_ROFF, TOT_TFALL )
  
  cable_runtime%um_hydrology =.FALSE.
  
  !-------- End Unique subroutine body -----------
  
!H!  fprintf_dir=trim(fprintf_dir_root)//trim(unique_subdir)//"/"
!H!  if(L_cable_fprint) then 
!H!    !basics to std output stream
!H!    if (knode_gl == 0 .and. ktau_gl == 1)  call cable_fprintf(subr_name, .true.) 
!H!    !more detailed output
!H!    vname=trim(subr_name//'_')
!H!    call cable_fprintf( cDiag00, vname, knode_gl, ktau_gl, .true. )
!H!  endif
!H!
!H!  if(L_cable_Pyfprint .and. ktau_gl == 1) then 
!H!    !vname='latitude'; dimx=mp
!H!    !call cable_Pyfprintf( cDiag1, vname, cable%lat, dimx, .true.)
!H!  endif

  first_call = .false.        

return

End subroutine cable_hyd_main
  
End module cable_hyd_main_mod

