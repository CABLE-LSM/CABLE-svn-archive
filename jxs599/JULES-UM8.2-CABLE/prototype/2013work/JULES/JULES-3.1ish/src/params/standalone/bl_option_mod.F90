#if !defined(UM_JULES)
! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

! **********************************************************************
! THIS VERSION OF bl_option_mod IS STRICTLY FOR USE WITH
! STANDALONE JULES!
!***********************************************************************

!+ Data module for switches/options concerned with the BL scheme.
! Description:
!   Module containing runtime options/data used by the boundary
!   layer scheme.

! Method:
!   Switches and associated data values used by the boundary layer
!   scheme are defined here and assigned default values. These may
!   be overridden by namelist input.

!   Any routine wishing to use these options may do so with the 'USE'
!   statement.

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Boundary Layer

! Code Description:
!   Language: FORTRAN 90

MODULE bl_option_mod

! Declarations:

IMPLICIT NONE

!======================================================================= 
!   Permissible settings for BL options.
!=======================================================================
INTEGER, PARAMETER :: off = 0  ! Switch disabled
INTEGER, PARAMETER :: on  = 1  ! Switch enabled


!     Options for form drag
INTEGER, PARAMETER :: No_drag         = 0
INTEGER, PARAMETER :: Effective_z0    = 1
INTEGER, PARAMETER :: Explicit_stress = 2

!     Options for marine boundary layers
INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory

!     Options for surface exchange
INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
INTEGER, PARAMETER :: Limit_expl_ustar = 2
!       Option under the COR_UST switch to limit the magnitude of the
!       explicitly calculated ustar
INTEGER, PARAMETER :: IP_SrfExWithCnv = 1
!       Option to include deep convective gustiness in the surface
!       transfer

INTEGER :: ISrfExCnvGust = 0
                       ! Switch to include the effect of convective
                       ! downdraughts on surface exchange
                       ! OFF (=0) => not used: only boundary-layer
                       !   gustiness is considered (original version)
                       ! IP_SrfExWithCnv (=1) the impact of gustiness
                       !   due to boundary layer eddies is reduced
                       !   relative to the above, but eddies driven
                       !   by convective downdraughts are included

INTEGER :: formdrag = No_drag
                       ! switch for orographic form drag
INTEGER :: FD_stab_dep = off
                       ! switch to implement stability
                       ! dependence of orographic form drag
INTEGER :: ISeaZ0T = Fixed_Z0T
                       ! option for specifying the thermal roughness length
                       !       at sea points
                       
REAL :: max_stress_grad = 0.05
                       ! Maximum implied stress gradient across the
                       ! boundary layer, used to limit the explicit
                       ! stress applied in non-local scalings (m/s2)

REAL :: orog_drag_param = 0.3
                       ! Drag coefficient for orographic form drag   
                       
REAL :: SeaSalinityFactor = 1.0
                       ! Scaling of qsat allowing for salinity of sea water

! CHARNOCK is a constant in the Charnock formula for sea-surface
!          roughness length for momentum (Z0MSEA).
REAL :: charnock = 0.011


END MODULE bl_option_mod

#endif
