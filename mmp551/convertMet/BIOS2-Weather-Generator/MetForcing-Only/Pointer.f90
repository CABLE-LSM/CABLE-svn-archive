MODULE PointerModule
!-------------------------------------------------------------------------------
! * This module defines and assigns pointers to specific variables in WaterDynStep
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! * Pointers for time variables (array TTime)
real(sp),pointer:: Time1y       ! [y]           ! 01 run time at end of step
real(sp),pointer:: TTDay        ! [d]           ! 02 current date: day
real(sp),pointer:: TTMonth      ! [mth]         ! 03 current date: month
real(sp),pointer:: TTYear       ! [y]           ! 04 current date: year
real(sp),pointer:: TTEndMth     ! [-]           ! 05 (0,1) = last day of month (N,Y)
! * Pointers for stores (array XX)
real(sp),pointer:: WRelA(:)     ! [0,1]         ! 01 rel soil water (layer A)
real(sp),pointer:: WRelB(:)     ! [0,1]         ! 02 rel soil water (layer 2B)
real(sp),pointer:: CLea(:)      ! [molC/m2]     ! 03 live leaf C
real(sp),pointer:: WRel1(:)     ! [0,1]         ! 04 rel soil water (layer 1)
real(sp),pointer:: WRel2(:)     ! [0,1]         ! 05 rel soil water (layer 2)
real(sp),pointer:: WRel3(:)     ! [0,1]         ! 06 rel soil water (layer 3)
real(sp),pointer:: WRel4(:)     ! [0,1]         ! 07 rel soil water (layer 4)
real(sp),pointer:: WRel5(:)     ! [0,1]         ! 08 rel soil water (layer 5)
real(sp),pointer:: WRel6(:)     ! [0,1]         ! 09 rel soil water (layer 6)
real(sp),pointer:: WRel7(:)     ! [0,1]         ! 10 rel soil water (layer 7)
real(sp),pointer:: WRel8(:)     ! [0,1]         ! 11 rel soil water (layer 8)
real(sp),pointer:: WRel9(:)     ! [0,1]         ! 12 rel soil water (layer 9)
real(sp),pointer:: WRel10(:)     ! [0,1]         ! 13 rel soil water (layer 10)
real(sp),pointer:: Tsoil1(:)    ! [0,1]         ! 14 soil T (layer 1)
real(sp),pointer:: Tsoil2(:)    ! [0,1]         ! 15 soil T (layer 2)
real(sp),pointer:: Tsoil3(:)    ! [0,1]         ! 16 soil T (layer 3)
real(sp),pointer:: Tsoil4(:)    ! [0,1]         ! 17 soil T (layer 4)
real(sp),pointer:: Tsoil5(:)    ! [0,1]         ! 18 soil T (layer 5)
real(sp),pointer:: Tsoil6(:)    ! [0,1]         ! 19 soil T (layer 6)
real(sp),pointer:: Tsoil7(:)    ! [0,1]         ! 20 soil T (layer 7)
real(sp),pointer:: Tsoil8(:)    ! [0,1]         ! 21 soil T (layer 8)
real(sp),pointer:: Tsoil9(:)    ! [0,1]         ! 22 soil T (layer 9)
real(sp),pointer:: Tsoil10(:)    ! [0,1]         ! 23 soil T (layer 10)
real(sp),pointer:: cplant1(:)   ! [gC/m2]       ! 24 leaf C
real(sp),pointer:: cplant2(:)   ! [gC/m2]       ! 25 wood C
real(sp),pointer:: cplant3(:)   ! [gC/m2]       ! 26 root C
real(sp),pointer:: csoil1(:)    ! [gC/m2]       ! 27 soil microbial biomass
real(sp),pointer:: csoil2(:)    ! [gC/m2]       ! 28 soil fast C pool
real(sp),pointer:: cansto(:)    ! [mm]           ! 29 water stored in canopy
real(sp),pointer:: theta1(:)     ! [0,1]         ! 30 vol soil water (layer 1)
real(sp),pointer:: theta2(:)     ! [0,1]         ! 31 vol soil water (layer 2)
real(sp),pointer:: theta3(:)     ! [0,1]         ! 32 vol soil water (layer 3)
real(sp),pointer:: theta4(:)     ! [0,1]         ! 33 vol soil water (layer 4)
real(sp),pointer:: theta5(:)     ! [0,1]         ! 34 vol soil water (layer 5)
real(sp),pointer:: theta6(:)     ! [0,1]         ! 35 vol soil water (layer 6)
real(sp),pointer:: theta7(:)     ! [0,1]         ! 36 vol soil water (layer 7)
real(sp),pointer:: theta8(:)     ! [0,1]         ! 37 vol soil water (layer 8)
real(sp),pointer:: theta9(:)     ! [0,1]         ! 38 vol soil water (layer 9)
real(sp),pointer:: theta10(:)     ! [0,1]         ! 39 vol soil water (layer 10)
real(sp),pointer:: wcol(:)     ! [m]         ! 240 total soil water column (layers 1-10)

! * Pointers for fluxes and dXXdt (array FF)
real(sp),pointer:: dWcolAdt(:)  ! [0,1]/d       ! 01 dXX/dt:  soil water (A)
real(sp),pointer:: dWcolBdt(:)  ! [0,1]/d       ! 02 dXX/dt:  soil water (B)
real(sp),pointer:: FWPrec(:)    ! [m/d]         ! 03 precipitation
real(sp),pointer:: FWTraA(:)    ! [m/d]         ! 04 transpiration (layer 1)
real(sp),pointer:: FWTraB(:)    ! [m/d]         ! 05 transpiration (layer 2)
real(sp),pointer:: FWSoil(:)    ! [m/d]         ! 06 soil evaporation
real(sp),pointer:: FWRun(:)     ! [m/d]         ! 07 runoff
real(sp),pointer:: FWLchA(:)    ! [m/d]         ! 08 leach ex layer 1
real(sp),pointer:: FWLchB(:)    ! [m/d]         ! 09 leach ex layer 2 (deep drainage)
real(sp),pointer:: FWE(:)       ! [m/d]         ! 10 total evaporation
real(sp),pointer:: FWTra(:)     ! [m/d]'        ! 11 total transpiration
real(sp),pointer:: FWDis(:)     ! [m/d]'        ! 12 total point discharge (runoff+lch)
real(sp),pointer:: FWThrough(:) ! [m/d]'        ! 13 throughfall
real(sp),pointer:: FWwc(:)      ! [m/d]'        ! 14 flux from wet canopy
real(sp),pointer:: dWCanopydt(:)! [m/d]'        ! 15 dXX/dt: stored water in canopy
real(sp),pointer:: FCGPP(:)     ! [molC/m2/d]   ! 16 photsynthetic uptake (GPP)
real(sp),pointer:: FCGro(:)     ! [molC/m2/d]   ! 17 Growth C flux (NPP)
real(sp),pointer:: FCNEE(:)     ! [molC/m2/d]   ! 18 Net Ecosystem Exchange (NEE)
real(sp),pointer:: phiE(:)       ! [W/m2]         ! 19 daytime latent heat flux
real(sp),pointer:: phiH(:)       ! [W/m2]         ! 20 daytime sensible heat flux
real(sp),pointer:: FWPt(:)     ! [m/d]         ! 21 Preistley-Taylor (not calc in CABLE)
real(sp),pointer:: dCLeadt(:)     ! [molC/m2/d]         ! 22
real(sp),pointer:: dWcoldt(:)  ! [0,1]/d       ! 23 dXX/dt: change in total soil water 

! * Pointers for derived quantities (array DD)
real(sp),pointer:: DayltFrac(:) ! [-]           ! 01 daylight fraction
real(sp),pointer:: RhoA(:)      ! [molA/m3]     ! 02 air density
real(sp),pointer:: FracVExt(:)  ! [-]           ! 03 veg cover fraction (from external file)
real(sp),pointer:: rLAIExt(:)   ! [-]           ! 04 leaf area index    (from external file)
real(sp),pointer:: FracVCLea(:) ! [-]           ! 05 veg cover fraction (from CLea)
real(sp),pointer:: rLAICLea(:)  ! [-]           ! 06 leaf area index    (from CLea)
real(sp),pointer:: AllocLea(:)  ! [-]           ! 07 leaf allocation coeff
real(sp),pointer:: FCGroL(:)    ! [molC/m2/d]   ! 08 light-lim NPP
real(sp),pointer:: FCGroW(:)    ! [molC/m2/d]   ! 09 water-lim NPP
real(sp),pointer:: ImBalA(:)    ! [m/d]         ! 10 dWA/dt-(FWPrec-FWTraA-FWSoil-FWRun-FWLchA)
real(sp),pointer:: ImBalB(:)    ! [m/d]         ! 11 dWB/dt-(FWLchA-FWTraB-FWLchB)
real(sp),pointer:: ImBalCanopy(:)! [m/d]        ! 12 Canopy Imbalance
real(sp),pointer:: ImBalSoil(:)  ! [m/d]        ! 13   Soil Imbalance
real(sp),pointer:: ImBal(:)      ! [m/d]        ! 14 Total moisture imbalance
real(sp),pointer:: fws(:)        ! [-]          ! 15 limitation to stomatal conductance due to soil moisture deficit
real(sp),pointer:: Tabar(:)      ! [degC]       ! 16 mean air temperature (for CASA-CNP)
real(sp),pointer:: Tsoilbar_g(:) ! [degC]       ! 17 root-weighted mean soil temp, grass (for CASA-CNP)
real(sp),pointer:: Tsoilbar_w(:) ! [degC]       ! 18 root-weighted mean soil temp, wood (for CASA-CNP)
real(sp),pointer:: Sbar_g(:)     ! [0;1]        ! 19 root-weighted soil moisture relative to saturation, grass (for CASA-CNP)
real(sp),pointer:: Sbar_w(:)     ! [0;1]        ! 20 root-weighted soil moisture relative to saturation, wood (for CASA-CNP)
real(sp),pointer:: btran_g(:)    ! [0;1]        ! 21 water limitation factor for carbon allocation, grass (for CASA-CNP)
real(sp),pointer:: btran_w(:)    ! [0;1]        ! 22 water limitation factor for carbon allocation, wood (for CASA-CNP)
real(sp),pointer:: FCGPP_g(:)    ! [molC/m2]/d  ! 23 GPP , grass (per unit tile area) (for CASA-CNP)
real(sp),pointer:: FCGPP_w(:)    ! [molC/m2]/d  ! 24 GPP , wood (per unit tile area) (for CASA-CNP)
real(sp),pointer:: FCLeafR_g(:)  ! [molC/m2]/d  ! 25 day-time leaf respiration , grass (per unit tile area) (for CASA-CNP)
real(sp),pointer:: FCLeafR_w(:)  ! [molC/m2]/d  ! 26 GPP , day-time leaf respiration, grass (per unit tile area) (for CASA-CNP)
real(sp),pointer:: LAI_g(:)  ! [] ! 27
real(sp),pointer:: LAI_w(:)  ! [] ! 28
real(sp),pointer:: fws_g(:)  ! [] ! 29
real(sp),pointer:: fws_w(:)  ! [] ! 30
! * Pointers for predicted point-scale observables (array ZZP)
real(sp),pointer:: ZNDVI(:)      ! [degC]        ! 01 NDVI
real(sp),pointer:: ZphiRnet(:)   ! [MJ m-2 d-1]  ! 02 day-time Rnet
real(sp),pointer:: ZphiH(:)      ! [MJ m-2 d-1]  ! 03 day-time phiH
real(sp),pointer:: ZphiE(:)      ! [MJ m-2 d-1]  ! 04 day-time phiE
real(sp),pointer:: ZphiNEE(:)    ! [moles CO2 m-2 d-1]  ! 05 day-time phiNEE
real(sp),pointer:: ZphiNPP(:)    ! [moles CO2 m-2 d-1]  ! 06 day-time phiNPP
real(sp),pointer:: ZphiGPP(:)    ! [moles CO2 m-2 d-1]  ! 07 day-time phiGPP
real(sp),pointer:: Zsm0008(:)    ! [m3 m-3]       !08 volumetric soil moisture content 0-8 cm
real(sp),pointer:: Zsm0090(:)    ! [m3 m-3]       !09 volumetric soil moisture content 0-90 cm
real(sp),pointer:: Zsm0030(:)    ! [m3 m-3]       !10 volumetric soil moisture content 0-30 cm
real(sp),pointer:: Zsm3060(:)    ! [m3 m-3]       !11 volumetric soil moisture content 30-60 cm
real(sp),pointer:: Zsm6090(:)    ! [m3 m-3]       !12 volumetric soil moisture content 60-90 cm
! * Pointers for actual point-scale observations (array AAP)
real(sp),pointer:: ANDVI(:)      ! [degC]        ! 01 NDVI
real(sp),pointer:: AphiRnet(:)   ! [MJ m-2 d-1]  ! 02 day-time Rnet
real(sp),pointer:: AphiH(:)      ! [MJ m-2 d-1]  ! 03 day-time phiH
real(sp),pointer:: AphiE(:)      ! [MJ m-2 d-1]  ! 04 day-time phiE
real(sp),pointer:: AphiNEE(:)    ! [moles CO2 m-2 d-1]  ! 05 day-time phiNEE
real(sp),pointer:: AphiNPP(:)    ! [moles CO2 m-2 d-1]  ! 06 day-time phiNPP
real(sp),pointer:: AphiGPP(:)    ! [moles CO2 m-2 d-1]  ! 07 day-time phiGPP
real(sp),pointer:: Asm0008(:)    ! [m3 m-3]       !08 volumetric soil moisture content 0-15 cm
real(sp),pointer:: Asm0090(:)    ! [m3 m-3]       !09 volumetric soil moisture content 0-90 cm
real(sp),pointer:: Asm0030(:)    ! [m3 m-3]       !10 volumetric soil moisture content 0-30 cm
real(sp),pointer:: Asm3060(:)    ! [m3 m-3]       !11 volumetric soil moisture content 30-60 cm
real(sp),pointer:: Asm6090(:)    ! [m3 m-3]       !12 volumetric soil moisture content 60-90 cm
real(sp),pointer:: ALAIg(:)       ![] ! 13 grassy leaf area index (per unit tile area)
real(sp),pointer:: ALAIw(:)		![] ! 14 woody leaf area index (per unit tile area)
real(sp),pointer:: AfWoody(:)	![] ! 15 fraction cover woody
real(sp),pointer:: AscattVIS(:)	![] ! 16 monthly scattering coefficient (vis)
real(sp),pointer:: AscattNIR(:)	![] ! 17 monthly scattering coefficient (NIR)

! * Pointers for predicted point-scale observables (array ZZPh)
real(sp),pointer:: ZLST(:,:)      ! [degC]       ! 01 LST at overpass time
real(sp),pointer:: ZLSTTime(:,:)  ! [hrSol]      ! 02 overpass time for LST
real(sp),pointer:: ZLSTAngle(:,:) ! [deg]        ! 03 view angle for LST
real(sp),pointer:: DelTsTaZ(:,:)  ! [degC]       ! 04 DelTsTaA = ZLST - TempAt
real(sp),pointer:: TaZ(:,:)       ! [degC]       ! 05 TempAt
real(sp),pointer:: ZphiRneth(:,:) ! [Wm-2]       ! 06 phiRneth
real(sp),pointer:: ZphiHh(:,:)    ! [Wm-2]       ! 07 phiHh
real(sp),pointer:: ZphiEh(:,:)    ! [Wm-2]       ! 08 phiEh
real(sp),pointer:: ZphiNEEh(:,:)  ! [umol m-2 s-1]! 09 phiNEEh
real(sp),pointer:: ZphiGh(:,:)    ! [Wm-2]       ! 10 phiGh
real(sp),pointer:: ZTsoilh(:,:)    ! [deg C]       ! 11 Tsoilh
real(sp),pointer:: ZphiHhtime(:,:)! [h]          ! 12 phiHhtime
real(sp),pointer:: TrVeg(:,:)     ! [degC]       ! 13 Radiative temperature of veg
real(sp),pointer:: TrSoil(:,:)    ! [degC]       ! 14 Radiative temperature of soil
real(sp),pointer:: TrEff(:,:)     ! [degC]       ! 15 Effective Radiative temperature of surface (depends on view angle)
real(sp),pointer:: TSoil1h(:,:)     ! [degC]       ! 16 Soil T
real(sp),pointer:: TSoil2h(:,:)     ! [degC]       ! 17 Soil T
real(sp),pointer:: TSoil3h(:,:)     ! [degC]       ! 18 Soil T
real(sp),pointer:: TSoil4h(:,:)     ! [degC]       ! 19 Soil T
real(sp),pointer:: TSoil5h(:,:)     ! [degC]       ! 20 Soil T
real(sp),pointer:: TSoil6h(:,:)     ! [degC]       ! 21 Soil T
real(sp),pointer:: TSoil7h(:,:)     ! [degC]       ! 22 Soil T
real(sp),pointer:: TSoil8h(:,:)     ! [degC]       ! 23 Soil T
real(sp),pointer:: TSoil9h(:,:)     ! [degC]       ! 24 Soil T
real(sp),pointer:: TSoil10h(:,:)     ! [degC]       ! 25 Soil T
real(sp),pointer:: theta1h(:,:)     ! []       ! 26 Soil moisture
real(sp),pointer:: theta2h(:,:)     ! []       ! 27 Soil moisture
real(sp),pointer:: theta3h(:,:)     ! []       ! 28 Soil moisture
real(sp),pointer:: theta4h(:,:)     ! []       ! 29 Soil moisture
real(sp),pointer:: theta5h(:,:)     ! []       ! 30 Soil moisture
real(sp),pointer:: theta6h(:,:)     ! []       ! 31 Soil moisture
real(sp),pointer:: theta7h(:,:)     ! []       ! 32 Soil moisture
real(sp),pointer:: theta8h(:,:)     ! []       ! 33 Soil moisture
real(sp),pointer:: theta9h(:,:)     ! []       ! 34 Soil moisture
real(sp),pointer:: theta10h(:,:)     ! []       ! 35 Soil moisture

! * Pointers for actual point-scale observations (array AAPh)
real(sp),pointer:: ALST(:,:)      ! [degC]        ! 01 LST at overpass time
real(sp),pointer:: ALSTTime(:,:)  ! [hrSol]       ! 02 overpass time for LST
real(sp),pointer:: ALSTAngle(:,:) ! [deg]         ! 03 view angle for LST
real(sp),pointer:: DelTsTaA(:,:)  ! [degC]        ! 04 DelTsTaA = ALST - TempAt
real(sp),pointer:: TempAt(:,:)    ! [degC]        ! 05 TempAt
real(sp),pointer:: AphiRneth(:,:) ! [Wm-2]        ! 06 AphiRneth
real(sp),pointer:: AphiHh(:,:)    ! [Wm-2]        ! 07 AphiHh
real(sp),pointer:: AphiEh(:,:)    ! [Wm-2]        ! 08 AphiEh
real(sp),pointer:: AphiNEEh(:,:)  ! [umol m-2 s-1]! 09 AphiNEEh
real(sp),pointer:: AphiGh(:,:)    ! [Wm-2]        ! 10 AphiGh
real(sp),pointer:: ATsoilh(:,:)    ! [degC]        ! 11 AphiTsoilh
real(sp),pointer:: AphiHhtime(:,:)! [h]           ! 12 AphiHhtime
! * Pointers for predicted catchment-scale observables (array ZZC)
real(sp),pointer:: ZDisCD(:)    ! [m/d]         ! 01 catchment daily outflow
real(sp),pointer:: ZDisCM(:)    ! [m/mth]       ! 02 catchment monthly outflow
real(sp),pointer:: ZRunCD(:)    ! [m/d]         ! 03 catchment daily runoff
real(sp),pointer:: ZLchCD(:)    ! [m/d]         ! 04 catchment daily leaching
! * Pointers for actual catchment-scale observations (array AAC)
real(sp),pointer:: ADisCD(:)    ! [m/d]         ! 01 catchment daily outflow
real(sp),pointer:: ADisCM(:)    ! [m/mth]       ! 02 catchment monthly outflow
real(sp),pointer:: ARunCD(:)    ! [m/d]         ! 03 catchment daily runoff
real(sp),pointer:: ALchCD(:)    ! [m/d]         ! 04 catchment daily leaching
! * Pointers for met forcing variables (array MM)
real(sp),pointer:: SolarMJ(:)   ! [MJ/d]        ! 01 incident solar radn
real(sp),pointer:: Precip(:)    ! [m/d]         ! 02 precipitation
real(sp),pointer:: TempMax(:)   ! [degC]        ! 03 24hr-max air temp
real(sp),pointer:: TempMin(:)   ! [degC]        ! 04 24hr-min air temp
real(sp),pointer:: vph09(:)     ! [mb]          ! 05 24hr-min air temp
real(sp),pointer:: vph15(:)     ! [mb]          ! 06 24hr-min air temp
! * Pointers for hourly met  forcing variables (array hMM)
real(sp),pointer:: hFsd(:,:)       ! [W m-2]       ! 01 incident solar radn
real(sp),pointer:: hFld(:,:)       ! [W m-2]       ! 02 down-welling lowg-wave
real(sp),pointer:: hPrecip(:,:)    ! [mm/dt]       ! 03 precip
real(sp),pointer:: hUa(:,:)        ! [ms-1]        ! 04 windspeed at ref height
real(sp),pointer:: hTc(:,:)        ! [degC]        ! 05 T at ref height
real(sp),pointer:: hqv(:,:)        ! [kg/kg]       ! 06 specific humidity at ref height
real(sp),pointer:: hpmb(:,:)       ! [mb]          ! 07 atmospheric pressure
real(sp),pointer:: hcoszen(:,:)       ! [mb]       ! 08 cos (solar zenith angle)
! * Pointers for remote sensing driver variables (array RR)
real(sp),pointer:: LAIg(:)   ! []        ! 01 grassy LAI 
real(sp),pointer:: LAIw(:)    ! []         ! 02 woody LAI
real(sp),pointer:: fWoody(:)   ! []        ! 03 fraction woody cover
real(sp),pointer:: scattVIS(:)   ! []        ! 04 leaf scattering coefficient (vis)
real(sp),pointer:: scattNIR(:)     ! []          ! 05 leaf scattering coeffcient (NIR)
! * Pointers for spatially uniform parameters (array UU)
real(sp),pointer:: CoeffPT      ! [-]           ! 01 Priestley-Taylor coeff
real(sp),pointer:: CoeffBeer    ! [-]           ! 02 Beer Law extinction coeff
real(sp),pointer:: CLea0        ! [-]           ! 03 CLea0 = RhoCLeaf*LeafThick
real(sp),pointer:: RateCLea     ! [1/d]         ! 04 rate constant for CLea decay
real(sp),pointer:: AllocLg      ! []         ! 05 allocation of C to leaves (grassy)
real(sp),pointer:: AllocLw     ! []         ! 06 allocation of C to leaves (woody)
real(sp),pointer:: HyConSat1    ! [m/d]         ! 07 sat hydraulic conductivity (1)
real(sp),pointer:: HyConSat2    ! [m/d]         ! 08 sat hydraulic conductivity (2)
real(sp),pointer:: PwrFWSoil    ! [-]           ! 09 FWSoil ~ WRel**(PwrFWSoil+1)
real(sp),pointer:: PwrFWLch     ! [-]           ! 10 FWLch  ~ WRel**(PwrFWLch+1)
real(sp),pointer:: alfaQ        ! [molC/molQ]   ! 11 Light Use Efficiency
real(sp),pointer:: alfaWpri     ! [molC/W]      ! 12 prior Water Use Eff   (use if > 0)
real(sp),pointer:: alfaWmul     ! [-]           ! 13 mult for WUE from def (use if > 0)
real(sp),pointer:: CO2A         ! [molC/A]      ! 14 air [CO2]
real(sp),pointer:: WRelA0       ! [-]           ! 15 scale for AllocLea(WRel)
real(sp),pointer:: ZSoil1Mult   ! [-]           ! 16 multiplier for ZSoil1
real(sp),pointer:: ZSoil2Mult   ! [-]           ! 17 multiplier for ZSoil2
real(sp),pointer:: rLAImax      ! [-]           ! 18 maximum rLAI
real(sp),pointer:: Gaero        ! [m/s]         ! 19 reference aerodynamic conductance
real(sp),pointer:: TimeTxFrac   ! [-]           ! 20 Tmax time as frac of (TDawn,TDusk)
real(sp),pointer:: cN0          ! [-]           ! 21 ZNDVI = Sum[cNi*FracV^i]
real(sp),pointer:: cN1          ! [-]           ! 22 ZNDVI = Sum[cNi*FracV^i]
real(sp),pointer:: cN2          ! [-]           ! 23 ZNDVI = Sum[cNi*FracV^i]
real(sp),pointer:: TZRunDef     ! [d]           ! 24 default catchment runoff timescale (use if >0)
real(sp),pointer:: TZLchDef     ! [d]           ! 25 default catchment leach  timescale (use if >0)
real(sp),pointer:: CoeffPAR     ! [-]           ! 26 FracV = CoeffPAR*FAPAR
real(sp),pointer:: RatioJV      ! [-]           ! 27 RatioJV=Jmax/Vcmax
real(sp),pointer:: za           ! [m]           ! 28 reference height
real(sp),pointer:: ratecp1      ! [1/year]		! 29 rate constant: plant carbon pool 1
real(sp),pointer:: ratecp2      ! [1/year]		! 30 rate constant: plant carbon pool 2
real(sp),pointer:: ratecp3      ! [1/year]		! 31 rate constant: plant carbon pool 3
real(sp),pointer:: ratecs1      ! [1/year]		! 32 rate constant: soil carbon pool 1
real(sp),pointer:: ratecs2      ! [1/year]		! 33 rate constant: soil carbon pool 2
real(sp),pointer:: zeta      ! [-]		    ! 34 macropore flow parameter
real(sp),pointer:: fsatmax      ! [-]		! 35 Multiplier for litter depth
real(sp),pointer:: B1Mult      ! [-]		    ! 36 Multiplier for b1
real(sp),pointer:: Psie1Mult      ! [-]		! 37 Multiplier for psie1
real(sp),pointer:: Ksat1Mult      ! [-]		! 38 Multiplier for Ksat
real(sp),pointer:: B2Mult      ! [-]		    ! 39 Multiplier for b1
real(sp),pointer:: Psie2Mult      ! [-]		! 40 Multiplier for psie1
real(sp),pointer:: Ksat2Mult      ! [-]		! 41 Multiplier for Ksat
real(sp),pointer:: dleaf_g	    ![m]		! 42 leaf length (grass)
real(sp),pointer:: vcmax_g	    ![mol m-2 s-1] ! 43 maximum RuBP carboxylation rate top leaf (grass)
real(sp),pointer:: hc_g	    ![m] ! 44 canopy ht (grass)
real(sp),pointer:: xfang_g	    ![-] ! 45 leaf angle dist parameter (grass)
real(sp),pointer:: rp20_g	    ![m] ! 46 relative plant respiration coefficient (grass)
real(sp),pointer:: vbeta_g	    ![m] ! 47 stomatal sensitivity to soil water (grass)
real(sp),pointer:: rootbeta_g	    ![m] ! 48 parameter to describe root density distribution (grass)
real(sp),pointer:: F10_g	    ![m] ! 49 fraction of roots in top 10 cm (grass)
real(sp),pointer:: ZR_g	    ![m] ! 50 maximum rooting depth (grass)
real(sp),pointer:: loggamma_g	    ![-] ! 51 parameter in root efficiency function (grass)
real(sp),pointer:: dleaf_w	    ![m]		! 52 leaf length (woody)
real(sp),pointer:: vcmax_w	    ![mol m-2 s-1] ! 53 maximum RuBP carboxylation rate top leaf (woody)
real(sp),pointer:: hc_w	    ![m] ! 54 canopy ht (woody)
real(sp),pointer:: xfang_w	    ![-] ! 55 leaf angle dist parameter (woody)
real(sp),pointer:: rp20_w	    ![m] ! 56 relative plant respiration coefficient (woody)
real(sp),pointer:: vbeta_w	    ![m] ! 57 stomatal sensitivity to soil water (woody)
real(sp),pointer:: rootbeta_w	    ![m] ! 58 parameter to describe root density distribution (woody)
real(sp),pointer:: F10_w	    ![m] ! 59 fraction of roots in top 10 cm (woody)
real(sp),pointer:: ZR_w	    ![m] ! 60 maximum rooting depth (woody)
real(sp),pointer:: loggamma_w	    ![-] ! 61 parameter in root efficiency function (woody)
real(sp),pointer:: a1	    ![-] ! 62 parameter in stomatal conductance function
real(sp),pointer:: ds0	    ![-] ! 63 sensitivity of stomatal conduxtance to VPD


! * Pointers for spatially variable parameters (array VV)
real(sp),pointer:: CatchMap(:)  ! [-]           ! 01 CRCCH catchment ID map
real(sp),pointer:: Albedo(:)    ! [-]           ! 02 albedo
real(sp),pointer:: LatDeg(:)    ! [degN]        ! 03 latitude
real(sp),pointer:: LongDeg(:)    ! [degE]        ! 04 longitude
real(sp),pointer:: Altitude(:)  ! [m]           ! 05 altitude = elevation
real(sp),pointer:: FAPAR(:,:)   ! [-]           ! 06-17 observed FAPAR: Jan-Dec
real(sp),pointer:: canst1(:)    ! [mm/LAI]      ! 18 max intercepted water by canopy (mm/LAI)
real(sp),pointer:: dleaf(:)     ! [m]           ! 19 leaf length (m)
real(sp),pointer:: vcmax(:)     ! [mol m-2 s-1] ! 20 maximum RuBP carboxylation rate top leaf (mol/m2/s)
real(sp),pointer:: hc(:)        ! [m]           ! 21 canopy height (m)
real(sp),pointer:: xfang(:)     ! [-]           ! 22 leaf angle dist
real(sp),pointer:: rp20(:)      ! [-]           ! 23 relative plant respiration coefficient
real(sp),pointer:: rpcoef(:)    ![1/degC]       ! 24 temperature coef nonleaf plant respiration 
real(sp),pointer:: rs20(:)     ! [-]           ! 25 relative soil respiration at 20C (tuning parameter)
real(sp),pointer:: shelrb(:)    ! [-]           ! 26 sheltering factor (dimensionless)
real(sp),pointer:: frac4(:)     ! [-]           ! 27 fraction of c4 plant
real(sp),pointer:: tminvj(:)    ! [degC]        ! 28 min temperature of the start of photosynthesis
real(sp),pointer:: tmaxvj(:)    ! [degC]        ! 29 min temperature of the start of photosynthesis
real(sp),pointer:: vbeta(:)     ! [-]           ! 30 stomatal sensitivity to soil water
real(sp),pointer:: albsoil(:)   ! [-]           ! 31 soil albedo
real(sp),pointer:: swilt1(:)    ! [-]           ! 32 vol H2O @ wilting (A horizon)
real(sp),pointer:: sfc1(:)      ! [-]           ! 33 vol H2O @ field capacity (A horizon)
real(sp),pointer:: WVolSat1(:)  ! [-]           ! 34 vol H2O @ saturation (A horizon)
real(sp),pointer:: bch1(:)      ! [-]           ! 35 parameter b in Campbell equation (A horizon)
real(sp),pointer:: hyds1(:)     ! [m/s]			! 36 hydraulic conductivity @ saturation (A horizon)
real(sp),pointer:: sucs1(:)     ! [m]			! 37 suction at saturation (A horizon)
real(sp),pointer:: rhosoil1(:)  ! [kg/m3]		! 38 soil density (A horizon)
real(sp),pointer:: css1(:)      ! [kJ/kg/K]		! 39 thermal conductivity of dry soil (A horizon)
real(sp),pointer:: clay1(:)    ! [-]			! 40 fraction clay (A horizon)
real(sp),pointer:: silt1(:)    ! [%]			! 41 fraction silt(A horizon)
real(sp),pointer:: ZSoil1(:)    ! [m]			! 42 water-holding soil depth (A horizon)
real(sp),pointer:: swilt2(:)    ! [-]           ! 43 vol H2O @ wilting (B horizon)
real(sp),pointer:: sfc2(:)      ! [-]           ! 44 vol H2O @ field capacity (B horizon)
real(sp),pointer:: WVolSat2(:)  ! [-]           ! 45 vol H2O @ saturation (B horizon)
real(sp),pointer:: bch2(:)      ! [-]           ! 46 parameter b in Campbell equation (B horizon)
real(sp),pointer:: hyds2(:)     ! [m/s]         ! 47 hydraulic conductivity @ saturation (B horizon)
real(sp),pointer:: sucs2(:)     ! [m]			! 48 suction at saturation (B horizon)
real(sp),pointer:: rhosoil2(:)  ! [kg/m3]		! 49 soil density (B horizon)
real(sp),pointer:: css2(:)      ! [kJ/kg/K]		! 50 thermal conductivity of dry soil (B horizon)
real(sp),pointer:: ZSoil2(:)    ! [m]			! 51 water-holding soil depth (layer 1) (B horizon)
real(sp),pointer:: clay2(:)     ! [m]			! 52 fraction clay (B horizon)
real(sp),pointer:: silt2(:)     ! [m]			! 53 fraction silt (B horizon)
real(sp),pointer:: rootbeta(:)	    ![-]		! 54 parameter to describe root density distribution
real(sp),pointer:: extkn(:)		![-]		    ! 55 extinction coefficient for nitrogen with canopy depth
real(sp),pointer:: vegcf(:)		![-]		    ! 56  biome-specific soil respiration rate
real(sp),pointer:: iveg(:)		![-]		    ! 57 IGBP vegetation type
real(sp),pointer:: LAI(:,:)   ! [-]             ! 58-69 LAI: Jan-Dec
real(sp),pointer:: clittequil(:)    ! [tC/ha]        ! 70 litter pool (from BIOS)
real(sp),pointer:: albsoilVIS(:) ! [-]  ! 71 soil albedo (vis)
real(sp),pointer:: albsoilNIR(:) ! [-]  ! 72 soil albedo (NIR)

CONTAINS

!*******************************************************************************

SUBROUTINE PointAll (TargetTTime, TargetXX, TargetFF, TargetDD,     &
                     TargetZZP, TargetAAP,TargetZZPh, TargetAAPh,   &
					 TargetZZC, TargetAAC,    &
                     TargetMM,TargetRR,TargetUU, TargetVV,TargethMM)
!-------------------------------------------------------------------------------
! Assign all pointers from generic to specific array names
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
! Input variables
real(sp),intent(in),target,optional:: TargetTTime(:)
real(sp),intent(in),target,optional:: TargetXX(:,:)
real(sp),intent(in),target,optional:: TargetFF(:,:)
real(sp),intent(in),target,optional:: TargetDD(:,:)
real(sp),intent(in),target,optional:: TargetZZP(:,:)
real(sp),intent(in),target,optional:: TargetAAP(:,:)
real(sp),intent(in),target,optional:: TargetZZPh(:,:,:)
real(sp),intent(in),target,optional:: TargetAAPh(:,:,:)
real(sp),intent(in),target,optional:: TargetZZC(:,:)
real(sp),intent(in),target,optional:: TargetAAC(:,:)
real(sp),intent(in),target,optional:: TargetMM(:,:)
real(sp),intent(in),target,optional:: TargethMM(:,:,:)
real(sp),intent(in),target,optional:: TargetRR(:,:)
real(sp),intent(in),target,optional:: TargetUU(:)
real(sp),intent(in),target,optional:: TargetVV(:,:)
!-------------------------------------------------------------------------------
! * Pointers for time variables (array TTime)
if (present(TargetTTime)) then
  if (size(TargetTTime) /= 5) stop "PointAll: wrong TTime dimension"
  Time1y     => TargetTTime(01)
  TTDay      => TargetTTime(02)
  TTMonth    => TargetTTime(03)
  TTYear     => TargetTTime(04)
  TTEndMth   => TargetTTime(05)
end if
! * Pointers for stores (array XX)
if (present(TargetXX)) then
  if (size(TargetXX,2) /= 40) stop "PointAll: wrong nxx in XX"
  WRelA      => TargetXX(:,01)
  WRelB      => TargetXX(:,02)
  CLea       => TargetXX(:,03)
  WRel1      => TargetXX(:,04)
  WRel2      => TargetXX(:,05)
  WRel3      => TargetXX(:,06)
  WRel4      => TargetXX(:,07)
  WRel5      => TargetXX(:,08)
  WRel6      => TargetXX(:,09)
  WRel7      => TargetXX(:,10)
  WRel8      => TargetXX(:,11)
  WRel9      => TargetXX(:,12)
  WRel10      => TargetXX(:,13)
  Tsoil1     => TargetXX(:,14)
  Tsoil2     => TargetXX(:,15)
  Tsoil3     => TargetXX(:,16)
  Tsoil4     => TargetXX(:,17)
  Tsoil5     => TargetXX(:,18)
  Tsoil6     => TargetXX(:,19)
  Tsoil7     => TargetXX(:,20)
  Tsoil8     => TargetXX(:,21)
  Tsoil9     => TargetXX(:,22)
  Tsoil10     => TargetXX(:,23)
  cplant1    => TargetXX(:,24)
  cplant2    => TargetXX(:,25)
  cplant3    => TargetXX(:,26)
  csoil1     => TargetXX(:,27)
  csoil2     => TargetXX(:,28)
  cansto     => TargetXX(:,29)
  theta1	 => TargetXX(:,30)
  theta2	 => TargetXX(:,31)
  theta3	 => TargetXX(:,32)
  theta4	 => TargetXX(:,33)
  theta5	 => TargetXX(:,34)
  theta6	 => TargetXX(:,35)
  theta7	 => TargetXX(:,36)
  theta8	 => TargetXX(:,37)
  theta9	 => TargetXX(:,38)
  theta10	 => TargetXX(:,39)
  wcol		 => TargetXX(:,40)
end if
! * Pointers for fluxes and dXXdt (array FF)
if (present(TargetFF)) then
  if (size(TargetFF,2) /= 23) stop "PointAll: wrong nff"
  dWColAdt   => TargetFF(:,01)
  dWColBdt   => TargetFF(:,02)
  FWPrec     => TargetFF(:,03)
  FWTraA     => TargetFF(:,04)
  FWTraB    => TargetFF(:,05)
  FWSoil     => TargetFF(:,06)
  FWRun      => TargetFF(:,07)
  FWLchA     => TargetFF(:,08)
  FWLchB     => TargetFF(:,09)
  FWE        => TargetFF(:,10)
  FWTra      => TargetFF(:,11)
  FWDis      => TargetFF(:,12)
  FWThrough  => TargetFF(:,13)
  FWwc       => TargetFF(:,14)
  dWcanopydt => TargetFF(:,15)
  FCGPP      => TargetFF(:,16)
  FCGro      => TargetFF(:,17)
  FCNEE      => TargetFF(:,18)
  phiE       => TargetFF(:,19)
  phiH       => TargetFF(:,20)
  FWPT       => TargetFF(:,21)
  dCLeadt       => TargetFF(:,22)
  dWColdt   => TargetFF(:,23)
end if
if (present(TargetDD)) then
  if (size(TargetDD,2) /= 30) stop "PointAll: wrong ndd"
! * Pointers for derived quantities (array DD)
  DayltFrac   => TargetDD(:,01)
  RhoA        => TargetDD(:,02)
  FracVExt    => TargetDD(:,03)
  rLAIExt     => TargetDD(:,04)
  FracVCLea   => TargetDD(:,05)
  rLAICLea    => TargetDD(:,06)
  AllocLea    => TargetDD(:,07)
  FCGroL      => TargetDD(:,08)
  FCGroW      => TargetDD(:,09)
  ImBalA      => TargetDD(:,10)
  ImBalB      => TargetDD(:,11)
  ImBalCanopy => TargetDD(:,12)
  ImBalSoil   => TargetDD(:,13)
  ImBal       => TargetDD(:,14)
  fws      => TargetDD(:,15)
  Tabar     => TargetDD(:,16)   
  Tsoilbar_g    => TargetDD(:,17)
  Tsoilbar_w    => TargetDD(:,18)
  Sbar_g      => TargetDD(:,19)
  Sbar_w      => TargetDD(:,20) 
  btran_g     => TargetDD(:,21) 
  btran_w     => TargetDD(:,22)
  FCGPP_g     => TargetDD(:,23) 
  FCGPP_w     => TargetDD(:,24)  
  FCLeafR_g    => TargetDD(:,25)
  FCLeafR_w    => TargetDD(:,26)
  LAI_g => TargetDD(:,27)
  LAI_w => TargetDD(:,28)
  fws_g => TargetDD(:,29)
  fws_w => TargetDD(:,30)
end if
! * Pointers for predicted point observables (array ZZP)
if (present(TargetZZP)) then
  if (size(TargetZZP,2) /= 12) stop "PointAll: wrong nzzP"
  ZNDVI       => TargetZZP(:,01)
  ZphiRnet    => TargetZZP(:,02)
  ZphiH       => TargetZZP(:,03)
  ZphiE       => TargetZZP(:,04)
  ZphiNEE     => TargetZZP(:,05)
  ZphiNPP     => TargetZZP(:,06)
  ZphiGPP     => TargetZZP(:,07)
  Zsm0008    => TargetZZP(:,08)
  Zsm0090     => TargetZZP(:,09)
  Zsm0030     => TargetZZP(:,10)
  Zsm3060     => TargetZZP(:,11)
  Zsm6090     => TargetZZP(:,12)
end if
  ! * Pointers for predicted point observables (array AAP)
if (present(TargetAAP)) then
  if (size(TargetAAP,2) /= 12) stop "PointAll: wrong naaP"
  ANDVI       => TargetAAP(:,01)
  AphiRnet    => TargetAAP(:,02)
  AphiH       => TargetAAP(:,03)
  AphiE       => TargetAAP(:,04)
  AphiNEE     => TargetAAP(:,05)
  AphiNPP     => TargetAAP(:,06)
  AphiGPP     => TargetAAP(:,07)
  Asm0008     => TargetAAP(:,08)
  Asm0090     => TargetAAP(:,09)
  Asm0030     => TargetAAP(:,10)
  Asm3060     => TargetAAP(:,11)
  Asm6090     => TargetAAP(:,12)
  ALAIg     => TargetAAP(:,13)
  ALAIw     => TargetAAP(:,14)
  AfWoody     => TargetAAP(:,15)
  AscattVIS     => TargetAAP(:,16)
  AscattNIR     => TargetAAP(:,17)

end if
! * Pointers for predicted point observables (array ZZPh)
if (present(TargetZZPh)) then
  if (size(TargetZZPh,2) /= 35) stop "PointAll: wrong nzzPh"
  ZLST        => TargetZZPh(:,01,:)
  ZLSTTime    => TargetZZPh(:,02,:)
  ZLSTAngle   => TargetZZPh(:,03,:)
  DelTsTaZ    => TargetZZPh(:,04,:)
  TaZ         => TargetZZPh(:,05,:)
  ZphiRneth   => TargetZZPh(:,06,:)
  ZphiHh      => TargetZZPh(:,07,:)
  ZphiEh      => TargetZZPh(:,08,:)
  ZphiNEEh    => TargetZZPh(:,09,:)
  ZphiGh      => TargetZZPh(:,10,:)
  ZTsoilh      => TargetZZPh(:,11,:)
  ZphiHhTime  => TargetZZPh(:,12,:)
  TrVeg       => TargetZZPh(:,13,:)
  TrSoil      => TargetZZPh(:,14,:)
  TrEff       => TargetZZPh(:,15,:)
  TSoil1h       => TargetZZPh(:,16,:)
  TSoil2h       => TargetZZPh(:,17,:)
  TSoil3h       => TargetZZPh(:,18,:)
  TSoil4h       => TargetZZPh(:,19,:)
  TSoil5h       => TargetZZPh(:,20,:)
  TSoil6h       => TargetZZPh(:,21,:)
  TSoil7h       => TargetZZPh(:,22,:)
  TSoil8h       => TargetZZPh(:,23,:)
  TSoil9h       => TargetZZPh(:,24,:)
  TSoil10h       => TargetZZPh(:,25,:)
  theta1h       => TargetZZPh(:,26,:)
  theta2h       => TargetZZPh(:,27,:)
  theta3h       => TargetZZPh(:,28,:)
  theta4h       => TargetZZPh(:,29,:)
  theta5h       => TargetZZPh(:,30,:)
  theta6h       => TargetZZPh(:,31,:)
  theta7h       => TargetZZPh(:,32,:)
  theta8h       => TargetZZPh(:,33,:)
  theta9h       => TargetZZPh(:,34,:)
  theta10h       => TargetZZPh(:,35,:)
end if
! * Pointers for actual point observations (array AAPh)
if (present(TargetAAPh)) then
  if (size(TargetAAPh,2) /= 12) stop "PointAll: wrong naaPh"
  ALST       => TargetAAPh(:,01,:)
  ALSTTime   => TargetAAPh(:,02,:)
  ALSTAngle  => TargetAAPh(:,03,:)
  DelTsTaA   => TargetAAPh(:,04,:)
  TempAt     => TargetAAPh(:,05,:)
  AphiRneth  => TargetAAPh(:,06,:)
  AphiHh     => TargetAAPh(:,07,:)
  AphiEh     => TargetAAPh(:,08,:)
  AphiNEEh   => TargetAAPh(:,09,:)
  AphiGh     => TargetAAPh(:,10,:)
  ATsoilh     => TargetAAPh(:,11,:)
  AphiHhTime => TargetAAPh(:,12,:)
end if
! * Pointers for predicted catchment observables (array ZZC)
if (present(TargetZZC)) then
  if (size(TargetZZC,2) /= 4) stop "PointAll: wrong nzzC"
  ZDisCD     => TargetZZC(:,01)
  ZDisCM     => TargetZZC(:,02)  
  ZRunCD     => TargetZZC(:,03)
  ZLchCD     => TargetZZC(:,04)
end if
! * Pointers for actual catchment observables (array AAC)
if (present(TargetAAC)) then
  if (size(TargetAAC,2) /= 4) stop "PointAll: wrong naaC"
  ADisCD     => TargetAAC(:,01)
  ADisCM     => TargetAAC(:,02)  
  ARunCD     => TargetAAC(:,03)
  ALchCD     => TargetAAC(:,04)
end if
! * Pointers for met forcing variables (array MM)
if (present(TargetMM)) then
  if (size(TargetMM,2) /= 6) stop "PointAll: wrong nmm"
  SolarMJ    => TargetMM(:,01)
  Precip     => TargetMM(:,02)
  TempMax    => TargetMM(:,03)
  TempMin    => TargetMM(:,04)
  vph09      => TargetMM(:,05)
  vph15      => TargetMM(:,06)
end if
! * Pointers for met forcing variables (array hMM)
if (present(TargethMM)) then
  if (size(TargethMM,2) /= 8) stop "PointAll: wrong nhMM"
  hFsd	    => TargethMM(:,01,:)
  hFld      => TargethMM(:,02,:)
  hPrecip	=> TargethMM(:,03,:)
  hUa       => TargethMM(:,04,:)
  hTc	    => TargethMM(:,05,:)
  hqv	    => TargethMM(:,06,:)
  hpmb	    => TargethMM(:,07,:)
  hcoszen	=> TargethMM(:,08,:)
end if
! * Pointers for remote seinsing forcing variables (array RR)
if (present(TargetRR)) then
  if (size(TargetRR,2) /= 5) stop "PointAll: wrong nRR"
  LAIg		=> TargetRR(:,01)
  LAIw		=> TargetRR(:,02)
  fWoody    => TargetRR(:,03)
  scattVIS  => TargetRR(:,04)
  scattNIR  => TargetRR(:,05)
end if
! * Pointers for spatially uniform parameters (array UU)
if (present(TargetUU)) then
  if (size(TargetUU) /= 63) stop "PointAll: wrong nuu"
  CoeffPT    => TargetUU(01)
  CoeffBeer  => TargetUU(02)
  CLea0      => TargetUU(03)
  RateCLea   => TargetUU(04)
  AllocLg    => TargetUU(05)
  AllocLw    => TargetUU(06)
  HyConSat1  => TargetUU(07)
  HyConSat2  => TargetUU(08)
  PwrFWSoil  => TargetUU(09)
  PwrFWLch   => TargetUU(10)
  alfaQ      => TargetUU(11)
  alfaWpri   => TargetUU(12)
  alfaWmul   => TargetUU(13)
  CO2A       => TargetUU(14)
  WRelA0     => TargetUU(15)
  ZSoil1Mult => TargetUU(16)
  ZSoil2Mult => TargetUU(17)
  rLAImax    => TargetUU(18)
  Gaero      => TargetUU(19)
  TimeTxFrac => TargetUU(20)
  cN0        => TargetUU(21)
  cN1        => TargetUU(22)
  cN2        => TargetUU(23)
  TZRunDef   => TargetUU(24)
  TZLchDef   => TargetUU(25)
  CoeffPAR   => TargetUU(26)
  RatioJV    => TargetUU(27)
  za         => TargetUU(28)
  ratecp1    => TargetUU(29)
  ratecp2    => TargetUU(30)
  ratecp3    => TargetUU(31)
  ratecs1    => TargetUU(32)
  ratecs2    => TargetUU(33)
  zeta    => TargetUU(34)
  fsatmax    => TargetUU(35)
  B1Mult    => TargetUU(36)
  Psie1Mult    => TargetUU(37)
  Ksat1Mult    => TargetUU(38)
  B2Mult    => TargetUU(39)
  Psie2Mult    => TargetUU(40)
  Ksat2Mult => TargetUU(41)
  dleaf_g    => TargetUU(42)
  vcmax_g    => TargetUU(43)
  hc_g       => TargetUU(44)
  xfang_g    => TargetUU(45)
  rp20_g     => TargetUU(46)
  vbeta_g    => TargetUU(47)
  rootbeta_g => TargetUU(48)
  F10_g      => TargetUU(49)
  ZR_g       => TargetUU(50)
  loggamma_g    => TargetUU(51)
  dleaf_w    => TargetUU(52)
  vcmax_w    => TargetUU(53)
  hc_w       => TargetUU(54)
  xfang_w    => TargetUU(55)
  rp20_w     => TargetUU(56)
  vbeta_w    => TargetUU(56)
  rootbeta_w => TargetUU(58)
  F10_w      => TargetUU(59)
  ZR_w       => TargetUU(60)
  loggamma_w    => TargetUU(61)
  a1		    => TargetUU(62)
  ds0		    => TargetUU(63)
end if
! * Pointers for spatially variable parameters (array VV)
if (present(TargetVV)) then
  if (size(TargetVV,2) /= 72) stop "PointAll: wrong nvv"
  CatchMap   => TargetVV(:,01)
  Albedo     => TargetVV(:,02)
  LatDeg     => TargetVV(:,03)
  LongDeg     => TargetVV(:,04)
  Altitude   => TargetVV(:,05)
  FAPAR      => TargetVV(:,06:17)
  canst1     => TargetVV(:,18)
  dleaf		 => TargetVV(:,19)
  vcmax		 => TargetVV(:,20)
  hc		 => TargetVV(:,21)
  xfang		 => TargetVV(:,22)
  rp20		 => TargetVV(:,23)
  rpcoef	 => TargetVV(:,24)
  rs20		 => TargetVV(:,25)
  shelrb	 => TargetVV(:,26)
  frac4		 => TargetVV(:,27)
  tminvj	 => TargetVV(:,28)
  tmaxvj	 => TargetVV(:,29)
  vbeta		 => TargetVV(:,30)
  albsoil	 => TargetVV(:,31)
  swilt1	 => TargetVV(:,32)
  sfc1		 => TargetVV(:,33)
  WVolSat1   => TargetVV(:,34)
  bch1		 => TargetVV(:,35)
  hyds1		 => TargetVV(:,36)
  sucs1		 => TargetVV(:,37)
  rhosoil1	 => TargetVV(:,38)
  css1		 => TargetVV(:,39)
  ZSoil1     => TargetVV(:,40)
  clay1      => TargetVV(:,41)
  silt1      => TargetVV(:,42)
  swilt2	 => TargetVV(:,43)
  sfc2		 => TargetVV(:,44)
  WVolSat2   => TargetVV(:,45)
  bch2		 => TargetVV(:,46)
  hyds2		 => TargetVV(:,47)
  sucs2		 => TargetVV(:,48)
  rhosoil2	 => TargetVV(:,49)
  css2		 => TargetVV(:,50)
  ZSoil2     => TargetVV(:,51)
  clay2      => TargetVV(:,52)
  silt2      => TargetVV(:,53)
  rootbeta	 => TargetVV(:,54)
  extkn      => TargetVV(:,55)
  vegcf      => TargetVV(:,56)
  iveg       => TargetVV(:,57)
  LAI         => TargetVV(:,58:69)
  clittequil     => TargetVV(:,70)
  albsoilVIS     => TargetVV(:,71)
  albsoilNIR     => TargetVV(:,72)
end if

END SUBROUTINE PointAll

!*******************************************************************************

END MODULE PointerModule
