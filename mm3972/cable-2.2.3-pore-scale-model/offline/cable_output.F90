!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Output module for CABLE offline 
!
! Contact: Bernard.Pak@csiro.au
!
! History: Developed by Gab Abramowitz
!          Output of additional variables and parameters relative to v1.4b
!
!
! ==============================================================================
! CALLed from:    cable_driver.F90
!
! MODULEs used:   cable_abort_module
!                 cable_common_module
!                 cable_checks_module
!                 cable_def_types_mod
!                 cable_IO_vars_module
!                 cable_write_module
!                 netcdf
!
! CALLs:          open_output_file
!                 write_output
!                 close_output_file
!                 create_restart
!
MODULE cable_output_module


  USE cable_abort_module, ONLY: abort, nc_abort
  USE cable_def_types_mod
  USE cable_IO_vars_module
  USE cable_checks_module, ONLY: mass_balance, energy_balance, ranges
  USE cable_write_module
  USE netcdf
  USE cable_common_module, ONLY: filename,cable_user
  IMPLICIT NONE
  PRIVATE
  PUBLIC open_output_file, write_output, close_output_file, create_restart
  INTEGER :: ncid_out ! output data netcdf file ID
  REAL :: missing_value = -999999.0 ! for netcdf output
  TYPE out_varID_type ! output variable IDs in netcdf file
    INTEGER :: SWdown, LWdown, Wind, Wind_E, PSurf,                            &
                    Tair, Qair, Rainf, Snowf, CO2air,                          &
                    Qle, Qh, Qg, NEE, SWnet,                                   &
                    LWnet, SoilMoist, SoilTemp, Albedo, Qs,                    &
                    Qsb, Evap, BaresoilT, SWE, SnowT,                          &
                    RadT, VegT, Ebal, Wbal, AutoResp,                          &
                    LeafResp, HeteroResp, GPP, NPP, LAI,                       &
                    ECanop, TVeg, ESoil, CanopInt, SnowDepth,                  &
                    HVeg, HSoil, Rnet, tvar,cancd, gswx_1, gswx_2, gswmin_1,   &
                    gswmin_2,                                                  &
                    !MD
                    WatTable,GWMoist,SoilMatPot,EqSoilMatPot,EqSoilMoist,      &
                    EqGWMoist,EqGWSoilMatPot,Qinfl,GWSoilMatPot,fldcap,forg,   &
                    wiltp,SoilIce,SatFrac,Qrecharge,                           &
                    VISalbedo,NIRalbedo, &
                    rtevap,sublayer_dz,rtevap_sat,z0soil,rtsoil,beta_rtsoil,&
                    sv_rtevap,bl_rtevap
  END TYPE out_varID_type
  TYPE(out_varID_type) :: ovid ! netcdf variable IDs for output variables
  TYPE(parID_type) :: opid ! netcdf variable IDs for output variables
  TYPE output_temporary_type
    REAL(KIND=4), POINTER, DIMENSION(:) :: SWdown ! 6 downward short-wave
                                                  ! radiation [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: LWdown ! 7 downward long-wave
                                                  ! radiation [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Rainf  ! 8 rainfall [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Snowf  ! 9 snowfall [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: PSurf  ! 10 surface pressure [Pa]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Tair   ! 11 surface air temperature
                                                  ! [K]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Qair   ! 12 specific humidity [kg/kg]
    REAL(KIND=4), POINTER, DIMENSION(:) :: CO2air ! 13 CO2 concentration [ppmv]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Wind   ! 14 windspeed [m/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Wind_N ! 15 surface wind speed, N
                                                  ! component [m/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Wind_E ! 16 surface wind speed, E
                                                  ! component [m/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: LAI
    REAL(KIND=4), POINTER, DIMENSION(:) :: Qh     ! 17 sensible heat flux [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Qle    ! 18 latent heat flux [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Qg     ! 19 ground heat flux [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: SWnet  ! 20 net shortwave [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: LWnet  ! 21 net longwave [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Evap   ! 22 total evapotranspiration
                                                  ! [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Ewater ! 23 evap. from surface water
                                                  ! storage [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: ESoil  ! 24 bare soil evaporation
                                                  ! [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: TVeg   ! 25 vegetation transpiration
                                                  ! [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: ECanop ! 26 interception evaporation
                                                  ! [kg/m2/s]
    ! 27 potential evapotranspiration [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: PotEvap
    REAL(KIND=4), POINTER, DIMENSION(:) :: ACond   ! 28 aerodynamic conductance
                                                   ! [m/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: SoilWet ! 29 total soil wetness [-]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Albedo  ! 30 albedo [-]
    REAL(KIND=4), POINTER, DIMENSION(:) :: VegT    ! 31 vegetation temperature
                                                   ! [K]
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilTemp  ! 32 av.layer soil
                                                       ! temperature [K]
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilMoist ! 33 av.layer soil
                                                       ! moisture [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Qs  ! 34 surface runoff [kg/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Qsb ! 35 subsurface runoff [kg/m2/s]
    ! 36 change in soilmoisture (sum layers) [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: DelSoilMoist
    ! 37 change in snow water equivalent [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: DelSWE
    ! 38 change in interception storage [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: DelIntercept
    REAL(KIND=4), POINTER, DIMENSION(:) :: SnowT     ! 39 snow surface temp [K]
    REAL(KIND=4), POINTER, DIMENSION(:) :: BaresoilT ! 40 surface bare soil
                                                     ! temp [K]
    REAL(KIND=4), POINTER, DIMENSION(:) :: AvgSurfT  ! 41 Average surface
                                                     ! temperature [K]
    REAL(KIND=4), POINTER, DIMENSION(:) :: RadT      ! 42 Radiative surface
                                                     ! temperature [K]
    REAL(KIND=4), POINTER, DIMENSION(:) :: SWE       ! 43 snow water equivalent
                                                     ! [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: RootMoist ! 44 root zone soil
                                                     ! moisture [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: CanopInt  ! 45 total canopy water
                                                     ! storage [kg/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: NEE       ! 46 net ecosystem exchange
                                                     ! [umol/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: NPP       ! 47 net primary production
                                                     ! of C by veg [umol/m2/s]
    ! 48 gross primary production C by veg [umol/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: GPP
    REAL(KIND=4), POINTER, DIMENSION(:) :: AutoResp   ! 49 autotrophic
                                                      ! respiration [umol/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: LeafResp   ! 51 autotrophic
                                                      ! respiration [umol/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: HeteroResp ! 50 heterotrophic
                                                      ! respiration [umol/m2/s]
    REAL(KIND=4), POINTER, DIMENSION(:) :: SnowDepth  ! actual depth of snow in
                                                      ! [m]
    ! Non-Alma variables
    REAL(KIND=4), POINTER, DIMENSION(:) :: Rnet  ! net absorbed radiation [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: HVeg  ! sensible heat from vegetation
                                                 ! [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: HSoil ! sensible heat from soil
                                                 ! [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Ebal  ! cumulative energy balance
                                                 ! [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: Wbal  ! cumulative water balance
                                                 ! [W/m2]
    REAL(KIND=4), POINTER, DIMENSION(:) :: cancd ! Canopy cond (m/s) jtk561
    REAL(KIND=4), POINTER, DIMENSION(:) :: gswx_1 ! sunlit cond (dunnounits),jtk561
    REAL(KIND=4), POINTER, DIMENSION(:) :: gswx_2 ! shaded cond (dunno units),jtk561
    REAL(KIND=4), POINTER, DIMENSION(:) :: gswmin_1 ! min sunlit cond, jtk561
    REAL(KIND=4), POINTER, DIMENSION(:) :: gswmin_2 ! min shaded cond, jtk561

    !MD GW
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilMatPot    ! soil matric potential [mm]
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: EqSoilMatPot  ! equilibirum soil matric potential [mm]
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: EqSoilMoist   ! equilibirum soil moisture [mm3/mm3]    
    REAL(KIND=4), POINTER, DIMENSION(:)   :: GWMoist       ! water balance of aquifer [mm3/mm3]
    REAL(KIND=4), POINTER, DIMENSION(:)   :: WatTable      ! water table depth [m]
    REAL(KIND=4), POINTER, DIMENSION(:)   :: EqGWMoist     ! equilibrium water of aquifer
    REAL(KIND=4), POINTER, DIMENSION(:)   :: EqGWSoilMatPot    ! equilibrium soil matric potential of aquifer [mm3/mm3]
    REAL(KIND=4), POINTER, DIMENSION(:)   :: GWSoilMatPot    ! equilibrium soil matric potential of aquifer [mm3/mm3]     
    REAL(KIND=4), POINTER, DIMENSION(:)   :: Qinfl         !infiltration rate into first soil layer [mm/s] 
    REAL(KIND=4), POINTER, DIMENSION(:)   :: SatFrac         !Saturated Fraction of Grid Cell
    REAL(KIND=4), POINTER, DIMENSION(:)   :: Qrecharge         !Saturated Fraction of Grid Cell

    REAL(KIND=4), POINTER, DIMENSION(:,:) :: wiltp         !wilt pnt inc forg
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: fldcap        !field capcaicty adj for organic content
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: Forg          !organic carbon frac.soil
    REAL(KIND=4), POINTER, DIMENSION(:,:) :: SoilIce       !SOil Ice volume [mm3/mm3]

    REAL(KIND=4), POINTER, DIMENSION(:) :: VISalbedo
    REAL(KIND=4), POINTER, DIMENSION(:) :: NIRalbedo

    REAL(KIND=4), POINTER, DIMENSION(:) :: rtevap,rtevap_sat,sublayer_dz,rtsoil,z0soil,beta_rtsoil,&
                                           sv_rtevap,bl_rtevap


  END TYPE output_temporary_type
  TYPE(output_temporary_type), SAVE :: out
  INTEGER :: ok   ! netcdf error status

CONTAINS

  SUBROUTINE open_output_file(dels, soil, veg, bgc, rough)
    ! Creates netcdf output file, defines all variables 
    ! and writes parameters to it if requested by user.
    REAL, INTENT(IN) :: dels ! time step size
    TYPE (soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)  :: veg  ! vegetation parameters
    TYPE (bgc_pool_type), INTENT(IN)       :: bgc
    TYPE (roughness_type), INTENT(IN)      :: rough
    ! REAL, POINTER,DIMENSION(:,:) :: surffrac ! fraction of each surf type

    INTEGER :: xID, yID, zID, radID, soilID, soilcarbID,                  &
                    plantcarbID, tID, landID, patchID ! dimension IDs
    INTEGER :: latID, lonID ! time,lat,lon variable ID
    INTEGER :: xvID, yvID   ! coordinate variable IDs for GrADS readability
    !    INTEGER :: surffracID         ! surface fraction varaible ID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
  
    ! Create output file:
    ok = NF90_CREATE(filename%out, NF90_CLOBBER, ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating output file '       &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ! Put the file in define mode:
    ok = NF90_REDEF(ncid_out)
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_out, 'x', xdimsize, xID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                          (ok, 'Error defining x dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'y', ydimsize, yID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                          (ok, 'Error defining y dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ! Define patch dimension, whether it's used or not:
    ok = NF90_DEF_DIM(ncid_out, 'patch', max_vegpatches, patchID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                       (ok,'Error defining patch dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    !   ! Define surftype dimension (currently only used for surffrac variable):
    !    ok = NF90_DEF_DIM(ncid_out,'surftype',4,surftypeID) 
    !    IF (ok /= NF90_NOERR) CALL nc_abort &
    !         (ok,'Error defining syrftype dimension in output file. '// &
    !         '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
              (ok, 'Error defining vertical soil dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'rad', nrb, radID) ! number of radiation bands
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                  (ok, 'Error defining radiation dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'soil_carbon_pools', ncs, soilcarbID) ! # pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
           (ok, 'Error defining soil carbon pool dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out,'plant_carbon_pools',ncp,plantcarbID) 
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
           (ok,'Error defining plant carbon pool dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_DIM(ncid_out, 'time', NF90_UNLIMITED, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                        (ok,'Error defining time dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    IF(output%grid == 'mask' .OR. output%grid == 'ALMA' .OR.                   &
       (metGrid == 'mask' .AND. output%grid == 'default')) THEN
       ! for land/sea mask type grid:
       ! Atmospheric 'z' dim of size 1 to comply with ALMA grid type:
       ok = NF90_DEF_DIM(ncid_out, 'z', 1, zID)
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
                          (ok, 'Error defining z dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ELSE IF(output%grid == 'land' .OR.                                         &
            (metGrid == 'land' .AND. output%grid == 'default')) THEN
       ! For land only compression grid:
       ok = NF90_DEF_DIM(ncid_out, 'land', mland, landID) ! number of land
                                                          ! points
       IF (ok /= NF90_NOERR) CALL nc_abort                                     &
                       (ok, 'Error defining land dimension in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    END IF
    ! Define "time" variable and its attributes:
    ok = NF90_DEF_VAR(ncid_out, 'time', NF90_DOUBLE, (/tID/), ovid%tvar)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                        (ok, 'Error defining time variable in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
             (ok, 'Error defining time variable attributes in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, ovid%tvar, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
             (ok, 'Error defining time variable attributes in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ! Define latitude and longitude variable (ALMA):
    ok = NF90_DEF_VAR(ncid_out, 'latitude', NF90_FLOAT, (/xID, yID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                    (ok, 'Error defining latitude variable in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, latID, 'units', 'degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining latitude variable attributes in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_VAR(ncid_out, 'longitude', NF90_FLOAT, (/xID, yID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                   (ok, 'Error defining longitude variable in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
        (ok, 'Error defining longitude variable attributes in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ! Write "cordinate variables" to enable reading by GrADS:
    ok = NF90_DEF_VAR(ncid_out, 'x', NF90_FLOAT, (/xID/), xvID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
             (ok, 'Error defining "x" variable (for GrADS) in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, xvID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error writing x coordinate variable (GrADS) units in output '// &
                                          'file. (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, xvID, 'comment',                               &
                                'x coordinate variable for GrADS compatibility')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                   (ok, 'Error writing x variables comment in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_DEF_VAR(ncid_out, 'y', NF90_FLOAT, (/yID/), yvID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
             (ok, 'Error defining "y" variable (for GrADS) in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, yvID, 'units', 'degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
        (ok, 'Error writing y coordinate variable (GrADS) units in output '//  &
                                          'file. (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, yvID, 'comment',                               &
                                'y coordinate variable for GrADS compatibility')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                   (ok, 'Error writing y variables comment in output file. '// &
                                                '(SUBROUTINE open_output_file)')
    !   ! Define fraction of each surface type:
    !   CALL define_ovar(ncid_out,surffracID,'surffrac','-', &
    !       'Fraction of each surface type: vegetated; urban; lake; land ice', &
    !       .FALSE.,surftypeID,'surftype',xID,yID,zID,landID,patchID)

    !=============DEFINE OUTPUT VARIABLES=======================================
    ! Define met forcing variables in output file and allocate temp output vars:
    IF(output%met .OR. output%SWdown) THEN
       CALL define_ovar(ncid_out,                                              &
            ovid%SWdown, 'SWdown', 'W/m^2', 'Downward shortwave radiation',    &
            patchout%SWdown, 'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SWdown(mp))
       out%SWdown = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%LWdown) THEN
       CALL define_ovar(ncid_out, ovid%LWdown, 'LWdown', 'W/m^2',              &
                        'Downward longwave radiation', patchout%LWdown,        &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LWdown(mp))
       out%LWdown = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%Tair) THEN
       CALL define_ovar(ncid_out, ovid%Tair,                                   &
                        'Tair', 'K', 'Surface air temperature', patchout%Tair, &
                        'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Tair(mp))
       out%Tair = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%Rainf) THEN
       CALL define_ovar(ncid_out, ovid%Rainf, 'Rainf',                         &
                        'kg/m^2/s', 'Rainfall+snowfall', patchout%Rainf,       &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Rainf(mp))
       out%Rainf = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%Snowf) THEN
       CALL define_ovar(ncid_out, ovid%Snowf, 'Snowf',                         &
                        'kg/m^2/s', 'Snowfall', patchout%Snowf,                &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Snowf(mp))
       out%Snowf = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%Qair) THEN
       CALL define_ovar(ncid_out, ovid%Qair, 'Qair',                           &
                        'kg/kg', 'Surface specific humidity', patchout%Qair,   &
                        'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qair(mp))
       out%Qair = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%Wind) THEN
       CALL define_ovar(ncid_out, ovid%Wind, 'Wind',                           &
                        'm/s', 'Scalar surface wind speed', patchout%Wind,     &
                        'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Wind(mp))
       out%Wind = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%PSurf) THEN
       CALL define_ovar(ncid_out, ovid%PSurf, 'PSurf',                         &
                        'hPa', 'Surface air pressure', patchout%PSurf,         &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%PSurf(mp))
       out%PSurf = 0.0 ! initialise
    END IF
    IF(output%met .OR. output%CO2air) THEN
       CALL define_ovar(ncid_out, ovid%CO2air, 'CO2air', 'ppmv',               &
                        'Surface air CO2 concentration', patchout%CO2air,      &
                        'ALMA', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CO2air(mp))
       out%CO2air = 0.0 ! initialise
    END IF
    ! Define surface flux variables in output file and allocate temp output
    ! vars:
    IF(output%flux .OR. output%Qle) THEN
       CALL define_ovar(ncid_out, ovid%Qle, 'Qle', 'W/m^2',                    &
                        'Surface latent heat flux',patchout%Qle,'dummy',       &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qle(mp))
       out%Qle = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%Qh) THEN
       CALL define_ovar(ncid_out,ovid%Qh,'Qh', 'W/m^2',                        &
                        'Surface sensible heat flux',patchout%Qh,'dummy',      &
                        xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%Qh(mp))
       out%Qh = 0.0 ! initialise
    END IF
    ! jtk561
    IF(output%flux .OR. output%cancd) THEN 
       CALL define_ovar(ncid_out,ovid%cancd,'cancd', 'm/s^1',&
                        'Canopy Conductance',patchout%cancd,'dummy',      &    
                        xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%cancd(mp))
       out%cancd = 0.0 ! initialise
    END IF
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswx_1)) THEN   !gswx_1 no in
                                                                   !mpi code!
       CALL define_ovar(ncid_out,ovid%gswx_1,'gswx_1', 'NA',&
                        'Sunlit Conductance',patchout%gswx_1,'dummy',      &    
                        xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%gswx_1(mp))
       out%gswx_1 = 0.0 ! initialise
    END IF
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswx_2)) THEN 
       CALL define_ovar(ncid_out,ovid%gswx_2,'gswx_2', 'NA',&
                        'Shaded Conductance',patchout%gswx_2,'dummy',      &    
                        xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%gswx_2(mp))
       out%gswx_2 = 0.0 ! initialise
    END IF
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswmin_1)) THEN 
       CALL define_ovar(ncid_out,ovid%gswmin_1,'gswmin_1', 'NA',&
                        'Min Sunlit Conductance',patchout%gswmin_1,'dummy',&    
                        xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%gswmin_1(mp))
       out%gswmin_1 = 0.0 ! initialise
    END IF
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswmin_2)) THEN 
       CALL define_ovar(ncid_out,ovid%gswmin_2,'gswmin_2', 'NA',&
                        'Min Shaded Conductance',patchout%gswmin_2,'dummy',&    
                        xID,yID,zID,landID,patchID,tID)
       ALLOCATE(out%gswmin_2(mp))
       out%gswmin_2 = 0.0 ! initialise
    END IF
    ! end modifs jtk561

    IF(output%flux .OR. output%Qg) THEN
       CALL define_ovar(ncid_out, ovid%Qg, 'Qg', 'W/m^2',                      &
                        'Surface ground heat flux', patchout%Qg, 'dummy',      &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qg(mp))
       out%Qg = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%Qs) THEN
       CALL define_ovar(ncid_out, ovid%Qs, 'Qs',                               &
                        'kg/m^2/s', 'Surface runoff', patchout%Qs, 'dummy',    &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qs(mp))
       out%Qs = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%Qsb) THEN
       CALL define_ovar(ncid_out, ovid%Qsb, 'Qsb', 'kg/m^2/s',                 &
                        'Subsurface runoff', patchout%Qsb, 'dummy',            &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qsb(mp))
       out%Qsb = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%Evap) THEN
       CALL define_ovar(ncid_out, ovid%Evap,'Evap', 'kg/m^2/s',                &
                        'Total evapotranspiration', patchout%Evap, 'dummy',    &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Evap(mp))
       out%Evap = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%ECanop) THEN
       CALL define_ovar(ncid_out, ovid%Ecanop, 'ECanop', 'kg/m^2/s',           &
                        'Wet canopy evaporation', patchout%ECanop, 'dummy',    &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ECanop(mp))
       out%ECanop = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%TVeg) THEN
       CALL define_ovar(ncid_out, ovid%TVeg, 'TVeg', 'kg/m^2/s',               &
                        'Vegetation transpiration', patchout%TVeg, 'dummy',    &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%TVeg(mp))
       out%TVeg = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%ESoil) THEN
       CALL define_ovar(ncid_out, ovid%ESoil, 'ESoil', 'kg/m^2/s',             &
                        'Evaporation from soil', patchout%ESoil, 'dummy',      &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%ESoil(mp))
       out%ESoil = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%HVeg) THEN
       CALL define_ovar(ncid_out, ovid%HVeg, 'HVeg', 'W/m^2',                  &
                        'Sensible heat from vegetation', patchout%HVeg,        &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HVeg(mp))
       out%HVeg = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%HSoil) THEN
       CALL define_ovar(ncid_out, ovid%HSoil, 'HSoil', 'W/m^2',                &
                        'Sensible heat from soil', patchout%HSoil, 'dummy',    &
                        xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HSoil(mp))
       out%HSoil = 0.0 ! initialise
    END IF
    IF(output%flux .OR. output%carbon .OR. output%NEE) THEN
       CALL define_ovar(ncid_out, ovid%NEE, 'NEE', 'umol/m^2/s',               &
                        'Net ecosystem exchange of CO2', patchout%NEE,         &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NEE(mp))
       out%NEE = 0.0 ! initialise
    END IF
    ! Define soil state variables in output file and allocate temp output vars:
    IF(output%soil .OR. output%SoilMoist) THEN
       CALL define_ovar(ncid_out, ovid%SoilMoist, 'SoilMoist', 'm^3/m^3',      &
                        'Average layer soil moisture', patchout%SoilMoist,     &
                        'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilMoist(mp,ms))
       out%SoilMoist = 0.0 ! initialise
    END IF
    IF(output%soil .OR. output%SoilTemp) THEN
       CALL define_ovar(ncid_out, ovid%SoilTemp, 'SoilTemp', 'K',              &
                        'Average layer soil temperature', patchout%SoilTemp,   &
                        'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilTemp(mp,ms))
       out%SoilTemp = 0.0 ! initialise
    END IF
    IF(output%soil .OR. output%BaresoilT) THEN
       CALL define_ovar(ncid_out, ovid%BaresoilT, 'BaresoilT',                 &
                        'K', 'Bare soil temperature', patchout%BaresoilT,      &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%BaresoilT(mp))
       out%BaresoilT = 0.0 ! initialise
    END IF
    ! Define snow state variables in output file and allocate temp output vars:
    IF(output%snow .OR. output%SWE) THEN
       CALL define_ovar(ncid_out, ovid%SWE, 'SWE', 'kg/m^2',                   &
                        'Snow water equivalent', patchout%SWE,                 &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SWE(mp))
       out%SWE = 0.0 ! initialise
    END IF
    IF(output%snow .OR. output%SnowT) THEN
       CALL define_ovar(ncid_out, ovid%SnowT, 'SnowT', 'K',                    &
                        'Snow surface temperature', patchout%SnowT,            &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowT(mp))
       out%SnowT = 0.0 ! initialise
    END IF
    IF(output%snow .OR. output%SnowDepth) THEN
       CALL define_ovar(ncid_out, ovid%SnowDepth, 'SnowDepth',                 &
                        'm', 'Snow depth', patchout%SnowDepth,                 &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SnowDepth(mp))
       out%SnowDepth = 0.0 ! initialise
    END IF
    ! Define radiative variables in output file and allocate temp output vars:
    IF(output%radiation .OR. output%SWnet) THEN
       CALL define_ovar(ncid_out, ovid%SWnet, 'SWnet', 'W/m^2',                &
                        'Net shortwave radiation absorbed by surface',         &
                         patchout%SWnet, 'dummy', xID, yID, zID, landID,       &
                         patchID, tID)
       ALLOCATE(out%SWnet(mp))
       out%SWnet = 0.0 ! initialise
    END IF
    IF(output%radiation .OR. output%LWnet) THEN
       CALL define_ovar(ncid_out, ovid%LWnet, 'LWnet', 'W/m^2',                &
                        'Net longwave radiation absorbed by surface',          &
                        patchout%LWnet, 'dummy', xID, yID, zID, landID,        &
                        patchID, tID)
       ALLOCATE(out%LWnet(mp))
       out%LWnet = 0.0 ! initialise
    END IF
    IF(output%radiation .OR. output%Rnet) THEN
       CALL define_ovar(ncid_out, ovid%Rnet, 'Rnet', 'W/m^2',                  &
                        'Net radiation absorbed by surface', patchout%Rnet,    &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Rnet(mp))
       out%Rnet = 0.0 ! initialise
    END IF
    IF(output%radiation .OR. output%Albedo) THEN
       CALL define_ovar(ncid_out, ovid%Albedo, 'Albedo', '-',                  &
                        'Surface albedo', patchout%Albedo,                     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Albedo(mp))
       out%Albedo = 0.0 ! initialise
    END IF
    IF(output%radiation .OR. output%RadT) THEN
       CALL define_ovar(ncid_out, ovid%RadT, 'RadT', 'K',                      &
                        'Radiative surface temperature', patchout%RadT,        &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%RadT(mp))
       out%RadT = 0.0 ! initialise
    END IF
    ! Define vegetation variables in output file and allocate temp output vars:
    IF(output%veg .OR. output%VegT) THEN
       CALL define_ovar(ncid_out, ovid%VegT, 'VegT', 'K',                      &
                        'Average vegetation temperature', patchout%VegT,       &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%VegT(mp))
       out%VegT = 0.0 ! initialise
    END IF
    IF(output%veg .OR. output%CanopInt) THEN
       CALL define_ovar(ncid_out, ovid%CanopInt, 'CanopInt', 'kg/m^2',         &
                        'Canopy intercepted water storage', patchout%CanopInt, &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%CanopInt(mp))
       out%CanopInt = 0.0 ! initialise
    END IF
    IF(output%veg .OR. output%LAI) THEN
       CALL define_ovar(ncid_out, ovid%LAI, 'LAI', '-',                        &
                        'Leaf area index', patchout%LAI, 'dummy', xID,         &
                        yID, zID, landID, patchID, tID)
       ALLOCATE(out%LAI(mp))
       out%LAI = 0.0 ! initialise
    END IF
    ! Define balance variables in output file and allocate temp output vars:
    IF(output%balances .OR. output%Ebal) THEN
       CALL define_ovar(ncid_out, ovid%Ebal, 'Ebal', 'W/m^2',                  &
                        'Cumulative energy imbalance', patchout%Ebal,          &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Ebal(mp))
       out%Ebal = 0.0 ! initialise
    END IF
    IF(output%balances .OR. output%Wbal) THEN
       CALL define_ovar(ncid_out, ovid%Wbal, 'Wbal', 'kg/m^2',                 &
                        'Cumulative water imbalance', patchout%Wbal,           &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Wbal(mp))
       out%Wbal = 0.0 ! initialise
    END IF
    ! Define carbon variables in output file and allocate temp output vars:
    IF(output%carbon .OR. output%AutoResp) THEN
       CALL define_ovar(ncid_out, ovid%AutoResp, 'AutoResp', 'umol/m^2/s',     &
                        'Autotrophic respiration', patchout%AutoResp,          &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%AutoResp(mp))
       out%AutoResp = 0.0 ! initialise
    END IF
    IF(output%carbon .OR. output%LeafResp) THEN
       CALL define_ovar(ncid_out, ovid%LeafResp, 'LeafResp', 'umol/m^2/s',     &
                        'Leaf respiration', patchout%LeafResp,                 &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%LeafResp(mp))
       out%LeafResp = 0.0 ! initialise
    END IF
    IF(output%carbon .OR. output%HeteroResp) THEN
       CALL define_ovar(ncid_out, ovid%HeteroResp, 'HeteroResp', 'umol/m^2/s', &
                        'Heterotrophic respiration', patchout%HeteroResp,      &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%HeteroResp(mp))
       out%HeteroResp = 0.0 ! initialise
    END IF
    IF(output%carbon.OR.output%GPP) THEN
       CALL define_ovar(ncid_out, ovid%GPP, 'GPP', 'umol/m^2/s',               &
                        'Gross primary production', patchout%GPP,              &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GPP(mp))
       out%GPP = 0.0 ! initialise
    END IF
    IF(output%carbon .OR. output%NPP) THEN
       CALL define_ovar(ncid_out, ovid%NPP, 'NPP', 'umol/m^2/s',               &
                        'Net primary production', patchout%NPP,                &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NPP(mp))
       out%NPP = 0.0 ! initialise
    END IF

    !MD groundwater related variables
    IF(output%soil .OR. output%WatTable) THEN
       CALL define_ovar(ncid_out, ovid%WatTable, 'WatTable', 'm',      &
                        'Water Table Depth', patchout%WatTable,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%WatTable(mp))
       out%WatTable = 0.0 ! initialise
    END IF
    IF(output%soil .OR. output%GWMoist) THEN
       CALL define_ovar(ncid_out, ovid%GWMoist, 'GWMoist', 'mm3/mm3',      &
                        'Aquifer mositure content', patchout%GWMoist,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GWMoist(mp))
       out%GWMoist = 0.0 ! initialise
    END IF
    IF(output%soil .OR. output%EqGWMoist) THEN
       CALL define_ovar(ncid_out, ovid%EqGWMoist, 'EqGWMoist', 'mm3/mm3',      &
                        'Aquifer Equilibrium moisture content', patchout%EqGWMoist,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%EqGWMoist(mp))
       out%EqGWMoist = 0.0 ! initialise
    END IF    
    IF(output%soil .OR. output%EqGWSoilMatPot) THEN
       CALL define_ovar(ncid_out, ovid%EqGWSoilMatPot, 'EqGWSoilMatPot', 'm',      &
                        'Aquifer equilibiruim matric potential', patchout%EqGWSoilMatPot,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%EqGWSoilMatPot(mp))
       out%EqGWSoilMatPot = 0.0 ! initialise
    END IF      
    IF(output%soil .OR. output%GWSoilMatPot) THEN
       CALL define_ovar(ncid_out, ovid%GWSoilMatPot, 'GWSoilMatPot', 'm',      &
                        'Aquifer  pressure (matric potential)', patchout%GWSoilMatPot,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%GWSoilMatPot(mp))
       out%GWSoilMatPot = 0.0 ! initialise
    END IF         
    IF(output%soil .OR. output%SoilMatPot) THEN
       CALL define_ovar(ncid_out, ovid%SoilMatPot, 'SoilMatPot', 'm',      &
                        'Average layer soil matric potential', patchout%SoilMatPot,     &
                        'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilMatPot(mp,ms))
       out%SoilMatPot = 0.0 ! initialise
    END IF    
    IF(output%soil .OR. output%EqSoilMatPot) THEN
       CALL define_ovar(ncid_out, ovid%EqSoilMatPot, 'EqSoilMatPot', 'm',      &
                        'Average layer soil equilibiruim matric potential', patchout%EqSoilMatPot,     &
                        'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%EqSoilMatPot(mp,ms))
       out%EqSoilMatPot = 0.0 ! initialise
    END IF    
    IF(output%soil .OR. output%EqSoilMoist) THEN
       CALL define_ovar(ncid_out, ovid%EqSoilMoist, 'EqSoilMoist', 'mm3/mm3',      &
                        'Average layer eq soil moisture', patchout%EqSoilMoist,     &
                        'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%EqSoilMoist(mp,ms))
       out%EqSoilMoist = 0.0 ! initialise
    END IF     

    IF(output%soil .OR. output%SatFrac) THEN
       CALL define_ovar(ncid_out, ovid%SatFrac, 'SatFrac', 'unitless',      &
                        'Saturated Fraction of Gridcell', patchout%SatFrac,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%SatFrac(mp))
       out%SatFrac = 0.0 ! initialise
    END IF         

    IF(output%soil .OR. output%Qrecharge) THEN
       CALL define_ovar(ncid_out, ovid%Qrecharge, 'Qrecharge', 'mm/s',      &
                           'Recharge to or from Aquifer', patchout%Qrecharge,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qrecharge(mp))
       out%Qrecharge = 0.0 ! initialise
    END IF   
    
    IF(output%soil .OR. output%Qinfl) THEN
       CALL define_ovar(ncid_out, ovid%Qinfl, 'Qinfl', 'mm/s',      &
                        'infiltration rate into first soil layer', patchout%Qinfl,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%Qinfl(mp))
       out%Qinfl = 0.0 ! initialise
    END IF         

    IF(output%soil .OR. output%SoilIce) THEN
       CALL define_ovar(ncid_out, ovid%SoilIce, 'SoilIce', 'm^3/m^3',      &
                        'Average layer soil ice', patchout%SoilIce,     &
                        'soil', xID, yID, zID, landID, patchID, soilID, tID)
       ALLOCATE(out%SoilIce(mp,ms))
       out%SoilIce = 0.0 ! initialise
    END IF

    IF(output%radiation .OR. output%VISalbedo) THEN
       CALL define_ovar(ncid_out, ovid%VISalbedo, 'VISalbedo', '-',                  &
                        'Surface VIS albedo', patchout%VISalbedo,                     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%VISalbedo(mp))
       out%VISalbedo = 0.0 ! initialise
    END IF    

    IF(output%radiation .OR. output%NIRalbedo) THEN
       CALL define_ovar(ncid_out, ovid%NIRalbedo, 'NIRalbedo', '-',                  &
                        'Surface NIR albedo', patchout%NIRalbedo,                     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%NIRalbedo(mp))
       out%NIRalbedo = 0.0 ! initialise
    END IF   

    IF(output%rtevap) THEN
       CALL define_ovar(ncid_out, ovid%rtevap, 'rtevap', 's/m',      &
                        'Evap resistance', patchout%rtevap,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%rtevap(mp))
       out%rtevap = 0.0 ! initialise
    END IF
    IF(output%rtevap) THEN
       CALL define_ovar(ncid_out, ovid%sv_rtevap, 'sv_rtevap', 's/m',      &
                        'Evap resistance', patchout%sv_rtevap,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%sv_rtevap(mp))
       out%sv_rtevap = 0.0 ! initialise
    END IF
    IF(output%rtevap) THEN
       CALL define_ovar(ncid_out, ovid%bl_rtevap, 'bl_rtevap', 's/m',      &
                        'Evap resistance', patchout%bl_rtevap,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%bl_rtevap(mp))
       out%bl_rtevap = 0.0 ! initialise
    END IF

    IF(output%rtevap_sat) THEN
       CALL define_ovar(ncid_out, ovid%rtevap_sat, 'rtevap_sat', 's/m',      &
                        'Evap resistance sat', patchout%rtevap_sat,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%rtevap_sat(mp))
       out%rtevap_sat = 0.0 ! initialise
    END IF


    IF(output%sublayer_dz) THEN
       CALL define_ovar(ncid_out, ovid%sublayer_dz, 'sublayer_dz', 'm',      &
                        'Sublayer thickness', patchout%sublayer_dz,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%sublayer_dz(mp))
       out%sublayer_dz = 0.0 ! initialise
    END IF

    IF(output%z0soil) THEN
       CALL define_ovar(ncid_out, ovid%z0soil, 'z0soil', 'm',      &
                        'Soil surface resistance', patchout%z0soil,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%z0soil(mp))
       out%z0soil = 0.0 ! initialise
    END IF


    IF(output%rtsoil) THEN
       CALL define_ovar(ncid_out, ovid%rtsoil, 'rtsoil', 'm',      &
                        'soil resistance', patchout%rtsoil,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%rtsoil(mp))
       out%rtsoil = 0.0 ! initialise
    END IF

    IF(output%beta_rtsoil) THEN
       CALL define_ovar(ncid_out, ovid%beta_rtsoil, 'beta_rtsoil', 'm',      &
                        'soil resistance / beta', patchout%beta_rtsoil,     &
                        'dummy', xID, yID, zID, landID, patchID, tID)
       ALLOCATE(out%beta_rtsoil(mp))
       out%beta_rtsoil = 0.0 ! initialise
    END IF
    ! Define CABLE parameters in output file:
    IF(output%params .OR. output%iveg) CALL define_ovar(ncid_out, opid%iveg,   &
                     'iveg', '-', 'Vegetation type', patchout%iveg, 'integer', &
                                                 xID, yID, zID, landID, patchID)

    IF((output%params .OR. output%patchfrac)                                   &
         .AND. (patchout%patchfrac .OR. output%patch))                         &
         CALL define_ovar(ncid_out, opid%patchfrac, 'patchfrac', '-',          &
         'Fractional cover of vegetation patches', patchout%patchfrac, 'real', &
                          xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%isoil) CALL define_ovar(ncid_out, opid%isoil, &
                         'isoil', '-', 'Soil type', patchout%isoil, 'integer', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%bch) CALL define_ovar(ncid_out, opid%bch,     &
           'bch', '-', 'Parameter b, Campbell eqn 1985', patchout%bch, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%clay) CALL define_ovar(ncid_out, opid%clay,   &
         'clay', '-', 'Fraction of soil which is clay', patchout%clay, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%sand) CALL define_ovar(ncid_out, opid%sand,   &
         'sand', '-', 'Fraction of soil which is sand', patchout%sand, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%silt) CALL define_ovar(ncid_out, opid%silt,   &
         'silt', '-', 'Fraction of soil which is silt', patchout%silt, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%ssat) CALL define_ovar(ncid_out, opid%ssat,   &
           'ssat', '-', 'Fraction of soil volume which is water @ saturation', &
                          patchout%ssat, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%sfc) CALL define_ovar(ncid_out, opid%sfc,     &
        'sfc', '-', 'Fraction of soil volume which is water @ field capacity', &
                           patchout%sfc, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%swilt) CALL define_ovar(ncid_out, opid%swilt, &
       'swilt', '-', 'Fraction of soil volume which is water @ wilting point', &
                         patchout%swilt, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%hyds) CALL define_ovar(ncid_out, opid%hyds,   &
                         'hyds', 'm/s', 'Hydraulic conductivity @ saturation', &
                          patchout%hyds, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%sucs) CALL define_ovar(ncid_out, opid%sucs,   &
                                          'sucs', 'm', 'Suction @ saturation', &
                          patchout%sucs, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%css) CALL define_ovar(ncid_out, opid%css,     &
                            'css', 'J/kg/C', 'Heat capacity of soil minerals', &
                           patchout%css, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%rhosoil) CALL define_ovar(ncid_out,           &
                opid%rhosoil, 'rhosoil', 'kg/m^3', 'Density of soil minerals', &
                       patchout%rhosoil, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%rs20) CALL define_ovar(ncid_out, opid%rs20,   &
                           'rs20', '-', 'Soil respiration coefficient at 20C', &
                          patchout%rs20, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%albsoil) CALL define_ovar(ncid_out,           &
                                                 opid%albsoil, 'albsoil', '-', &
                              'Snow free shortwave soil reflectance fraction', &
           patchout%albsoil, radID, 'radiation', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%hc) CALL define_ovar(ncid_out, opid%hc,       &
                                   'hc', 'm', 'Height of canopy', patchout%hc, &
                                         'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%canst1) CALL define_ovar(ncid_out,            &
           opid%canst1, 'canst1', 'mm/LAI', 'Max water intercepted by canopy', &
                        patchout%canst1, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%dleaf) CALL define_ovar(ncid_out, opid%dleaf, &
                              'dleaf', 'm', 'Chararacteristic length of leaf', &
                         patchout%dleaf, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%frac4) CALL define_ovar(ncid_out, opid%frac4, &
                              'frac4', '-', 'Fraction of plants which are C4', &
                         patchout%frac4, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%ejmax) CALL define_ovar(ncid_out, opid%ejmax, &
       'ejmax', 'mol/m^2/s', 'Max potential electron transport rate top leaf', &
                         patchout%ejmax, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%vcmax) CALL define_ovar(ncid_out, opid%vcmax, &
             'vcmax', 'mol/m^2/s', 'Maximum RuBP carboxylation rate top leaf', &
                         patchout%vcmax, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%rp20) CALL define_ovar(ncid_out, opid%rp20,   &
                          'rp20', '-', 'Plant respiration coefficient at 20C', &
                          patchout%rp20, 'real', xID, yID, zID, landID, patchID)
    ! Ticket #56
    IF(output%params .OR. output%g0c3) CALL define_ovar(ncid_out, opid%g0c3,   &
                          'g0c3', '-', 'g0c3 term in Medlyn Stom Cond. Param', &
                          patchout%g0c3, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%g0c4) CALL define_ovar(ncid_out, opid%g0c4,   &
                          'g0c4', '-', 'g0c4 term in Medlyn Stom Cond. Param', &
                          patchout%g0c4, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%g1c3) CALL define_ovar(ncid_out, opid%g1c3,   &
                          'g1c3', '-', 'g1c3 term in Medlyn Stom Cond. Param', &
                          patchout%g1c3, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%g1c4) CALL define_ovar(ncid_out, opid%g1c4,   &
                          'g1c4', '-', 'g1c4 term in Medlyn Stom Cond. Param', &
                          patchout%g1c4, 'real', xID, yID, zID, landID, patchID)
    ! end Ticket #56 

    IF(output%params .OR. output%rpcoef) CALL define_ovar(ncid_out,            &
                                                 opid%rpcoef, 'rpcoef', '1/C', &
                                 'Temperature coef nonleaf plant respiration', &
                        patchout%rpcoef, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%shelrb) CALL define_ovar(ncid_out,            &
             opid%shelrb, 'shelrb', '-', 'Sheltering factor', patchout%shelrb, &
                                         'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%xfang) CALL define_ovar(ncid_out, opid%xfang, &
                  'xfang', '-', 'Leaf angle parameter',patchout%xfang, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%wai) CALL define_ovar(ncid_out, opid%wai,     &
                          'wai', '-', 'Wood area index', patchout%wai, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%vegcf) CALL define_ovar(ncid_out, opid%vegcf, &
                                'vegcf', '-', 'vegcf', patchout%vegcf, 'real', &
                                                 xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%extkn) CALL define_ovar(ncid_out, opid%extkn, &
            'extkn', '-', 'Nitrogen extinction coef for vert. canopy profile', &
                         patchout%extkn, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%tminvj) CALL define_ovar(ncid_out,            &
                                                   opid%tminvj, 'tminvj', 'C', &
                            'Min temperature for the start of photosynthesis', &
                        patchout%tminvj, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%tmaxvj) CALL define_ovar(ncid_out,            &
             opid%tmaxvj, 'tmaxvj', 'C', 'Max temperature for photosynthesis', &
                        patchout%tmaxvj, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%vbeta) CALL define_ovar(ncid_out, opid%vbeta, &
                           'vbeta', '-', 'Stomatal sensitivity to soil water', &
                         patchout%vbeta, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%xalbnir) CALL define_ovar(ncid_out,           &
          opid%xalbnir, 'xalbnir', '-', 'Modifier for albedo in near ir band', &
                       patchout%xalbnir, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%meth) CALL define_ovar(ncid_out, opid%meth,   &
                     'meth', '-', 'Canopy turbulence parameterisation choice', &
                          patchout%meth, 'real', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%za) THEN
      CALL define_ovar(ncid_out, opid%za_uv, 'za_uv', 'm',                     &
                    'Reference height (lowest atm. model layer) for momentum', &
                            patchout%za, 'real', xID, yID, zID, landID, patchID)
      CALL define_ovar(ncid_out, opid%za_tq, 'za_tq', 'm',                     &
                     'Reference height (lowest atm. model layer) for scalars', &
                            patchout%za, 'real', xID, yID, zID, landID, patchID)
    ENDIF
    IF(output%params .OR. output%ratecp) CALL define_ovar(ncid_out,            &
                opid%ratecp, 'ratecp', '1/year', 'Plant carbon rate constant', &
                   patchout%ratecp, plantcarbID, 'plantcarbon', xID, yID, zID, &
                                                                landID, patchID)
    IF(output%params .OR. output%ratecs) CALL define_ovar(ncid_out,            &
                 opid%ratecs, 'ratecs', '1/year', 'Soil carbon rate constant', &
                     patchout%ratecs, soilcarbID, 'soilcarbon', xID, yID, zID, &
                                                                landID, patchID)
    IF(output%params .OR. output%zse) CALL define_ovar(ncid_out, opid%zse,     &
                                       'zse', 'm', 'Depth of each soil layer', &
                   patchout%zse, soilID, 'soil', xID, yID, zID, landID, patchID)
    IF(output%params .OR. output%froot) CALL define_ovar(ncid_out, opid%froot, &
                         'froot', '-', 'Fraction of roots in each soil layer', &
                 patchout%froot, soilID, 'soil', xID, yID, zID, landID, patchID)
                 
         
    !MD
    IF(output%params .OR. output%WatSat) CALL define_ovar(ncid_out, opid%WatSat, &
                         'WatSat', '-', 'Max water content in soil layer', &
                 patchout%WatSat, soilID, 'soil', xID, yID, zID, landID, patchID)
                 
    !IF(output%params .OR. output%GWWatSat) CALL define_ovar(ncid_out, opid%GWWatSat, &
    !                     'GWWatSat', '-', 'Max water content in aquifer', &
    !             patchout%GWWatSat, 'real', xID, yID, zID, landID, patchID)    

    IF(output%params .OR. output%GWWatSat) &
              CALL define_ovar(ncid_out, opid%GWWatSat, 'GWWatSat', 'mm3mm-3', &
                   'Aquifer Maximum porosity',                                 &
                   patchout%GWWatSat, 'real', xID, yID, zID, landID, patchID)

  
                 
    IF(output%params .OR. output%Watr) CALL define_ovar(ncid_out, opid%Watr, &
                         'Watr', '-', 'residual water content in soil layer', &
                 patchout%Watr, soilID, 'soil', xID, yID, zID, landID, patchID)
                 
    IF(output%params .OR. output%GWWatr) CALL define_ovar(ncid_out, opid%GWWatr, &
                         'GWWatr', '-', 'residual water content in aquifer', &
                 patchout%GWWatr, 'real', xID, yID, zID, landID, patchID)                  
                 
    IF(output%params .OR. output%SoilMatPotSat) CALL define_ovar(ncid_out, opid%SoilMatPotSat, &
                         'SoilMatPotSat', '-', 'soil matric potent at saturation content in soil layer', &
                 patchout%SoilMatPotSat, soilID, 'soil', xID, yID, zID, landID, patchID)    
                 
    IF(output%params .OR. output%GWSoilMatPotSat) CALL define_ovar(ncid_out, opid%GWSoilMatPotSat, &
                         'GWSoilMatPotSat', '-', 'soil matric potent at saturation content in aquifer', &
                 patchout%GWSoilMatPotSat, 'real', xID, yID, zID, landID, patchID)
                 
    IF(output%params .OR. output%HkSat) CALL define_ovar(ncid_out, opid%HkSat, &
                         'HkSat', '-', 'Max hydraulic conductivity in soil layer', &
                 patchout%HkSat, soilID, 'soil', xID, yID, zID, landID, patchID)
                 
    IF(output%params .OR. output%GWHkSat) CALL define_ovar(ncid_out, opid%GWHkSat, &
                         'GWHkSat', '-', 'Max hydraulic conductivityin aquifer', &
                 patchout%GWHkSat, 'real', xID, yID, zID, landID, patchID)  
                 
    IF(output%params .OR. output%FrcSand) CALL define_ovar(ncid_out, opid%FrcSand, &
                         'FrcSand', '-', 'sand fraction  in soil layer', &
                 patchout%FrcSand, soilID, 'soil', xID, yID, zID, landID, patchID)       
                 
    IF(output%params .OR. output%FrcClay) CALL define_ovar(ncid_out, opid%FrcClay, &
                         'FrcClay', '-', 'clay fraction  in soil layer', &
                 patchout%FrcClay, soilID, 'soil', xID, yID, zID, landID, patchID)    
                 
    IF(output%params .OR. output%ClappB) CALL define_ovar(ncid_out, opid%ClappB, &
                         'ClappB', '-', 'clapp and horn b param  in soil layer', &
                 patchout%ClappB, soilID, 'soil', xID, yID, zID, landID, patchID)                     

    IF(output%params) CALL define_ovar(ncid_out, opid%fldcap,&
                         'FieldCap', 'mm3/mm3', 'field capcaicty from hk',&
                 patchout%fldcap, soilID, 'soil', xID, yID, zID, landID,patchID)

    IF(output%params) CALL define_ovar(ncid_out, opid%wiltp,&
                         'WiltP', 'mm3/mm3', 'wilting point from hk',&
                 patchout%wiltp, soilID, 'soil', xID, yID, zID, landID,patchID)

    IF(output%params) CALL define_ovar(ncid_out, opid%forg,&
                         'FrcOrg', 'mm3/mm3', 'organic fraction',&
                 patchout%forg, soilID, 'soil', xID, yID, zID, landID,patchID)


    !MDeck
    write(*,*) 'write global attributes for th file'


    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Production",                     &
                      TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out,NF90_GLOBAL,"Source", &
         'CABLE LSM output file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "CABLE_input_file",               &
                      TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Determine output aggregation details:
    IF(output%averaging(1:4) == 'user') THEN
       ! User-specified aggregation interval for output:
       READ(output%averaging(5:7), *) output%interval
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
                         TRIM(output%averaging(5:7))//'-hourly output')
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
                                             'Error writing global detail to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Convert interval value from hours to time steps (for use in output
       ! write):
       output%interval = output%interval * 3600/INT(dels)
    ELSE IF(output%averaging(1:3) == 'all') THEN ! output all timesteps
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
                         TRIM(output%averaging)//' timesteps recorded')
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
                                             'Error writing global detail to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Set output interval to be one time step
       output%interval = 1
    ELSE IF(output%averaging(1:2) == 'mo') THEN ! monthly output
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
                         TRIM(output%averaging))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
                                             'Error writing global detail to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Output interval will be determined dynamically for monthly output
    ELSE IF(output%averaging(1:2) == 'da') THEN ! daily output
       ok = NF90_PUT_ATT(ncid_out, NF90_GLOBAL, "Output_averaging",            &
                         TRIM(output%averaging))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
                                             'Error writing global detail to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
       ! Set output interval to be # time steps in 24 hours:
       output%interval = 3600*24/INT(dels)
    ELSE
       CALL abort ('Unknown output averaging interval specified '//            &
            'in namelist file. (SUBROUTINE open_output_file)')
    END IF

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating output file '       &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Write latitude and longitude variables:
    ok = NF90_PUT_VAR(ncid_out, latID, REAL(lat_all, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
                                    (ok, 'Error writing latitude variable to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_VAR(ncid_out, lonID, REAL(lon_all, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
                                   (ok, 'Error writing longitude variable to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ! Write GrADS coordinate variables
    ok = NF90_PUT_VAR(ncid_out, xvID, REAL(lon_all(:, 1), 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
                          (ok, 'Error writing GrADS x coordinate variable to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')
    ok = NF90_PUT_VAR(ncid_out, yvID, REAL(lat_all(1, :), 4))
    IF(ok /= NF90_NOERR) CALL nc_abort                                         &
                          (ok, 'Error writing GrADS y coordinate variable to ' &
                        //TRIM(filename%out)// ' (SUBROUTINE open_output_file)')

    ! Write model parameters if requested:
    IF(output%params .OR. output%iveg) CALL write_ovar(ncid_out, opid%iveg,    &
               'iveg', REAL(veg%iveg, 4), ranges%iveg, patchout%iveg, 'integer')
    IF((output%params .OR. output%patchfrac)                                   &
       .AND. (patchout%patchfrac .OR. output%patch))                           &
       CALL write_ovar(ncid_out, opid%patchfrac, 'patchfrac'                   &
               ,REAL(patch(:)%frac, 4), (/0.0, 1.0/), patchout%patchfrac, 'real')
    IF(output%params .OR. output%isoil) CALL write_ovar(ncid_out, opid%isoil,  &
          'isoil', REAL(soil%isoilm, 4), ranges%isoil, patchout%isoil,'integer')
    IF(output%params .OR. output%bch) CALL write_ovar(ncid_out, opid%bch,      &
                     'bch', REAL(soil%bch, 4), ranges%bch, patchout%bch, 'real')
    IF(output%params .OR. output%clay) CALL write_ovar(ncid_out, opid%clay,    &
                 'clay', REAL(soil%clay, 4), ranges%clay, patchout%clay, 'real')
    IF(output%params .OR. output%sand) CALL write_ovar(ncid_out, opid%sand,    &
                 'sand', REAL(soil%sand, 4), ranges%sand, patchout%sand, 'real')
    IF(output%params .OR. output%silt) CALL write_ovar(ncid_out, opid%silt,    &
                 'silt', REAL(soil%silt, 4), ranges%silt, patchout%silt, 'real')
    IF(output%params .OR. output%css) CALL write_ovar(ncid_out, opid%css,      &
                     'css', REAL(soil%css, 4), ranges%css, patchout%css, 'real')
    IF(output%params .OR. output%rhosoil) CALL write_ovar(ncid_out,            &
                                 opid%rhosoil, 'rhosoil',REAL(soil%rhosoil,4), &
                                       ranges%rhosoil, patchout%rhosoil, 'real')
    IF(output%params .OR. output%hyds) CALL write_ovar(ncid_out, opid%hyds,    &
                 'hyds', REAL(soil%hyds, 4), ranges%hyds, patchout%hyds, 'real')
    IF(output%params .OR. output%sucs) CALL write_ovar(ncid_out, opid%sucs,    &
                 'sucs', REAL(soil%sucs, 4), ranges%sucs, patchout%sucs, 'real')
    IF(output%params .OR. output%rs20) CALL write_ovar(ncid_out, opid%rs20,    &
                  'rs20', REAL(veg%rs20, 4), ranges%rs20, patchout%rs20, 'real')
!         'rs20',REAL(soil%rs20,4),ranges%rs20,patchout%rs20,'real')
    IF(output%params .OR. output%ssat) CALL write_ovar(ncid_out, opid%ssat,    &
                 'ssat', REAL(soil%ssat, 4), ranges%ssat, patchout%ssat, 'real')
    IF(output%params .OR. output%sfc) CALL write_ovar(ncid_out, opid%sfc,      &
                     'sfc', REAL(soil%sfc, 4), ranges%sfc, patchout%sfc, 'real')
    IF(output%params .OR. output%swilt) CALL write_ovar(ncid_out, opid%swilt,  &
             'swilt', REAL(soil%swilt, 4), ranges%swilt, patchout%swilt, 'real')
    IF(output%params .OR. output%albsoil) CALL write_ovar(ncid_out,            &
                               opid%albsoil, 'albsoil', REAL(soil%albsoil, 4), &
                                  ranges%albsoil, patchout%albsoil, 'radiation')
    IF(output%params .OR. output%canst1) CALL write_ovar(ncid_out,             &
                                   opid%canst1, 'canst1', REAL(veg%canst1, 4), &
                                         ranges%canst1, patchout%canst1, 'real')
    IF(output%params .OR. output%dleaf) CALL write_ovar(ncid_out, opid%dleaf,  &
              'dleaf', REAL(veg%dleaf, 4), ranges%dleaf, patchout%dleaf, 'real')
    IF(output%params .OR. output%ejmax) CALL write_ovar(ncid_out, opid%ejmax,  &
              'ejmax', REAL(veg%ejmax, 4), ranges%ejmax, patchout%ejmax, 'real')
    IF(output%params .OR. output%vcmax) CALL write_ovar(ncid_out, opid%vcmax,  &
              'vcmax', REAL(veg%vcmax, 4), ranges%vcmax, patchout%vcmax, 'real')
    IF(output%params .OR. output%frac4) CALL write_ovar(ncid_out, opid%frac4,  &
              'frac4', REAL(veg%frac4, 4), ranges%frac4, patchout%frac4, 'real')
    IF(output%params .OR. output%hc) CALL write_ovar(ncid_out, opid%hc,        &
                          'hc', REAL(veg%hc, 4), ranges%hc, patchout%hc, 'real')
    IF(output%params .OR. output%rp20) CALL write_ovar(ncid_out, opid%rp20,    &
                   'rp20', REAL(veg%rp20, 4),ranges%rp20, patchout%rp20, 'real')

    ! Ticket #56
    IF(output%params .OR. output%g0c3) CALL write_ovar(ncid_out, opid%g0c3,    &
                   'g0c3', REAL(veg%g0c3, 4),ranges%g0c3, patchout%g0c3, 'real')
    IF(output%params .OR. output%g0c4) CALL write_ovar(ncid_out, opid%g0c4,    &
                   'g0c4', REAL(veg%g0c4, 4),ranges%g0c4, patchout%g0c4, 'real')
    IF(output%params .OR. output%g1c3) CALL write_ovar(ncid_out, opid%g1c3,    &
                   'g1c3', REAL(veg%g1c3, 4),ranges%g1c3, patchout%g1c3, 'real')
    IF(output%params .OR. output%g1c4) CALL write_ovar(ncid_out, opid%g1c4,    &
                   'g1c4', REAL(veg%g1c4, 4),ranges%g1c4, patchout%g1c4, 'real')
    ! End Ticket #56

    IF(output%params .OR. output%rpcoef) CALL write_ovar(ncid_out,             &
                                   opid%rpcoef, 'rpcoef', REAL(veg%rpcoef, 4), &
                                         ranges%rpcoef, patchout%rpcoef, 'real')
    IF(output%params .OR. output%shelrb) CALL write_ovar(ncid_out,             &
                                   opid%shelrb, 'shelrb', REAL(veg%shelrb, 4), &
                                         ranges%shelrb, patchout%shelrb, 'real')
    IF(output%params .OR. output%xfang) CALL write_ovar(ncid_out, opid%xfang,  &
              'xfang', REAL(veg%xfang, 4), ranges%xfang, patchout%xfang, 'real')
    IF(output%params .OR. output%wai) CALL write_ovar(ncid_out, opid%wai,      &
                      'wai', REAL(veg%wai, 4), ranges%wai, patchout%wai, 'real')
    IF(output%params .OR. output%vegcf) CALL write_ovar(ncid_out, opid%vegcf,  &
              'vegcf', REAL(veg%vegcf, 4), ranges%vegcf, patchout%vegcf, 'real')
    IF(output%params .OR. output%extkn) CALL write_ovar(ncid_out, opid%extkn,  &
              'extkn', REAL(veg%extkn, 4), ranges%extkn, patchout%extkn, 'real')
    IF(output%params .OR. output%tminvj) CALL write_ovar(ncid_out,             &
                                   opid%tminvj, 'tminvj', REAL(veg%tminvj, 4), &
                                         ranges%tminvj, patchout%tminvj, 'real')
    IF(output%params .OR. output%tmaxvj) CALL write_ovar(ncid_out,             &
                                   opid%tmaxvj, 'tmaxvj', REAL(veg%tmaxvj, 4), &
                                         ranges%tmaxvj, patchout%tmaxvj, 'real')
    IF(output%params .OR. output%vbeta) CALL write_ovar(ncid_out, opid%vbeta,  &
              'vbeta', REAL(veg%vbeta, 4), ranges%vbeta, patchout%vbeta, 'real')
    IF(output%params .OR. output%xalbnir) CALL write_ovar(ncid_out,            &
                                opid%xalbnir, 'xalbnir', REAL(veg%xalbnir, 4), &
                                       ranges%xalbnir, patchout%xalbnir, 'real')
    IF(output%params .OR. output%meth) CALL write_ovar(ncid_out, opid%meth,    &
               'meth', REAL(veg%meth, 4), ranges%meth, patchout%meth, 'integer')
    IF(output%params .OR. output%za) THEN
      CALL write_ovar(ncid_out, opid%za_uv,                                    &
                  'za_uv', REAL(rough%za_uv, 4), ranges%za, patchout%za, 'real')
      CALL write_ovar(ncid_out, opid%za_tq,                                    &
                  'za_tq', REAL(rough%za_tq, 4), ranges%za, patchout%za, 'real')
    ENDIF
    IF(output%params .OR. output%ratecp) CALL write_ovar(ncid_out,             &
         opid%ratecp, 'ratecp',SPREAD(REAL(bgc%ratecp,4),1,mp), ranges%ratecp, &
                       patchout%ratecp,'plantcarbon')! no spatial dim at present
    IF(output%params .OR. output%ratecs) CALL write_ovar(ncid_out,             &
    opid%ratecs, 'ratecs', SPREAD(REAL(bgc%ratecs, 4), 1, mp), ranges%ratecs,  &
                       patchout%ratecs, 'soilcarbon')! no spatial dim at present
    IF(output%params .OR. output%froot) CALL write_ovar (ncid_out, opid%froot, &
              'froot', REAL(veg%froot, 4), ranges%froot, patchout%froot, 'soil')
    IF(output%params .OR. output%zse) CALL write_ovar(ncid_out, opid%zse,      &
                           'zse', SPREAD(REAL(soil%zse, 4), 1, mp),ranges%zse, &
                                patchout%zse, 'soil')! no spatial dim at present
                                
    !MD
    IF(output%params .OR. output%WatSat) CALL write_ovar (ncid_out, opid%WatSat, &
              'WatSat', REAL(soil%watsat, 4), ranges%WatSat, patchout%WatSat, 'soil')
    IF(output%params .OR. output%Watr) CALL write_ovar (ncid_out, opid%Watr, &
              'Watr', REAL(soil%watr, 4), ranges%Watr, patchout%Watr, 'soil')

    IF(output%params .OR. output%GWWatSat) CALL write_ovar (ncid_out, opid%GWWatSat, &
              'GWWatSat', REAL(soil%GWwatsat, 4), ranges%GWWatSat,                   &
                                                            patchout%GWWatSat, 'real')

    IF(output%params .OR. output%GWWatr) CALL write_ovar (ncid_out, opid%GWWatr, &
              'GWWatr', REAL(soil%GWwatr, 4), ranges%GWWatr, patchout%GWWatr, 'real')     
    IF(output%params .OR. output%SoilMatPotSat) CALL write_ovar (ncid_out, opid%SoilMatPotSat, &
              'SoilMatPotSat', REAL(soil%smpsat, 4), ranges%SoilMatPotSat, patchout%SoilMatPotSat, 'soil')   

    IF(output%params .OR. output%GWSoilMatPotSat) &
               CALL write_ovar (ncid_out, opid%GWSoilMatPotSat, &
              'GWSoilMatPotSat', REAL(soil%GWsmpsat, 4), ranges%GWSoilMatPotSat, &
              patchout%GWSoilMatPotSat, 'real')    
 
    IF(output%params .OR. output%HkSat) CALL write_ovar (ncid_out, opid%HkSat, &
              'HkSat', REAL(soil%hksat, 4), ranges%HkSat, patchout%HkSat, 'soil')       
    IF(output%params .OR. output%GWHkSat) CALL write_ovar (ncid_out, opid%GWHkSat, &
              'GWHkSat', REAL(soil%GWhksat, 4), ranges%GWHkSat, patchout%GWHkSat, 'real')                  
    IF(output%params .OR. output%FrcSand) CALL write_ovar (ncid_out, opid%FrcSand, &
              'FrcSand', REAL(soil%FSand, 4), ranges%FrcSand, patchout%FrcSand, 'soil')          
    IF(output%params .OR. output%FrcClay) CALL write_ovar (ncid_out, opid%FrcClay, &
              'FrcClay', REAL(soil%FClay, 4), ranges%FrcClay, patchout%FrcClay, 'soil')          
    IF(output%params .OR. output%ClappB) CALL write_ovar (ncid_out, opid%ClappB, &
              'ClappB', REAL(soil%clappB, 4), ranges%ClappB, patchout%ClappB, 'soil')    

    IF(output%params .and. cable_user%GW_MODEL) CALL write_ovar (ncid_out,opid%fldcap, &
              'FieldCap', REAL(soil%fldcap, 4), (/0.,1./), patchout%fldcap,'soil')
    IF(output%params .and. cable_user%GW_MODEL) CALL write_ovar (ncid_out,opid%wiltp, &
              'WiltPoint', REAL(soil%wiltp, 4), (/0.,1./), patchout%wiltp,'soil')
    IF(output%params .and. cable_user%GW_MODEL) CALL write_ovar (ncid_out,opid%forg, &
              'Frcorg', REAL(soil%forg, 4), (/0.,1./), patchout%forg,'soil')

              
  END SUBROUTINE open_output_file
  !=============================================================================
  SUBROUTINE write_output(dels, ktau, met, canopy, ssnow,                       &
                          rad, bal, air, soil, veg, rough,SBOLTZ, EMLEAF, EMSOIL)
    ! Writes model output variables and, if requested, calls
    ! energy and mass balance routines. This subroutine is called 
    ! each timestep, but may only write to the output file periodically,
    ! depending on whether the user has specified that output should be 
    ! aggregated, e.g. to monthly or 6-hourly averages.
    REAL, INTENT(IN)              :: dels ! time step size
    INTEGER, INTENT(IN)           :: ktau ! timestep number in loop which include spinup 
    REAL, INTENT(IN) :: SBOLTZ, EMLEAF, EMSOIL
    TYPE(met_type), INTENT(IN)         :: met  ! met data
    TYPE(canopy_type), INTENT(IN)      :: canopy ! canopy variable data
    TYPE(soil_snow_type), INTENT(IN)   :: ssnow ! soil data
    TYPE(soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE(radiation_type), INTENT(IN)  :: rad   ! radiation data
    TYPE(air_type), INTENT(IN)        :: air
    TYPE(veg_parameter_type), INTENT(IN) :: veg ! vegetation parameters
    TYPE (roughness_type), INTENT(IN)      :: rough
    TYPE(balances_type), INTENT(INOUT) :: bal

    REAL(r_2), DIMENSION(1) :: timetemp ! temporary variable for storing time
                                        ! value
    LOGICAL :: writenow ! write to output file during this time step?
    INTEGER, SAVE :: out_timestep ! counter for output time steps
    INTEGER, SAVE :: out_month ! counter for output month
    INTEGER, DIMENSION(mp) :: realyear ! fix problem for yr b4 leap yr
    INTEGER :: backtrack  ! modify timetemp for averaged output


    !write(*,*) ' in write_output '    !MDeck

    ! IF asked to check mass/water balance:
    IF(check%mass_bal) CALL mass_balance(dels, ktau, ssnow, soil, canopy,            &
                                         met,air,bal)

    ! IF asked to check energy balance:
    IF(check%energy_bal) CALL energy_balance(dels,met,rad,                     &
                                             canopy,bal,ssnow,                 &
                                             SBOLTZ, EMLEAF, EMSOIL ) 

    ! Initialise output time step counter and month counter:
    IF(ktau == 1) THEN
       out_timestep = 0
       out_month = 0
    END IF

    ! Decide on output averaging regime:
    IF(output%averaging(1:3) == 'all') THEN ! write every time step to file
       ! Set flag to write data for current time step:
       writenow = .TRUE.
       ! Set output time step to be current model time step:
       out_timestep = ktau
       backtrack = 0
    ELSE IF(output%averaging(1:4) == 'user' .OR. output%averaging(1:2)=='da')  &
            THEN
       ! user defined output interval or daily output
       IF(MOD(ktau, output%interval) == 0) THEN ! i.e.ktau divisible by
                                                   ! interval
          ! write to output file this time step 
          writenow = .TRUE.
          ! increment output time step counter:
          out_timestep = out_timestep + 1
          backtrack = output%interval / 2
       ELSE
          writenow = .FALSE.
       END IF
    ELSE IF(output%averaging(1:2) == 'mo') THEN ! write monthly averages to file
       realyear = met%year
       IF(ktau >= 365*24*3600/INT(dels)) THEN
         WHERE(met%doy == 1) realyear = realyear - 1   ! last timestep of year
       END IF
       
       ! Are we using leap year calendar?
       IF(leaps) THEN
          ! If currently a leap year:
          IF(((ANY(MOD(realyear,4)==0).AND.ANY(MOD(realyear,100)/=0)).OR. &
               (ANY(MOD(realyear,4)==0).AND.ANY(MOD(realyear,400)==0)))) THEN
             
             IF(ANY((lastdayl * 24 * 3600 / INT(dels)) == ktau)) THEN
                ! increment output month counter
                out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
                ! write to output file this time step 
                writenow = .TRUE.
                ! increment output time step counter:
                out_timestep = out_timestep + 1
                ! set numbr of time steps in output period
                output%interval = daysml(out_month) * 24 * 3600 / INT(dels)
             ELSE
                writenow = .FALSE.
             END IF
          ELSE ! not currently a leap year
             ! last time step of month
             IF(ANY((lastday * 24 * 3600 / INT(dels)) == ktau)) THEN
                ! increment output month counter
                out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
                ! write to output file this time step 
                writenow = .TRUE.
                ! increment output time step counter:
                out_timestep = out_timestep + 1
                ! set numbr of time steps in output period
                output%interval = daysm(out_month) * 24 * 3600 / INT(dels)
             ELSE
                writenow = .FALSE.
             END IF
          END IF
       ELSE ! not using leap year timing in this run
          IF(ANY((lastday*24*3600/INT(dels))==ktau)) THEN ! last time step of
                                                             ! month
             ! increment output month counter
             out_month = MOD(out_month, 12) + 1 ! can only be 1 - 12
             ! write to output file this time step 
             writenow = .TRUE.
             ! increment output time step counter:
             out_timestep = out_timestep + 1
             ! set numbr of time steps in output period
             output%interval = daysm(out_month) * 24 * 3600 / INT(dels)
          ELSE
             writenow = .FALSE.
          END IF
       END IF ! using leap year timing or not
       backtrack = output%interval / 2
    ELSE ! type of output aggregation
       CALL abort('Unknown output averaging request in namelist file.'//       &
                  '(SUBROUTINE write_output)')
    END IF

    ! Note that size of averaging interval, output%interval, is set when opening 
    ! output file unless output is monthly (in which case it's set above)

    ! If this time step is an output time step:
    IF(writenow) THEN
       ! Write to temporary time variable:
       timetemp(1) = DBLE(REAL(ktau-backtrack)*dels)
       ! Write time variable for this output time step:
       ok = NF90_PUT_VAR(ncid_out, ovid%tvar, timetemp,                        &
                                        start = (/out_timestep/), count = (/1/))
       IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                  &
                                             'Error writing time variable to ' &
                             //TRIM(filename%out)// '(SUBROUTINE write_output)')
    END IF

    ! Arguments to write_ovar: current time step; output file netcdf file ID;
    ! netcdf variable ID; variable name; variable data; variable ranges; 
    ! non-land fill value; include patch info for this var; any specific 
    ! formatting info; met variables for reporting in case of abort.

    ! SWdown:  downward short-wave radiation [W/m^2]
    IF(output%met .OR. output%SWdown) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SWdown = out%SWdown + REAL(met%fsd(:, 1) + met%fsd(:, 2), 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SWdown = out%SWdown/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SWdown, 'SWdown',       &
                     out%SWdown, ranges%SWdown, patchout%SWdown, 'default', met)
          ! Reset temporary output variable:
          out%SWdown = 0.0
       END IF
    END IF
    ! LWdown: downward long-wave radiation [W/m^2]
    IF(output%met .OR. output%LWdown) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LWdown = out%LWdown + REAL(met%fld, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LWdown = out%LWdown/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LWdown, 'LWdown',       &
                     out%LWdown, ranges%LWdown, patchout%LWdown, 'default', met)
          ! Reset temporary output variable:
          out%LWdown = 0.0
       END IF
    END IF
    ! Tair: surface air temperature [K]
    IF(output%met .OR. output%Tair) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Tair = out%Tair + REAL(met%tk, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Tair = out%Tair/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Tair, 'Tair', out%Tair, &
                                        ranges%Tair, patchout%Tair, 'ALMA', met)
          ! Reset temporary output variable:
          out%Tair = 0.0
       END IF
    END IF
    ! Rainf: rainfall [kg/m^2/s]
    IF(output%met .OR. output%Rainf) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Rainf = out%Rainf + REAL(met%precip / dels, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Rainf = out%Rainf/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Rainf, 'Rainf',         &
                        out%Rainf, ranges%Rainf, patchout%Rainf, 'default', met)
          ! Reset temporary output variable:
          out%Rainf = 0.0
       END IF
    END IF
    ! Snowf: snowfall [kg/m^2/s]
    IF(output%met .OR. output%Snowf) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Snowf = out%Snowf + REAL(met%precip_sn / dels, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Snowf = out%Snowf/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Snowf, 'Snowf',         &
                        out%Snowf, ranges%Snowf, patchout%Snowf, 'default', met)
          ! Reset temporary output variable:
          out%Snowf = 0.0
       END IF
    END IF
    ! PSurf: surface pressure [Pa]
    IF(output%met .OR. output%PSurf) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%PSurf = out%PSurf + REAL(met%pmb, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%PSurf = out%PSurf / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%PSurf, 'PSurf',         &
                        out%PSurf, ranges%PSurf, patchout%PSurf, 'default', met)
          ! Reset temporary output variable:
          out%PSurf = 0.0
       END IF
    END IF
    ! Qair: specific humidity [kg/kg]
    IF(output%met .OR. output%Qair) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qair = out%Qair + REAL(met%qv, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qair = out%Qair / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qair, 'Qair', out%Qair, &
                          ranges%Qair, patchout%Qair, 'ALMA', met)
          ! Reset temporary output variable:
          out%Qair = 0.0
       END IF
    END IF
    ! Wind: windspeed [m/s]
    IF(output%met .OR. output%Wind) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Wind = out%Wind + REAL(met%ua, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Wind = out%Wind/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Wind, 'Wind', out%Wind, &
                          ranges%Wind, patchout%Wind, 'ALMA', met)
          ! Reset temporary output variable:
          out%Wind = 0.0
       END IF
    END IF
    ! CO2air: CO2 concentration [ppmv]
    IF(output%met .OR. output%CO2air) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%CO2air = out%CO2air + REAL(met%ca * 1000000.0, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%CO2air = out%CO2air / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%CO2air, 'CO2air',       &
                        OUT%CO2air, ranges%CO2air, patchout%CO2air, 'ALMA', met)
          ! Reset temporary output variable:
          out%CO2air = 0.0
       END IF
    END IF
    !-----------------------WRITE FLUX DATA-------------------------------------
    ! Qle: latent heat flux [W/m^2]
    IF(output%flux .OR. output%Qle) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qle = out%Qle + REAL(canopy%fe, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qle = out%Qle / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qle, 'Qle', out%Qle,    &
                          ranges%Qle, patchout%Qle, 'default', met)
          ! Reset temporary output variable:
          out%Qle = 0.0
       END IF
    END IF
    ! Qh: sensible heat flux [W/m^2]
    IF(output%flux .OR. output%Qh) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qh = out%Qh + REAL(canopy%fh, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qh = out%Qh / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qh, 'Qh', out%Qh,       &
                          ranges%Qh, patchout%Qh, 'default', met)
          ! Reset temporary output variable:
          out%Qh = 0.0
       END IF
    END IF
    ! cancd: canopy conductance [m/s^1]
    IF(output%flux .OR. output%cancd) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%cancd = out%cancd + REAL(canopy%gswx_T, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%cancd = out%cancd / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%cancd, 'cancd',out%cancd,       &
                          ranges%cancd, patchout%cancd, 'default', met)
          ! Reset temporary output variable:
          out%cancd = 0.0
       END IF
    END IF
    ! gswx_1: sunlit conductance [NA]
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswx_1)) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%gswx_1 = out%gswx_1 + REAL(canopy%gswx_1, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%gswx_1 = out%gswx_1 / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%gswx_1,'gswx_1',out%gswx_1,       &
                          ranges%gswx_1, patchout%gswx_1, 'default', met)
          ! Reset temporary output variable:
          out%gswx_1 = 0.0
       END IF
    END IF
    ! gswx_2: sunlit conductance [NA]
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswx_2)) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%gswx_2 = out%gswx_2 + REAL(canopy%gswx_2, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%gswx_2 = out%gswx_2 / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out,ovid%gswx_2,'gswx_2',out%gswx_2,       &
                          ranges%gswx_2, patchout%gswx_2, 'default', met)
          ! Reset temporary output variable:
          out%gswx_2 = 0.0
       END IF
    END IF
    ! gswmin_1: sunlit conductance [NA]
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswmin_1)) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%gswmin_1 = out%gswmin_1 + REAL(canopy%gswmin_1, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%gswmin_1 = out%gswmin_1 / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out,ovid%gswmin_1,'gswmin_1',out%gswmin_1,       &
                          ranges%gswmin_1, patchout%gswmin_1, 'default', met)
          ! Reset temporary output variable:
          out%gswmin_1 = 0.0
       END IF
    END IF
    ! gswmin_2: sunlit conductance [NA]
    if ((mp .eq. 1) .and. (output%flux .OR. output%gswmin_2)) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%gswmin_2 = out%gswmin_2 + REAL(canopy%gswmin_2, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%gswmin_2 = out%gswmin_2 / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep,ncid_out,ovid%gswmin_2,'gswmin_2',out%gswmin_2,       &
                          ranges%gswmin_2, patchout%gswmin_2, 'default', met)
          ! Reset temporary output variable:
          out%gswmin_2 = 0.0
       END IF
    END IF

    ! Qg: ground heat flux [W/m^2]
    IF(output%flux .OR. output%Qg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qg = out%Qg + REAL(canopy%ga, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qg = out%Qg / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qg, 'Qg', out%Qg,       &
                          ranges%Qg, patchout%Qg, 'default', met)
          ! Reset temporary output variable:
          out%Qg = 0.0
       END IF
    END IF
    ! Qs: surface runoff [kg/m^2/s]
    IF(output%flux .OR. output%Qs) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qs = out%Qs + ssnow%rnof1 / dels
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qs = out%Qs / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qs, 'Qs', out%Qs,       &
                          ranges%Qs, patchout%Qs, 'default', met)
          ! Reset temporary output variable:
          out%Qs = 0.0
       END IF
    END IF
    ! Qsb: subsurface runoff [kg/m^2/s]
    IF(output%flux .OR. output%Qsb) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qsb = out%Qsb + REAL(ssnow%rnof2 / dels, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qsb = out%Qsb / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qsb, 'Qsb', out%Qsb,    &
                          ranges%Qsb, patchout%Qsb, 'default', met)
          ! Reset temporary output variable:
          out%Qsb = 0.0
       END IF
    END IF
    ! Evap: total evapotranspiration [kg/m^2/s]
    IF(output%flux .OR. output%Evap) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Evap = out%Evap + REAL(canopy%fe / air%rlam, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Evap = out%Evap / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Evap, 'Evap', out%Evap, &
                          ranges%Evap, patchout%Evap, 'default', met)
          ! Reset temporary output variable:
          out%Evap = 0.0
       END IF
    END IF
    ! ECanop: interception evaporation [kg/m^2/s]
    IF(output%flux .OR. output%ECanop) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%ECanop = out%ECanop + REAL(canopy%fevw / air%rlam, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%ECanop = out%ECanop / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%ECanop, 'ECanop',       &
                     out%ECanop, ranges%ECanop, patchout%ECanop, 'default', met)
          ! Reset temporary output variable:
          out%ECanop = 0.0
       END IF
    END IF
    ! TVeg: vegetation transpiration [kg/m^2/s]
    IF(output%flux .OR. output%TVeg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%TVeg = out%TVeg + REAL(canopy%fevc / air%rlam, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%TVeg = out%TVeg / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%TVeg, 'TVeg', out%TVeg, &
                          ranges%TVeg, patchout%TVeg, 'default', met)
          ! Reset temporary output variable:
          out%TVeg = 0.0
       END IF
    END IF
    ! ESoil: bare soil evaporation [kg/m^2/s]
    IF(output%flux .OR. output%ESoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%ESoil = out%ESoil + REAL(canopy%fes / air%rlam, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%ESoil = out%ESoil / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%ESoil, 'ESoil',         &
                        out%ESoil, ranges%ESoil, patchout%ESoil, 'default', met)
          ! Reset temporary output variable:
          out%ESoil = 0.0
       END IF
    END IF
    ! HVeg: sensible heat from vegetation [W/m^2]
    IF(output%flux .OR. output%HVeg) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%HVeg = out%HVeg + REAL(canopy%fhv, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%HVeg = out%HVeg/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%HVeg, 'HVeg', out%HVeg, &
                          ranges%HVeg, patchout%HVeg, 'default', met)
          ! Reset temporary output variable:
          out%HVeg = 0.0
       END IF
    END IF
    ! HSoil: sensible heat from soil [W/m^2]
    IF(output%flux .OR. output%HSoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%HSoil = out%HSoil + REAL(canopy%fhs, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%HSoil = out%HSoil / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%HSoil, 'HSoil',         &
                        out%HSoil, ranges%HSoil, patchout%HSoil, 'default', met)
          ! Reset temporary output variable:
          out%HSoil = 0.0
       END IF
    END IF
    ! NEE: net ecosystem exchange [umol/m^2/s]
    IF(output%flux .OR. output%carbon .OR. output%NEE) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%NEE = out%NEE + REAL(canopy%fnee / 1.201E-5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%NEE = out%NEE / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%NEE, 'NEE', out%NEE,    &
                          ranges%NEE, patchout%NEE, 'default', met)
          ! Reset temporary output variable:
          out%NEE = 0.0
       END IF
    END IF
    !-----------------------WRITE SOIL STATE DATA-------------------------------
    ! SoilMoist: av.layer soil moisture [kg/m^2]
    IF(output%soil .OR. output%SoilMoist) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SoilMoist = out%SoilMoist + REAL(ssnow%wb, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilMoist = out%SoilMoist / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilMoist, 'SoilMoist', &
               out%SoilMoist, ranges%SoilMoist, patchout%SoilMoist, 'soil', met)
          ! Reset temporary output variable:
          out%SoilMoist = 0.0
       END IF
    END IF
    ! SoilTemp: av.layer soil temperature [K]
    IF(output%soil .OR. output%SoilTemp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SoilTemp = out%SoilTemp + REAL(ssnow%tgg, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilTemp = out%SoilTemp/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilTemp, 'SoilTemp',   &
                  out%SoilTemp, ranges%SoilTemp, patchout%SoilTemp, 'soil', met)
          ! Reset temporary output variable:
          out%SoilTemp = 0.0
       END IF
    END IF
    ! BaresoilT: surface bare soil temp [K]
    IF(output%soil .OR. output%BaresoilT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%BaresoilT = out%BaresoilT + REAL(ssnow%tgg(:, 1), 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%BaresoilT = out%BaresoilT / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%BaresoilT, 'BaresoilT', &
            out%BaresoilT, ranges%BaresoilT, patchout%BaresoilT, 'default', met)
          ! Reset temporary output variable:
          out%BaresoilT = 0.0
       END IF
    END IF
    !----------------------WRITE SNOW STATE DATA--------------------------------
    ! SWE: snow water equivalent [kg/m^2]
    IF(output%snow .OR. output%SWE) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SWE = out%SWE + REAL(ssnow%snowd, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SWE = out%SWE / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SWE, 'SWE', out%SWE,    &
                          ranges%SWE, patchout%SWE, 'default', met)
          ! Reset temporary output variable:
          out%SWE = 0.0
       END IF
    END IF
    ! SnowT: snow surface temp [K]
    IF(output%snow .OR. output%SnowT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SnowT = out%SnowT + REAL(ssnow%tggsn(:, 1), 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SnowT = out%SnowT / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SnowT, 'SnowT',         &
                        out%SnowT, ranges%SnowT, patchout%SnowT, 'default', met)
          ! Reset temporary output variable:
          out%SnowT = 0.0
       END IF
    END IF
    ! SnowDepth: actual depth of snow in [m]
    IF(output%snow .OR. output%SnowDepth) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SnowDepth = out%SnowDepth + REAL(SUM(ssnow%sdepth, 2), 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SnowDepth = out%SnowDepth/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SnowDepth, 'SnowDepth', &
            out%SnowDepth, ranges%SnowDepth, patchout%SnowDepth, 'default', met)
          ! Reset temporary output variable:
          out%SnowDepth = 0.0
       END IF
    END IF
    !-------------------------WRITE RADIATION DATA------------------------------
    ! SWnet: net shortwave [W/m^2]
    IF(output%radiation .OR. output%SWnet) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SWnet = out%SWnet + REAL(SUM(rad%qcan(:, :, 1), 2) +                &
                                      SUM(rad%qcan(:, :, 2), 2) + rad%qssabs, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SWnet = out%SWnet / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SWnet, 'SWnet',         &
                        out%SWnet, ranges%SWnet, patchout%SWnet, 'default', met)
          ! Reset temporary output variable:
          out%SWnet = 0.0
       END IF
    END IF
    ! LWnet: net longwave [W/m^2]
    IF(output%radiation .OR. output%LWnet) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LWnet = out%LWnet +                                                 &
          REAL(met%fld - sboltz * emleaf * canopy%tv ** 4 * (1 - rad%transd) - &
               rad%flws * rad%transd, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LWnet = out%LWnet / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LWnet, 'LWnet',         &
                        out%LWnet, ranges%LWnet, patchout%LWnet, 'default', met)
          ! Reset temporary output variable:
          out%LWnet = 0.0
       END IF
    END IF
    ! Rnet: net absorbed radiation [W/m^2]
    IF(output%radiation .OR. output%Rnet) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Rnet = out%Rnet + REAL(met%fld - sboltz * emleaf * canopy%tv ** 4 * &
                                  (1 - rad%transd) -rad%flws * rad%transd +    &
                                  SUM(rad%qcan(:, :, 1), 2) +                  &
                                  SUM(rad%qcan(:, :, 2), 2) + rad%qssabs, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Rnet = out%Rnet / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Rnet, 'Rnet', out%Rnet, &
                          ranges%Rnet, patchout%Rnet, 'default', met)
          ! Reset temporary output variable:
          out%Rnet = 0.0
       END IF
    END IF
    ! Albedo:
    IF(output%radiation .OR. output%Albedo) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Albedo = out%Albedo + REAL((rad%albedo(:, 1) + rad%albedo(:, 2))    &
                                       * 0.5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Albedo = out%Albedo / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Albedo, 'Albedo',       &
                     out%Albedo, ranges%Albedo, patchout%Albedo, 'default', met)
          ! Reset temporary output variable:
          out%Albedo = 0.0
       END IF
    END IF
    ! VISalbedo:
    IF(output%radiation .OR. output%VISalbedo) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%VISalbedo = out%VISalbedo + REAL(rad%albedo(:, 1), 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%VISalbedo = out%VISalbedo / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%VISalbedo, 'VISalbedo',       &
                     out%VISalbedo, ranges%Albedo, patchout%VISalbedo, 'default', met)
          ! Reset temporary output variable:
          out%VISalbedo = 0.0
       END IF
    END IF
    ! NIRalbedo:
    IF(output%radiation .OR. output%NIRalbedo) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%NIRalbedo = out%NIRalbedo + REAL(rad%albedo(:, 2), 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%NIRalbedo = out%NIRalbedo / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%NIRalbedo, 'NIRalbedo',       &
                     out%NIRalbedo, ranges%Albedo, patchout%NIRalbedo, 'default', met)
          ! Reset temporary output variable:
          out%NIRalbedo = 0.0
       END IF
    END IF
    ! RadT: Radiative surface temperature [K]
    IF(output%radiation .OR. output%RadT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%RadT = out%RadT + REAL((((1.0 - rad%transd) * emleaf * sboltz *     &
       canopy%tv ** 4 + rad%transd * emsoil * sboltz * ((1 - ssnow%isflag) *   &
       ssnow%tgg(:, 1) + ssnow%isflag * ssnow%tggsn(:, 1)) ** 4) / sboltz)     &
       ** 0.25, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%RadT = out%RadT/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%RadT, 'RadT', out%RadT, &
                          ranges%RadT, patchout%RadT, 'default', met)
          ! Reset temporary output variable:
          out%RadT = 0.0
       END IF
    END IF
    !------------------------WRITE VEGETATION DATA------------------------------
    ! VegT: vegetation temperature [K]
    IF(output%veg .OR. output%VegT) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%VegT = out%VegT + REAL(canopy%tv, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%VegT = out%VegT / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%VegT, 'VegT', out%VegT, &
                          ranges%VegT, patchout%VegT, 'default', met)
          ! Reset temporary output variable:
          out%VegT = 0.0
       END IF
    END IF
    ! CanopInt: total canopy water storage [kg/m^2]
    IF(output%veg .OR. output%CanopInt) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%CanopInt = out%CanopInt + REAL(canopy%cansto, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%CanopInt = out%CanopInt / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%CanopInt, 'CanopInt',   &
               out%CanopInt, ranges%CanopInt, patchout%CanopInt, 'default', met)
          ! Reset temporary output variable:
          out%CanopInt = 0.0
       END IF
    END IF
    ! LAI:
    IF(output%veg .OR. output%LAI) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LAI = out%LAI + REAL(veg%vlai, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LAI = out%LAI/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LAI, 'LAI', out%LAI,    &
                          ranges%LAI, patchout%LAI, 'default', met)
          ! Reset temporary output variable:
          out%LAI = 0.0
       END IF
    END IF
    !------------------------WRITE BALANCES DATA--------------------------------
    ! Ebal: cumulative energy balance [W/m^2]
    IF(output%balances .OR. output%Ebal) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Ebal = out%Ebal + REAL(bal%ebal_tot, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Ebal = out%Ebal / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Ebal, 'Ebal', out%Ebal, &
                          ranges%Ebal, patchout%Ebal, 'default', met)
          ! Reset temporary output variable:
          out%Ebal = 0.0
       END IF
    END IF
    ! Wbal: cumulative water balance  [kg/m^2/s]
    IF(output%balances .OR. output%Wbal) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Wbal = out%Wbal + REAL(bal%wbal_tot, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Wbal = out%Wbal / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Wbal, 'Wbal',           &
                          out%Wbal, ranges%Wbal, patchout%Wbal, 'default', met)
          ! Reset temporary output variable:
          out%Wbal = 0.0
       END IF
    END IF
    !------------------------WRITE CARBON DATA----------------------------------
    ! GPP: gross primary production C by veg [umol/m^2/s]
    !      added frday in the calculation of GPP (BP may08)
    IF(output%carbon .OR. output%GPP) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%GPP = out%GPP + REAL((-1.0 * canopy%fpn + canopy%frday)             &
                                / 1.201E-5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%GPP = out%GPP/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%GPP, 'GPP', out%GPP,    &
               ranges%GPP, patchout%GPP, 'default', met)
          ! Reset temporary output variable:
          out%GPP = 0.0
       END IF
    END IF
    ! NPP: net primary production of C by veg [umol/m^2/s]
    IF(output%carbon .OR. output%NPP) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%NPP = out%NPP + REAL((-1.0 * canopy%fpn - canopy%frp) / 1.201E-5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%NPP = out%NPP / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%NPP, 'NPP', out%NPP,    &
                          ranges%NPP, patchout%NPP, 'default', met)
          ! Reset temporary output variable:
          out%NPP = 0.0
       END IF
    END IF
    ! AutoResp: autotrophic respiration [umol/m^2/s]
    IF(output%carbon .OR. output%AutoResp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%AutoResp = out%AutoResp + REAL((canopy%frp + canopy%frday)          &
                                           / 1.201E-5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%AutoResp = out%AutoResp/REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%AutoResp, 'AutoResp',   &
               out%AutoResp, ranges%AutoResp, patchout%AutoResp, 'default', met)
          ! Reset temporary output variable:
          out%AutoResp = 0.0
       END IF
    END IF
    ! LeafResp: Leaf respiration [umol/m^2/s]
    IF(output%carbon .OR. output%LeafResp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%LeafResp = out%LeafResp + REAL(canopy%frday / 1.201E-5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%LeafResp = out%LeafResp / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%LeafResp, 'LeafResp',   &
               out%LeafResp, ranges%LeafResp, patchout%LeafResp, 'default', met)
          ! Reset temporary output variable:
          out%LeafResp = 0.0
       END IF
    END IF
    ! HeteroResp: heterotrophic respiration [umol/m^2/s]
    IF(output%carbon .OR. output%HeteroResp) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%HeteroResp = out%HeteroResp + REAL(canopy%frs / 1.201E-5, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%HeteroResp = out%HeteroResp / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%HeteroResp,             &
                          'HeteroResp', out%HeteroResp, ranges%HeteroResp,     &
                          patchout%HeteroResp, 'default', met)
          ! Reset temporary output variable:
          out%HeteroResp = 0.0
       END IF
    END IF
    
    
    
    !MD Write the hydrology output data from the groundwater module calculations
    !water table depth
    IF((output%soil .OR. output%WatTable) .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'wtd'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%WatTable = out%WatTable + REAL(ssnow%wtd/1000.0, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%WatTable = out%WatTable / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%WatTable, 'WatTable', &
               out%WatTable, ranges%WatTable, patchout%WatTable, 'default', met)
          ! Reset temporary output variable:
          out%WatTable = 0.0
       END IF
    END IF    
    !aquifer water content
    IF((output%soil .OR. output%GWMoist)  .and. cable_user%GW_MODEL) THEN
       out%GWMoist = out%GWMoist + REAL(ssnow%GWwb, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%GWMoist = out%GWMoist / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%GWMoist, 'GWMoist', &
               out%GWMoist, ranges%GWMoist, patchout%GWMoist, 'default', met)
          ! Reset temporary output variable:
          out%GWMoist = 0.0
       END IF
    END IF      

    !aquifer equilibrium water content
    IF((output%soil .or. output%EqGWMoist) .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'EQ-GWmoist'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%EqGWMoist = out%EqGWMoist + REAL(ssnow%GWwbeq, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%EqGWMoist = out%EqGWMoist / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%EqGWMoist, 'EqGWMoist', &
               out%EqGWMoist, ranges%EqGWMoist, patchout%EqGWMoist, 'default', met)
          ! Reset temporary output variable:
          out%EqGWMoist = 0.0
       END IF
    END IF        
    !aquifer soil matric potential
    IF((output%soil .or. output%GWSoilMatPot) .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'GW smp'  !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%GWSoilMatPot = out%GWSoilMatPot + REAL(ssnow%GWsmp/1000.0, 4)
       !write(*,*) minval(out%GWSoilMatPot),'  ',minval(ssnow%GWsmp)  !MDeck
       !write(*,*) maxval(out%GWSoilMatPot),'  ',maxval(ssnow%GWsmp)  !MDeck
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%GWSoilMatPot = out%GWSoilMatPot / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%GWSoilMatPot, 'GWSoilMatPot', &
               out%GWSoilMatPot, ranges%GWSoilMatPot, patchout%GWSoilMatPot, 'default', met)
          ! Reset temporary output variable:
          out%GWSoilMatPot = 0.0
       END IF
    END IF         
    !equilibrium aquifer soil matric potential
    IF((output%soil .or. output%EqGWSoilMatPot) .and. cable_user%GW_MODEL) THEN
!       write(*,*) 'EQ-GW smp'   !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%EqGWSoilMatPot = out%EqGWSoilMatPot + REAL(ssnow%GWzq/1000.0, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%EqGWSoilMatPot = out%EqGWSoilMatPot / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%EqGWSoilMatPot, 'EqGWSoilMatPot', &
               out%EqGWSoilMatPot, ranges%EqGWSoilMatPot, patchout%EqGWSoilMatPot, 'default', met)
          ! Reset temporary output variable:
          out%EqGWSoilMatPot = 0.0
       END IF
    END IF     
    ! soil matric potential
    IF((output%soil .OR. output%SoilMatPot) .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'smp'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%SoilMatPot = out%SoilMatPot + REAL(ssnow%smp/1000.0, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilMatPot = out%SoilMatPot / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilMatPot, 'SoilMatPot', &
               out%SoilMatPot, ranges%SoilMatPot, patchout%SoilMatPot, 'soil', met)
          ! Reset temporary output variable:
          out%SoilMatPot = 0.0
       END IF
    END IF     
    ! equilibrium soil matric potential
    IF((output%soil .or. output%EqSoilMatPot)  .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'EQ smp'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%EqSoilMatPot = out%EqSoilMatPot + REAL(ssnow%zq/1000.0, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%EqSoilMatPot = out%EqSoilMatPot / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%EqSoilMatPot, 'EqSoilMatPot', &
               out%EqSoilMatPot, ranges%EqSoilMatPot, patchout%EqSoilMatPot, 'soil', met)
          ! Reset temporary output variable:
          out%EqSoilMatPot = 0.0
       END IF
    END IF  
    ! equilibrium soil water content
    IF((output%soil .or. output%EqSoilMoist)  .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'EQ soilmoist'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%EqSoilMoist = out%EqSoilMoist + REAL(ssnow%wbeq, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%EqSoilMoist = out%EqSoilMoist / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%EqSoilMoist, 'EqSoilMoist', &
               out%EqSoilMoist, ranges%EqSoilMoist, patchout%EqSoilMoist, 'soil', met)
          ! Reset temporary output variable:
          out%EqSoilMoist = 0.0
       END IF
    END IF      
    ! infiltration rate
    IF(output%soil .OR. output%Qinfl) THEN
       !write(*,*) 'Qinfl'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%Qinfl = out%Qinfl + REAL(ssnow%fwtop/dels, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qinfl = out%Qinfl / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qinfl, 'Qinfl', &
               out%Qinfl, ranges%Qinfl, patchout%Qinfl, 'default', met)
          ! Reset temporary output variable:
          out%Qinfl = 0.0
       END IF
    END IF      
    ! infiltration rate
    IF((output%soil .OR. output%SatFrac)  .and. cable_user%GW_MODEL) THEN
       !write(*,*) 'Qinfl'    !MDeck
       ! Add current timestep's value to total of temporary output variable:
       out%SatFrac = out%SatFrac + REAL(ssnow%satfrac, 4)
       IF(writenow) THEN
          out%SatFrac = out%SatFrac / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SatFrac, 'SatFrac', &
               out%SatFrac, ranges%SatFrac, patchout%SatFrac, 'default', met)
          ! Reset temporary output variable:
          out%SatFrac = 0.0
       END IF
    END IF      

    ! recharge rate
    IF(output%soil .OR. output%Qrecharge) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%Qrecharge = out%Qrecharge + REAL(ssnow%Qrecharge, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%Qrecharge = out%Qrecharge / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%Qrecharge, 'Qrecharge', &
                       out%Qrecharge, ranges%Qrecharge, patchout%Qrecharge, 'default', met)
          ! Reset temporary output variable:
          out%Qrecharge = 0.0
       END IF
    END IF   

     ! SoilIce: av.layer soil moisture [kg/m^2]
    IF(output%soil .OR. output%SoilIce) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%SoilIce = out%SoilIce + REAL(ssnow%wbice, 4)
       IF(writenow) THEN
          ! Divide accumulated variable by number of accumulated time steps:
          out%SoilIce = out%SoilIce / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%SoilIce, 'SoilIce', &
               out%SoilIce, ranges%SoilMoist, patchout%SoilIce, 'soil', met)
          ! Reset temporary output variable:
          out%SoilIce = 0.0
       END IF
    END IF
   
    IF(output%rtevap) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%rtevap = out%rtevap + REAL(ssnow%rtevap_unsat, 4)
       IF(writenow) THEN
          out%rtevap = out%rtevap / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%rtevap, 'rtevap', &
               out%rtevap, (/-9999.0,9999.0/), patchout%rtevap, 'default', met)
          ! Reset temporary output variable:
          out%rtevap = 0.0
       END IF

       ! Add current timestep's value to total of temporary output variable:
       out%sv_rtevap = out%sv_rtevap + REAL(ssnow%sv_rtevap, 4)
       IF(writenow) THEN
          out%sv_rtevap = out%sv_rtevap / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%sv_rtevap, 'sv_rtevap', &
               out%sv_rtevap, (/-9999.0,9999.0/), patchout%sv_rtevap, 'default', met)
          ! Reset temporary output variable:
          out%sv_rtevap = 0.0
       END IF

       out%bl_rtevap = out%bl_rtevap + REAL(ssnow%bl_rtevap, 4)
       IF(writenow) THEN
          out%bl_rtevap = out%bl_rtevap / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%bl_rtevap, 'bl_rtevap', &
               out%bl_rtevap, (/-9999.0,9999.0/), patchout%bl_rtevap, 'default', met)
          ! Reset temporary output variable:
          out%bl_rtevap = 0.0
       END IF

    END IF      

    IF(output%rtevap_sat) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%rtevap_sat = out%rtevap_sat + REAL(ssnow%rtevap_sat, 4)
       IF(writenow) THEN
          out%rtevap_sat = out%rtevap_sat / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%rtevap_sat, 'rtevap_sat', &
               out%rtevap_sat, (/-9999.0,9999.0/), patchout%rtevap_sat, 'default', met)
          ! Reset temporary output variable:
          out%rtevap_sat = 0.0
       END IF
    END IF      
    IF(output%sublayer_dz) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%sublayer_dz = out%sublayer_dz + REAL(canopy%sublayer_dz, 4)
       IF(writenow) THEN
          out%sublayer_dz = out%sublayer_dz / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%sublayer_dz, 'sublayer_dz', &
               out%sublayer_dz, (/-9999.0,9999.0/), patchout%sublayer_dz, 'default', met)
          ! Reset temporary output variable:
          out%sublayer_dz = 0.0
       END IF
    END IF      

    IF(output%z0soil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%z0soil = out%z0soil + REAL(rough%z0soilsn, 4)
       IF(writenow) THEN
          out%z0soil = out%z0soil / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%z0soil, 'z0soil', &
               out%z0soil, (/-9999.0,9999.0/), patchout%z0soil, 'default', met)
          ! Reset temporary output variable:
          out%z0soil = 0.0
       END IF
    END IF      

    IF(output%rtsoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%rtsoil = out%rtsoil + REAL(ssnow%rtsoil, 4)
       IF(writenow) THEN
          out%rtsoil = out%rtsoil / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%rtsoil, 'rtsoil', &
               out%rtsoil, (/-9999.0,9999.0/), patchout%rtsoil, 'default', met)
          ! Reset temporary output variable:
          out%rtsoil = 0.0
       END IF
    END IF      

    IF(output%beta_rtsoil) THEN
       ! Add current timestep's value to total of temporary output variable:
       out%beta_rtsoil = out%beta_rtsoil + REAL(ssnow%beta_rtsoil, 4)
       IF(writenow) THEN
          out%beta_rtsoil = out%beta_rtsoil / REAL(output%interval, 4)
          ! Write value to file:
          CALL write_ovar(out_timestep, ncid_out, ovid%beta_rtsoil, 'beta_rtsoil', &
               out%beta_rtsoil, (/-9999.0,9999.0/), patchout%beta_rtsoil, 'default', met)
          ! Reset temporary output variable:
          out%beta_rtsoil = 0.0
       END IF
    END IF      

    !write(*,*) ' at end of write_output '    !MDeck
 
  END SUBROUTINE write_output
  !=============================================================================
  SUBROUTINE close_output_file(bal, air, bgc, canopy, met,                     &
                               rad, rough, soil, ssnow, sum_flux, veg)
    ! Closes output file, reports cumulative mass and energy 
    ! balances, and deallocates variables.
    TYPE (met_type), INTENT(INOUT)       :: met
    TYPE (air_type), INTENT(INOUT)       :: air
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE (veg_parameter_type), INTENT(INOUT) :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (radiation_type),INTENT(INOUT)  :: rad
    TYPE (sum_flux_type), INTENT(INOUT)  :: sum_flux
    TYPE(balances_type),INTENT(INOUT)    :: bal

    INTEGER :: i ! do loop counter

    ! Close file
    ok = NF90_CLOSE(ncid_out)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing output file '        &
                        //TRIM(filename%out)// '(SUBROUTINE close_output_file)')

    ! Report balance info to log file if verbose writing is requested:
    IF(output%balances .AND. verbose) THEN
       WRITE(logn, *)
       DO i = 1, mland
          WRITE(logn, '(A51,I7,1X,A11,E12.4,A6)')                              &
                ' Cumulative energy balance for each patch in site #',         &
                i,'is (W/m^2):'
          WRITE(logn, *)                                                       &
                bal%ebal_tot(landpt(i)%cstart:landpt(i)%cstart +               &
                landpt(i)%nap - 1)
          WRITE(logn,'(A50,I7,1X,A8,E12.4,A3)')                                &
                ' Cumulative water balance for each patch in site #',          &
                i,'is (mm):'
          WRITE(logn, *)                                                       &
                bal%wbal_tot(landpt(i)%cstart:landpt(i)%cstart +               &
                landpt(i)%nap - 1)
          WRITE(logn, *)
       END DO
    END IF

    ! Successful run!
    WRITE(logn, *)
    WRITE(logn, *) 'Run finished and output file closed.'

  END SUBROUTINE close_output_file
  !=============================================================================
  SUBROUTINE create_restart(logn, dels, ktau, soil, veg, ssnow,                      &
                            canopy, rough, rad, bgc, bal)
    ! Creates a restart file for CABLE using a land only grid with mland
    ! land points and max_vegpatches veg/soil patches (some of which may
    ! not be active). It uses CABLE's internal variable names.
    INTEGER, INTENT(IN) :: logn ! log file number
    REAL, INTENT(IN) :: dels ! time step size
    INTEGER, INTENT(IN)           :: ktau ! timestep number in loop which include spinup 
    TYPE (soil_parameter_type), INTENT(IN) :: soil ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)  :: veg  ! vegetation parameters
    TYPE (soil_snow_type), INTENT(IN)      :: ssnow  ! soil and snow variables
    TYPE (bgc_pool_type), INTENT(IN)       :: bgc    ! carbon pool variables
    TYPE (canopy_type), INTENT(IN)         :: canopy ! vegetation variables
    TYPE (roughness_type), INTENT(IN)      :: rough  ! roughness varibles
    TYPE (radiation_type), INTENT(IN)  :: rad ! radiation variables
    TYPE (balances_type), INTENT(IN) :: bal ! energy and water balance variables
!    INTEGER, INTENT(IN) :: mvtype
!    INTEGER, INTENT(IN) :: mstype

    TYPE(parID_type) :: rpid ! parameter IDs for restart nc file
    INTEGER :: ncid_restart ! netcdf restart file ID
    ! REAL, POINTER,DIMENSION(:,:) :: surffrac ! fraction of each surf type
    INTEGER :: dummy ! dummy argument in subroutine call
    INTEGER :: mlandID, mpID, radID, soilID, napID,                       &
                    soilcarbID, plantcarbID, tID, snowID ! dimension IDs
    !    INTEGER :: mlandID, surftypeID, patchID, radID, soilID, &
    !         soilcarbID, plantcarbID, tID, snowID ! dimension IDs
    INTEGER :: tvarID, latID, lonID !,surffracID ! time,lat,lon variable ID
    INTEGER :: tggID, wbID, wbiceID, tssID, ssdnnID, ssdnID, osnowdID,    &
                    smassID, sdepthID, snageID, snowdID, rtsoilID, isflagID,   &
                    canstoID, albsoilsnID, gammzzID, tggsnID, sghfluxID,       &
                    ghfluxID, runoffID, rnof1ID, rnof2ID, gaID, dgdtgID,       &
                    fevID, fesID, fhsID, wbtot0ID, osnowd0ID, cplantID,        &
                    csoilID, tradID, albedoID, gwID
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp netcdf file
    dummy = 0 ! initialise

    WRITE(logn, '(A24)') ' Writing restart file...'
    ! Create output file:
    ok = NF90_CREATE(filename%restart_out, NF90_CLOBBER, ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file '      &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')
    ! Put the file in define mode:
    ok = NF90_REDEF(ncid_restart)
    ! Define dimensions:
    ok = NF90_DEF_DIM(ncid_restart, 'mland', mland, mlandID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                     (ok, 'Error defining mland dimension in restart file. '// &
                      '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'mp_patch', mp, mpID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                  (ok, 'Error defining mp_patch dimension in restart file. '// &
                   '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'soil', ms, soilID) ! number of soil layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
             (ok, 'Error defining vertical soil dimension in restart file. '// &
              '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'snow', 3, snowID) ! number of snow layers
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
             (ok, 'Error defining vertical snow dimension in restart file. '// &
              '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'rad', nrb, radID) ! number of rad. bands
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                 (ok, 'Error defining radiation dimension in restart file. '// &
                  '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'soil_carbon_pools', ncs, soilcarbID)
    ! number of soil carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
          (ok, 'Error defining soil carbon pool dimension in restart file. '// &
           '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'plant_carbon_pools', ncp, plantcarbID)
    ! number of plant carbon pools
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
         (ok, 'Error defining plant carbon pool dimension in restart file. '// &
          '(SUBROUTINE create_restart)')
    ok = NF90_DEF_DIM(ncid_restart, 'time', 1, tID)
    IF (ok /= NF90_NOERR) CALL nc_abort &
                      (ok, 'Error defining time dimension in restart file. '// &
                       '(SUBROUTINE create_restart)')
    ! Define "time" variable and its attributes:
    ok=NF90_DEF_VAR(ncid_restart,'time',NF90_DOUBLE,(/tID/),tvarID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                       (ok, 'Error defining time variable in restart file. '// &
                        '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'units', timeunits)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
            (ok, 'Error defining time variable attributes in restart file. '// &
             '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, tvarID, 'coordinate', time_coord)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
            (ok, 'Error defining time variable attributes in restart file. '// &
             '(SUBROUTINE create_restart)')
    ! Define latitude and longitude variable:
    ok=NF90_DEF_VAR(ncid_restart, 'latitude', NF90_FLOAT, (/mlandID/), latID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                   (ok, 'Error defining latitude variable in restart file. '// &
                    '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart,latID,'units','degrees_north')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
        (ok, 'Error defining latitude variable attributes in restart file. '// &
         '(SUBROUTINE create_restart)')
    ok=NF90_DEF_VAR(ncid_restart, 'longitude', NF90_FLOAT, (/mlandID/), lonID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                  (ok, 'Error defining longitude variable in restart file. '// &
                   '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, lonID, 'units', 'degrees_east')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
       (ok, 'Error defining longitude variable attributes in restart file. '// &
        '(SUBROUTINE create_restart)')
    ! Define number of active patches variable:
    ok = NF90_DEF_VAR(ncid_restart, 'nap', NF90_FLOAT, (/mlandID/), napID)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                        (ok, 'Error defining nap variable in restart file. '// &
                         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, napID, 'long_name',                        &
                      'Number of active patches')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
              (ok,'Error defining nap variable attributes in restart file. '// &
               '(SUBROUTINE create_restart)')
    ! Define patch fraction variable:
    ok=NF90_DEF_VAR(ncid_restart, 'patchfrac', NF90_FLOAT, (/mpID/),           &
                    rpid%patchfrac)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                  (ok, 'Error defining patchfrac variable in restart file. '// &
                   '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%patchfrac, 'long_name',               &
                      'Fraction of vegetated grid cell area occupied by a '//  &
                      'vegetation/soil patch')
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
       (ok, 'Error defining patchfrac variable attributes in restart file. '// &
        '(SUBROUTINE create_restart)')
    ! mvtype (Number of vegetation types):
    ok = NF90_DEF_VAR(ncid_restart, 'mvtype', NF90_INT, rpid%mvtype)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                     (ok, 'Error defining mvtype variable in restart file. '// &
                      '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%mvtype, "long_name",                  &
                      "Number of vegetation types")
    ! mstype (Number of soil types):
    ok = NF90_DEF_VAR(ncid_restart, 'mstype', NF90_INT, rpid%mstype)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                     (ok, 'Error defining mstype variable in restart file. '// &
                      '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%mstype, "long_name",                  &
                      "Number of soil types")

    !======begin defining state variables=======================================
    ! Interface arguments: netcdf file ID, variableID, variable name, variable 
    ! units, variable long name, YES to write patch info (as this is a restart
    ! file), OPTIONAL extra dimension ID (e.g. for soil dimensioned variables),
    ! dimension switch to indicate what extra dimension is real or integer for
    ! single dim variables, xdimID,ydimID, zdimID (all three not used here),
    ! land dim ID, patch dim ID, YES we're writing a restart file.
    !------------------define soil states---------------------------------------
    CALL define_ovar(ncid_restart, tggID, 'tgg', 'K',                          &
                     'Average layer soil temperature',                         &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wbID, 'wb', 'vol/vol',                      &
                     'Average layer volumetric soil moisture',                 &
                     .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, wbiceID, 'wbice', 'vol/vol',                &
                     'Average layer volumetric soil ice',                      &
                     .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, tssID, 'tss', 'K',                          &
                     'Combined soil/snow temperature',                         &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, albsoilsnID, 'albsoilsn', '-',              &
                     'Combined soil/snow albedo',                              &
                     .TRUE., radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rtsoilID, 'rtsoil', '??',                   &
                     'Turbulent resistance for soil', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gammzzID, 'gammzz', 'J/kg/C',               &
                     'Heat capacity for each soil layer',                      &
                     .TRUE., soilID, 'r2soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, runoffID, 'runoff', 'mm/timestep',          &
                   'Total runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rnof1ID, 'rnof1', 'mm/timestep',            &
                 'Surface runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rnof2ID, 'rnof2', 'mm/timestep',            &
              'Subsurface runoff', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)


    CALL define_ovar(ncid_restart, gwID, 'GWwb', 'mm3/mm3','GW water content',&
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)



    !---------------define snow states------------------------------------------
    CALL define_ovar(ncid_restart, tggsnID, 'tggsn', 'K',                      &
                     'Average layer snow temperature',                         &
                     .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ssdnnID, 'ssdnn', 'kg/m^3',                 &
           'Average snow density', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ssdnID, 'ssdn', 'kg/m^3',                   &
                     'Average layer snow density',                             &
                     .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, snowdID, 'snowd', 'mm',                     &
                     'Liquid water eqivalent snow depth',                      &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, snageID, 'snage', '??',                     &
                     'Snow age', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, smassID, 'smass', 'kg/m^2',                 &
                     'Average layer snow mass',                                &
                     .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, sdepthID, 'sdepth', 'm',                    &
       'Snow layer depth', .TRUE., snowID, 'snow', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, osnowdID, 'osnowd', 'mm',                   &
                     'Previous time step snow depth in water equivalent',      &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, isflagID, 'isflag', '-',                    &
      'Snow layer scheme flag', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    !----------------define canopy states----------------------------------
    CALL define_ovar(ncid_restart, canstoID, 'cansto', 'mm',                   &
                     'Canopy surface water storage', .TRUE., 'real', 0, 0, 0,  &
                     mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, ghfluxID, 'ghflux', 'W/m^2?',               &
                     '????', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, sghfluxID, 'sghflux', 'W/m^2?',             &
                     '????', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, gaID, 'ga', 'W/m^2',                        &
               'Ground heat flux', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, dgdtgID, 'dgdtg', 'W/m^2/K',                &
                'Derivative of ground heat flux wrt soil temperature', .TRUE., &
                     'r2', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fevID, 'fev', 'W/m^2',                      &
                     'Latent heat flux from vegetation', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fesID, 'fes', 'W/m^2',                      &
                     'Latent heat flux from soil',                             &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, fhsID, 'fhs', 'W/m^2',                      &
                     'Sensible heat flux from soil',                           &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !--------------biogeochemical variables------------------------
    CALL define_ovar(ncid_restart, cplantID, 'cplant', 'gC/m^2',               &
                     'Plant carbon stores',                                    &
               .TRUE., plantcarbID, 'plantcarbon', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, csoilID, 'csoil', 'gC/m^2',                 &
                     'Soil carbon stores',                                     &
                 .TRUE., soilcarbID, 'soilcarbon', 0, 0, 0, mpID, dummy, .TRUE.)
    !-------------------others---------------------------------
    CALL define_ovar(ncid_restart, wbtot0ID, 'wbtot0', 'mm',                   &
                     'Initial time step soil water total',                     &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, osnowd0ID, 'osnowd0', 'mm',                 &
                     'Initial time step snow water total',                     &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, albedoID, 'albedo', '-',                    &
                     'Albedo for shortwave and NIR radiation',                 &
                     .TRUE., radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, tradID, 'trad', 'K',                        &
                    'Surface radiative temperature (soil/snow/veg inclusive)', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    !---------------------MODEL PARAMETERS---------------------------------   
    WRITE(logn,'(A43)') '   Writing model parameters to restart file'
    CALL define_ovar(ncid_restart, rpid%iveg, 'iveg', '-',                     &
             'Vegetation type', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%isoil, 'isoil', '-',                   &
                   'Soil type', .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%clay, 'clay', '-',                     &
                     'Fraction of soil which is clay',                         &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%sand, 'sand', '-',                     &
                     'Fraction of soil which is sand',                         &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%silt, 'silt', '-',                     &
                     'Fraction of soil which is silt',                         &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%ssat, 'ssat', '-',                     &
                     'Fraction of soil volume which is water @ saturation',    &
                    .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%sfc, 'sfc', '-',                       &
                    'Fraction of soil volume which is water @ field capacity', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%swilt, 'swilt', '-',                   &
                     'Fraction of soil volume which is water @ wilting point', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    ! zse (depth of each soil layer):
    ok = NF90_DEF_VAR(ncid_restart, 'zse', NF90_FLOAT, (/soilID/), rpid%zse)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                        (ok, 'Error defining zse variable in restart file. '// &
                         '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%zse, "long_name",                     &
                      "Depth of each soil layer")
    ok = NF90_PUT_ATT(ncid_restart, rpid%zse, "units", "m")
    CALL define_ovar(ncid_restart, rpid%froot, 'froot', '-',                   &
                     'Fraction of roots in each soil layer',                   &
                      .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%bch, 'bch', '-',                       &
                     'Parameter b, Campbell eqn 1985',                         &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%hyds, 'hyds', 'm/s',                   &
                     'Hydraulic conductivity @ saturation',                    &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%sucs, 'sucs', 'm',                     &
                     'Suction @ saturation', .TRUE.,                           &
                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%css, 'css', 'J/kg/C',                  &
                     'Heat capacity of soil minerals',                         &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%rhosoil, 'rhosoil', 'kg/m^3',          &
                     'Density of soil minerals',                               &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%rs20, 'rs20', '-',                     &
                     'Soil respiration coefficient at 20C',                    &
                      .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%albsoil, 'albsoil', '-',               &
                     'Soil reflectance', .TRUE.,                               &
                     radID, 'radiation', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%hc, 'hc', 'm',                         &
                     'Height of canopy', .TRUE.,                               &
                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%canst1, 'canst1', 'mm/LAI',            &
                     'Max water intercepted by canopy',                        &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%dleaf, 'dleaf', 'm',                   &
                     'Chararacteristic length of leaf',                        &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%frac4, 'frac4', '-',                   &
                     'Fraction of plants which are C4',                        &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%ejmax, 'ejmax', 'mol/m^2/s',           &
                     'Max potential electron transport rate top leaf', .TRUE., &
                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%vcmax, 'vcmax', 'mol/m^2/s',           &
                     'Maximum RuBP carboxylation rate top leaf', .TRUE.,       &
                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%rp20, 'rp20', '-',                     &
                     'Plant respiration coefficient at 20C', .TRUE., 'real',   &
                     0, 0, 0, mpID, dummy, .TRUE.)

    CALL define_ovar(ncid_restart, rpid%g0c3, 'g0c3', '-',                     &
                     'g0c3 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
                     0, 0, 0, mpID, dummy, .TRUE.) ! Ticket #56
    CALL define_ovar(ncid_restart, rpid%g0c4, 'g0c4', '-',                     &
                     'g0c4 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
                     0, 0, 0, mpID, dummy, .TRUE.)  ! Ticket #56
    CALL define_ovar(ncid_restart, rpid%g1c3, 'g1c3', '-',                     &
                     'g1c3 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
                     0, 0, 0, mpID, dummy, .TRUE.)  ! Ticket #56
    CALL define_ovar(ncid_restart, rpid%g1c4, 'g1c4', '-',                     &
                     'g1c4 term in Medlyn Stomatal Cond. Param', .TRUE.,'real',&
                     0, 0, 0, mpID, dummy, .TRUE.) ! Ticket #56

    CALL define_ovar(ncid_restart, rpid%rpcoef, 'rpcoef', '1/C',               &
                     'Temperature coef nonleaf plant respiration', .TRUE.,     &
                     'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%shelrb, 'shelrb', '-',                 &
              'Sheltering factor', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%xfang, 'xfang', '-',                   &
           'Leaf angle parameter', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%wai, 'wai', '-',                       &
                'Wood area index', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%vegcf, 'vegcf', '-',                   &
                     'vegcf', .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%extkn, 'extkn', '-',                   &
                     'Extinction coef for vertical nitrogen profile',          &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%tminvj, 'tminvj', 'C',                 &
                     'Min temperature for the start of photosynthesis',        &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%tmaxvj, 'tmaxvj', 'C',                 &
                     'Max temperature for the start of photosynthesis',        &
                      .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%vbeta, 'vbeta', '-',                   &
                     'Stomatal sensitivity to soil water',                     &
                      .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%xalbnir, 'xalbnir', '-',               &
                     'modifier for albedo in near ir band',                    &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
                     
                     
                     
                     
    !MD
    CALL define_ovar(ncid_restart, rpid%WatSat, &
                     'WatSat', '-', 'Max water content in soil layer', &
                      .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)       
    !CALL define_ovar(ncid_restart, rpid%GWWatSat, &
    !                 'GWWatSat', '-', 'Max water content in aquifer', &
    !                 .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)      

    CALL define_ovar(ncid_restart, rpid%GWWatSat, 'GWWatSat', 'mm3mm-3', &
                     'Aquifer Maximum porosity', .TRUE.,                 &
                     'real', 0, 0, 0, mpID, dummy, .TRUE.)


    CALL define_ovar(ncid_restart, rpid%Watr, &
                     'Watr', '-', 'residual water content in soil layer', &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
                 
    CALL define_ovar(ncid_restart, rpid%GWWatr, &
                     'GWWatr', '-', 'residual water content in aquifer', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)                  
                 
    CALL define_ovar(ncid_restart, rpid%SoilMatPotSat, &
                     'SoilMatPotSat', '-', 'soil matric potent at saturation content in soil layer', &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
                 
    CALL define_ovar(ncid_restart, rpid%GWSoilMatPotSat, &
                     'GWSoilMatPotSat', '-', 'soil matric potent at saturation content in aquifer', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)  
                 
    CALL define_ovar(ncid_restart, rpid%HkSat, &
                     'HkSat', '-', 'Max hydraulic conductivity in soil layer', &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)
                 
    CALL define_ovar(ncid_restart, rpid%GWHkSat, &
                     'GWHkSat', '-', 'Max hydraulic conductivityin aquifer', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
                 
    CALL define_ovar(ncid_restart, rpid%FrcSand, &
                     'FrcSand', '-', 'sand fraction  in soil layer', &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)    
                 
    CALL define_ovar(ncid_restart, rpid%FrcClay, &
                     'FrcClay', '-', 'clay fraction  in soil layer', &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)  
                 
    CALL define_ovar(ncid_restart, rpid%ClappB, &
                     'ClappB', '-', 'clapp and horn b param  in soil layer', &
                     .TRUE., soilID, 'soil', 0, 0, 0, mpID, dummy, .TRUE.)          
                     
                     
    ! ratecp (Plant carbon rate constant):
    ok = NF90_DEF_VAR(ncid_restart, 'ratecp', NF90_FLOAT, (/plantcarbID/),     &
                      rpid%ratecp)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                     (ok, 'Error defining ratecp variable in restart file. '// &
                      '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecp, "long_name",                  &
                      "Plant carbon rate constant")
    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecp, "units", "1/year")
    ! ratecs (Soil carbon rate constant):
    ok = NF90_DEF_VAR(ncid_restart, 'ratecs', NF90_FLOAT, (/soilcarbID/),      &
                      rpid%ratecs)
    IF (ok /= NF90_NOERR) CALL nc_abort                                        &
                     (ok, 'Error defining ratecs variable in restart file. '// &
                      '(SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecs, "long_name",                  &
                      "Soil carbon rate constant")
    ok = NF90_PUT_ATT(ncid_restart, rpid%ratecs, "units", "1/year")
    CALL define_ovar(ncid_restart, rpid%meth, 'meth', '-',                     &
                     'Canopy turbulence parameterisation switch',              &
                     .TRUE., 'integer', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%za_uv, 'za_uv', 'm',                   &
                    'Reference height (lowest atm. model layer) for momentum', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
    CALL define_ovar(ncid_restart, rpid%za_tq, 'za_tq', 'm',                   &
                     'Reference height (lowest atm. model layer) for scalars', &
                     .TRUE., 'real', 0, 0, 0, mpID, dummy, .TRUE.)
                     
    ! Write global attributes for file:
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate = todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime = nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Production",                 &
                      TRIM(todaydate)//' at '//TRIM(nowtime))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
                  //TRIM(filename%restart_out)// ' (SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "Source",                     &
                      'CABLE LSM restart file')
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
         //TRIM(filename%restart_out)// ' (SUBROUTINE create_restart)')
    ok = NF90_PUT_ATT(ncid_restart, NF90_GLOBAL, "CABLE_input_file",           &
                      TRIM(filename%met))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing global detail to '   &
                  //TRIM(filename%restart_out)// ' (SUBROUTINE create_restart)')

    ! End netcdf define mode:
    ok = NF90_ENDDEF(ncid_restart)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error creating restart file '      &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')

    ! Write time variable:
    ok = NF90_PUT_VAR(ncid_restart, tvarID, DBLE(REAL(ktau) * dels))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error time variable to '           &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')

    ! Write latitude and longitude variables:
    ok = NF90_PUT_VAR(ncid_restart, latID, latitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing latitude variable to '   &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, lonID, longitude)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing longitude variable to '  &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')

    ! Write number of active patches for each land grid cell:
    ok = NF90_PUT_VAR(ncid_restart, napID, landpt(:)%nap)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing nap variable to '        &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')

    ! Write vegetated patch fractions
    ok = NF90_PUT_VAR(ncid_restart, rpid%patchfrac,                            &
                      patch(:)%frac, start = (/1/), count = (/mp/))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing patchfrac to '       &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')

    ! Write number of veg and soil types
    ok = NF90_PUT_VAR(ncid_restart, rpid%mvtype,mvtype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing mvtype parameter to '    &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rpid%mstype,mstype)
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing mstype parameter to '    &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')

    ! Write parameters:
    CALL write_ovar (ncid_restart, rpid%iveg, 'iveg', REAL(veg%iveg, 4),       &
                     ranges%iveg, .TRUE., 'integer', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%isoil, 'isoil', REAL(soil%isoilm, 4),  &
                     ranges%isoil, .TRUE., 'integer', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%bch, 'bch', REAL(soil%bch, 4),         &
                     ranges%bch, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%clay, 'clay', REAL(soil%clay, 4),      &
                     ranges%clay, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%sand, 'sand', REAL(soil%sand, 4),      &
                     ranges%sand, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%silt, 'silt', REAL(soil%silt, 4),      &
                     ranges%silt, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%css, 'css', REAL(soil%css, 4),         &
                     ranges%css, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%rhosoil, 'rhosoil'                    &
                     ,REAL(soil%rhosoil,4), ranges%rhosoil, .TRUE., 'real',     &
                     .TRUE.)
    CALL write_ovar (ncid_restart, rpid%hyds, 'hyds', REAL(soil%hyds, 4),      &
                     ranges%hyds, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%sucs, 'sucs', REAL(soil%sucs, 4),      &
                     ranges%sucs, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%rs20, 'rs20', REAL(veg%rs20, 4),       &
                     ranges%rs20, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%ssat, 'ssat', REAL(soil%ssat, 4),      &
                     ranges%ssat, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%sfc, 'sfc', REAL(soil%sfc, 4),         &
                     ranges%sfc, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%swilt, 'swilt', REAL(soil%swilt, 4),   &
                     ranges%swilt, .TRUE., 'real', .TRUE.)
    ! Soil dimensioned variables/parameters:
    CALL write_ovar (ncid_restart, rpid%froot, 'froot', REAL(veg%froot, 4),    &
                     ranges%froot, .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, tggID, 'tgg', REAL(ssnow%tgg, 4),           &
                     ranges%SoilTemp, .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, wbID, 'wb', ssnow%wb, ranges%SoilMoist,     &
                     .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, wbiceID, 'wbice', ssnow%wbice,              &
                     ranges%SoilMoist, .TRUE., 'soil', .TRUE.)
    CALL write_ovar (ncid_restart, gammzzID, 'gammzz', ssnow%gammzz,           &
                     (/-99999.0, 9999999.0/), .TRUE., 'soil', .TRUE.)
                     
   !MD
    CALL write_ovar (ncid_restart, rpid%WatSat, 'WatSat', REAL(soil%watsat, 4),    &
                     ranges%WatSat, .TRUE., 'soil', .TRUE.)   
    CALL write_ovar (ncid_restart, rpid%SoilMatPotSat, 'SoilMatPotSat', REAL(soil%smpsat, 4),    &
                     ranges%SoilMatPotSat, .TRUE., 'soil', .TRUE.)                
    CALL write_ovar (ncid_restart, rpid%HkSat, 'HkSat', REAL(soil%hksat, 4),    &
                     ranges%HkSat, .TRUE., 'soil', .TRUE.)   
    CALL write_ovar (ncid_restart, rpid%ClappB, 'ClappB', REAL(soil%clappB, 4),    &
                     ranges%ClappB, .TRUE., 'soil', .TRUE.)   
    CALL write_ovar (ncid_restart, rpid%FrcSand, 'FrcSand', REAL(soil%Fsand, 4),    &
                     ranges%FrcSand, .TRUE., 'soil', .TRUE.)                       
    CALL write_ovar (ncid_restart, rpid%FrcClay, 'FrcClay', REAL(soil%Fclay, 4),    &
                     ranges%FrcClay, .TRUE., 'soil', .TRUE.)                 

    !GW MD
    CALL write_ovar (ncid_restart, rpid%GWWatSat, 'GWWatSat', REAL(soil%GWwatsat, 4),    &
                     ranges%GWWatSat, .TRUE., 'real', .TRUE.)   
    CALL write_ovar (ncid_restart, rpid%GWSoilMatPotSat, 'GWSoilMatPotSat', REAL(soil%GWsmpsat, 4),    &
                     ranges%GWSoilMatPotSat, .TRUE., 'real', .TRUE.)                
    CALL write_ovar (ncid_restart, rpid%GWHkSat, 'GWHkSat', REAL(soil%GWhksat, 4),    &
                     ranges%GWHkSat, .TRUE., 'real', .TRUE.)       
                     
    ! Snow dimensioned variables/parameters:
    CALL write_ovar (ncid_restart, ssdnID, 'ssdn', REAL(ssnow%ssdn, 4),        &
                     (/0.0, 9999.0/), .TRUE., 'snow', .TRUE.)
    CALL write_ovar (ncid_restart, smassID, 'smass', REAL(ssnow%smass, 4),     &
                     (/0.0, 9999.0/), .TRUE., 'snow', .TRUE.)
    CALL write_ovar (ncid_restart, sdepthID, 'sdepth', REAL(ssnow%sdepth, 4),  &
                     (/0.0, 9999.0/), .TRUE., 'snow', .TRUE.)
    CALL write_ovar (ncid_restart, tggsnID, 'tggsn', REAL(ssnow%tggsn, 4),     &
                     (/100.0, 300.0/), .TRUE., 'snow', .TRUE.)
    ! Other dims
    CALL write_ovar (ncid_restart, albsoilsnID, 'albsoilsn',                   &
            REAL(ssnow%albsoilsn, 4), (/0.0, 1.0/), .TRUE., 'radiation', .TRUE.)
    CALL write_ovar (ncid_restart, cplantID, 'cplant', REAL(bgc%cplant, 4),    &
                     (/-99999.0, 9999999.0/), .TRUE., 'plantcarbon', .TRUE.)
    CALL write_ovar (ncid_restart, csoilID, 'csoil', REAL(bgc%csoil, 4),       &
                     (/-99999.0, 9999999.0/), .TRUE., 'soilcarbon', .TRUE.)
    ok = NF90_PUT_VAR(ncid_restart, rpid%zse, REAL(soil%zse, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok, 'Error writing zse parameter to '   &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')
    ! Single dim:

    !MD
    CALL write_ovar (ncid_restart, gwID, 'GWwb', REAL(ssnow%GWwb, 4),       &
                     ranges%GWwb, .TRUE., 'real', .TRUE.)

    CALL write_ovar (ncid_restart, rpid%albsoil, 'albsoil',                    &
                     REAL(soil%albsoil, 4), ranges%albsoil, .TRUE.,            &
                     'radiation', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%canst1, 'canst1', REAL(veg%canst1, 4), &
                     ranges%canst1, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%dleaf, 'dleaf', REAL(veg%dleaf, 4),    &
                     ranges%dleaf, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%ejmax, 'ejmax', REAL(veg%ejmax, 4),    &
                     ranges%ejmax, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%vcmax, 'vcmax', REAL(veg%vcmax, 4),    &
                     ranges%vcmax, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%frac4, 'frac4', REAL(veg%frac4, 4),    &
                     ranges%frac4, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%hc, 'hc', REAL(veg%hc, 4),             &
                     ranges%hc, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%rp20, 'rp20', REAL(veg%rp20, 4),       &
                     ranges%rp20, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%g0c3, 'g0c3', REAL(veg%g0c3, 4),       &
                     ranges%g0c3, .TRUE., 'real', .TRUE.) ! Ticket #56
    CALL write_ovar (ncid_restart, rpid%g0c4, 'g0c4', REAL(veg%g0c4, 4),       &
                     ranges%g0c4, .TRUE., 'real', .TRUE.) ! Ticket #56
    CALL write_ovar (ncid_restart, rpid%g1c3, 'g1c3', REAL(veg%g1c3, 4),       &
                     ranges%g1c3, .TRUE., 'real', .TRUE.) ! Ticket #56
    CALL write_ovar (ncid_restart, rpid%g1c4, 'g1c4', REAL(veg%g1c4, 4),       &
                     ranges%g1c4, .TRUE., 'real', .TRUE.) ! Ticket #56
    CALL write_ovar (ncid_restart, rpid%rpcoef, 'rpcoef', REAL(veg%rpcoef, 4), &
                     ranges%rpcoef, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%shelrb, 'shelrb', REAL(veg%shelrb, 4), &
                     ranges%shelrb, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%xfang, 'xfang', REAL(veg%xfang, 4),    &
                     ranges%xfang, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%wai, 'wai', REAL(veg%wai, 4),          &
                     ranges%wai, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%vegcf, 'vegcf', REAL(veg%vegcf, 4),    &
                     ranges%vegcf, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%extkn, 'extkn', REAL(veg%extkn, 4),    &
                     ranges%extkn, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%tminvj, 'tminvj', REAL(veg%tminvj, 4), &
                     ranges%tminvj, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%tmaxvj, 'tmaxvj', REAL(veg%tmaxvj, 4), &
                     ranges%tmaxvj, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%vbeta, 'vbeta', REAL(veg%vbeta, 4),    &
                     ranges%vbeta, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%xalbnir, 'xalbnir',                    &
                     REAL(veg%xalbnir, 4), ranges%xalbnir, .TRUE.,             &
                     'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%tmaxvj, 'tmaxvj', REAL(veg%tmaxvj, 4), &
                     ranges%tmaxvj, .TRUE., 'real', .TRUE.)
    ok = NF90_PUT_VAR(ncid_restart, rpid%ratecp, REAL(bgc%ratecp, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing ratecp parameter to '    &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')
    ok = NF90_PUT_VAR(ncid_restart, rpid%ratecs, REAL(bgc%ratecs, 4))
    IF(ok /= NF90_NOERR) CALL nc_abort(ok,                                     &
                                       'Error writing ratecs parameter to '    &
                   //TRIM(filename%restart_out)// '(SUBROUTINE create_restart)')
    CALL write_ovar (ncid_restart, rpid%meth, 'meth', REAL(veg%meth, 4),       &
                     ranges%meth, .TRUE., 'integer', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%za_uv, 'za_uv', REAL(rough%za_uv, 4),  &
                     ranges%za, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rpid%za_tq, 'za_tq', REAL(rough%za_tq, 4),  &
                     ranges%za, .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, tssID, 'tss', REAL(ssnow%tss, 4),           &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, ssdnnID, 'ssdnn', REAL(ssnow%ssdnn, 4),     &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, osnowdID, 'osnowd', REAL(ssnow%osnowd, 4),  &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, snageID, 'snage', REAL(ssnow%snage, 4),     &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, snowdID, 'snowd', REAL(ssnow%snowd, 4),     &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rtsoilID, 'rtsoil', REAL(ssnow%rtsoil, 4),  &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, isflagID, 'isflag', REAL(ssnow%isflag, 4),  &
                     (/-99999.0, 9999999.0/), .TRUE., 'integer', .TRUE.)
    CALL write_ovar (ncid_restart, canstoID, 'cansto', REAL(canopy%cansto, 4), &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, sghfluxID, 'sghflux',                       &
                     REAL(canopy%sghflux, 4), (/-99999.0, 9999999.0/),         &
                     .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, ghfluxID, 'ghflux', REAL(canopy%ghflux, 4), &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, runoffID, 'runoff', REAL(ssnow%runoff, 4),  &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rnof1ID, 'rnof1', REAL(ssnow%rnof1, 4),     &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, rnof2ID, 'rnof2', REAL(ssnow%rnof2, 4),     &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, gaID, 'ga', REAL(canopy%ga, 4),             &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, dgdtgID, 'dgdtg', canopy%dgdtg,             &
                     (/-99999.0, 9999999.0/), .TRUE., 'r2', .TRUE.)
    CALL write_ovar (ncid_restart, fevID, 'fev', REAL(canopy%fev, 4),          &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, fesID, 'fes', REAL(canopy%fes, 4),          &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, fhsID, 'fhs', REAL(canopy%fhs, 4),          &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, wbtot0ID, 'wbtot0', REAL(bal%wbtot0, 4),    &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, osnowd0ID, 'osnowd0', REAL(bal%osnowd0, 4), &
                     (/-99999.0, 9999999.0/), .TRUE., 'real', .TRUE.)
    CALL write_ovar (ncid_restart, albedoID, 'albedo', REAL(rad%albedo, 4),    &
                     ranges%Albedo, .TRUE., 'radiation', .TRUE.)
    CALL write_ovar (ncid_restart, tradID, 'trad',                             &
                     REAL(rad%trad, 4), ranges%RadT, .TRUE., 'real', .TRUE.)

    ! Close restart file
    ok = NF90_CLOSE(ncid_restart)

    WRITE(logn, '(A36)') '   Restart file complete and closed.'

  END SUBROUTINE create_restart


END MODULE cable_output_module
