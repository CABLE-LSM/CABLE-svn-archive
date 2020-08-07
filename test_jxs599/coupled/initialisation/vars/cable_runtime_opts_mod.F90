MODULE cable_runtime_opts_mod

IMPLICIT NONE

  ! hydraulic_redistribution parameters _soilsnow module
  REAL :: wiltParam=0.0, satuParam=0.0

  ! user switches turned on/off by the user thru namelists
  ! CABLE-2.0 user switches all in single namelist file cable.nml
  ! clean these up for new namelist(s) format
  TYPE kbl_user_switches
    !jhan:make this logical
     CHARACTER(LEN=3) :: DIAG_SOIL_RESP=''

     CHARACTER(LEN=20) :: FWSOIL_SWITCH=''

     ! Ticket #56
    !jhan:options?
     CHARACTER(LEN=20) :: GS_SWITCH=''

     !INH - new switch for revised coupling on implicit step of ACCESS-CM2 Ticket #132
     LOGICAL :: l_revised_coupling = .FALSE.

     !INH -apply revised sensitvity/correction terms to soilsnow energy balance
     LOGICAL :: L_REV_CORR = .FALSE.     !switch to revert to unchanged code

     !ticket#179
     LOGICAL :: soil_thermal_fix=.FALSE.

    !jhan:options?
     CHARACTER(LEN=3) :: SSNOW_POTEV=''
     
     
     
     
     
     !jhan: this is redundant now we all use filename%veg?
     CHARACTER(LEN=200) ::                                                    &
          VEG_PARS_FILE  !

     CHARACTER(LEN=20) ::                                                     &
          !H!FWSOIL_SWITCH, &     !
          PHENOLOGY_SWITCH = 'MODIS'   ! alternative is 'climate'
     !--- LN ------------------------------------------[

     CHARACTER(LEN=10) :: RunIden       = 'STANDARD'  !
     CHARACTER(LEN=6)  :: MetType       = ' ' !
     CHARACTER(LEN=20) :: SOIL_STRUC    = "default" ! 'default' or 'sli'
     CHARACTER(LEN=3)  :: POP_out       = 'rst' ! POP output type ('epi' or 'rst')
     CHARACTER(LEN=50) :: POP_rst       = ' ' !
     CHARACTER(LEN=8)  :: CASA_OUT_FREQ = 'annually' ! 'daily', 'monthly', 'annually'
     CHARACTER(LEN=10)  :: vcmax = 'standard' ! "standard" or "Walker2014"
     CHARACTER(LEN=10)  :: POPLUC_RunType = 'static' ! 'static', 'init', 'restart'

     LOGICAL ::                                                               &
          CALL_POP               = .FALSE., & !
          POP_fromZero           = .FALSE., &
          CALL_Climate           = .FALSE., &
          Climate_fromZero       = .FALSE., &
          CASA_fromZero          = .FALSE., &
          POPLUC                 = .FALSE.

     INTEGER  :: &
          CASA_SPIN_STARTYEAR = 1950, &
          CASA_SPIN_ENDYEAR   = 1960, &
          YEARSTART           = 0, &
          YEAREND             = 0, &
          CASA_NREP           = 1
     !--- LN ------------------------------------------]

     CHARACTER(LEN=5) ::                                                      &
          RUN_DIAG_LEVEL  !

     CHARACTER(LEN=3) ::                                                      &
          !H!DIAG_SOIL_RESP,   & ! either ON or OFF (jhan:Make Logical)
          LEAF_RESPIRATION    ! either ON or OFF (jhan:Make Logical)

     ! Custom soil respiration - see Ticket #42
     CHARACTER(LEN=10) ::                                                     &
          SMRF_NAME,   & ! Soil Moist Respiration Function
          STRF_NAME      ! Soil Temp Respiration Function

     LOGICAL ::                                                               &
          INITIALIZE_MAPPING    = .FALSE., & !
          CONSISTENCY_CHECK     = .FALSE., & !
          CASA_DUMP_READ        = .FALSE., & !
          CASA_DUMP_WRITE       = .FALSE., & !
          CABLE_RUNTIME_COUPLED = .TRUE. , & !
          LogWorker             = .TRUE. , & ! Write Output of each worker
                                ! L.Stevens - Test Switches
          L_NEW_ROUGHNESS_SOIL  = .FALSE., & !
          L_NEW_RUNOFF_SPEED    = .FALSE., & !
          L_NEW_REDUCE_SOILEVP  = .FALSE., & !

                                ! Switch for customized soil respiration - see Ticket #42
          SRF = .FALSE., &

                                !! vh_js !!
          litter = .FALSE.

     !MD
     LOGICAL :: GW_MODEL = .FALSE.
     LOGICAL :: alt_forcing = .FALSE.

     !using GSWP3 forcing?
     LOGICAL :: GSWP3 = .FALSE.
     LOGICAL :: or_evap = .FALSE.
     LOGICAL :: test_new_gw=.FALSE.
     LOGICAL :: sync_nc_file=.FALSE.
     INTEGER :: max_spins = -1
     LOGICAL :: fix_access_roots = .FALSE.  !use pft dependent roots in ACCESS
     !ACCESS roots
     LOGICAL :: access13roots = .FALSE.     !switch to use ACCESS1.3 %froot

     LOGICAL :: l_limit_labile = .FALSE.    ! #237: limit Labile in spinup
     LOGICAL :: NtilesThruMetFile = .FALSE. ! #199: Specify Ntiles thru met file 

  END TYPE kbl_user_switches

  ! instantiate internal switches
  TYPE(kbl_user_switches), SAVE :: cable_user



END MODULE cable_runtime_opts_mod

