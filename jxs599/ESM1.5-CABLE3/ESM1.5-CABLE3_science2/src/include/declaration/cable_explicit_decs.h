
      Integer ::                      &
      day   &
      ,OFF_X       & ! Size of small halo in i
      ,OFF_Y       & ! Size of small halo in j.
      ,ROW_LENGTH  & ! Local number of points on a row                
      ,ROWS        & ! Local number of rows in a theta field                  
      ,N_ROWS      & ! Local number of rows in a v field                      
      ,LAND_PTS    & ! IN No of land points being processed.                  
      ,LAND_PTS_TRIF               & 
      ,NTILES      & ! IN No. of land-surface tiles
      ,SM_LEVELS   &  ! IN No. of soil moisture levels
      ,NPFT_TRIF    
      
      INTEGER ::                     &                              
       LAND_INDEX(LAND_PTS)          ! IN LAND_INDEX(I)=J => the Jth         
!                                    ! point in ROW_LENGTH,ROWS is the    
!                                    ! land point.                        
                                                                       
      LOGICAL ::  l_cable
     
      LOGICAL,DIMENSION(LAND_PTS,NTILES) :: L_TILE_PTS
      LOGICAL ::                      &
       LAND_MASK(ROW_LENGTH,ROWS)   ! IN T if L_TILE_PTS, F elsewhere.

      REAL ::                         &                                   
      LW_DOWN(ROW_LENGTH,ROWS)    & ! IN Surface downward LW radiation (W/m2).
      ,surf_down_sw(row_length,rows,4)  &  ! IN Surface downward SW radiation
      ,cos_zenith_angle(row_length,rows)&
      ,sin_theta_latitude(row_length,rows)&
      ,SW_DOWN(ROW_LENGTH,ROWS)     & ! Surface downward SW radiation (W/m2).
      ,SW_TILE(LAND_PTS,NTILES)     &                                    
                                     ! IN Surface net SW radiation on
      ,TL_1(ROW_LENGTH,ROWS)       & ! IN Ice/liquid water temperature       
      ,LS_RAIN(ROW_LENGTH,ROWS)    & ! IN Large scale rain
      ,LS_SNOW(ROW_LENGTH,ROWS)    & ! IN Large scale snow
      ,CON_RAIN(ROW_LENGTH,ROWS)   & ! IN Convective rain
      ,CONV_SNOW(ROW_LENGTH,ROWS)  & ! IN Convective snow
                                    !    land tiles (W/m2).
!!!
      ,QW_1(ROW_LENGTH,ROWS)       & ! IN Total water content                
      ,PSTAR(ROW_LENGTH,ROWS)      & ! IN Surface pressure (Pascals).        
      ,VSHR_LAND(ROW_LENGTH,ROWS)  &
      ,LATITUDE(ROW_LENGTH,ROWS)   &
      ,LONGITUDE(ROW_LENGTH,ROWS)  &
      ,time_sec                    &   ! actual time of day in secs.
      ,Z1_UV(ROW_LENGTH,ROWS)      & ! IN Height of lowest uv level (m).     
      ,Z1_TQ(ROW_LENGTH,ROWS)       ! IN Height of lowest tq level (m).     
      

      REAL                 ::         &                                   
       albsoil(land_pts)            &
      ,CANOPY_TILE(LAND_PTS,NTILES) &  ! IN Surface/canopy water for           
!                                     !    snow-free land tiles (kg/m2)       
      ,CATCH(LAND_PTS,NTILES)       &  ! IN Surface/canopy water capacity      
!                                      !    of snow-free land tiles (kg/m2).   
      ,CATCH_SNOW(LAND_PTS)         &  ! IN Snow interception capacity of      
                                       !    NLT tile (kg/m2).                  
      ,HCON(LAND_PTS)               &  ! IN Soil thermal conductivity          
                                       !    (W/m/K).                           
      ,LAI_FT(LAND_PTS,NPFT)        &  ! IN Leaf area index
      ,TSOIL(LAND_PTS,SM_LEVELS)    &  ! IN Soil temperatures (K).
      ,SNOW_TILE(LAND_PTS,NTILES)   &  ! IN Lying snow on tiles (kg/m2)        
      ,SMVCCL(LAND_PTS)             &  ! IN Critical volumetric SMC            
                                       !    (cubic m per cubic m of soil).     
      ,SMVCST(LAND_PTS)             &  ! IN Volumetric saturation point        
                                       !    (m3/m3 of soil).                   
      ,SMVCWT(LAND_PTS)             &  ! IN Volumetric wilting point           
                                       !    (cubic m per cubic m of soil).     
      ,STHF(LAND_PTS,SM_LEVELS)     &  ! IN Frozen soil moisture content of    
                                       !    each layer as a fraction of        
                                       !    saturation.                        
      ,STHU(LAND_PTS,SM_LEVELS)     &  ! IN Unfrozen soil moisture content     
                                       !    of each layer as a fraction of     
                                       !    saturation.                        
      ,Z0_TILE(LAND_PTS,NTILES)   &  ! IN Tile roughness lengths (m).        
      ,FLAND(LAND_PTS)               ! IN Land fraction on land tiles.       

      INTEGER ::                       & 
      CO2_DIM_LEN &
                                   ! IN Length of a CO2 field row.
      ,CO2_DIM_ROW 

      REAL ::                        &                                
       TIMESTEP                     ! IN Timestep (seconds).                
      REAL ::                         &                                   
      CO2_MMR                          &
      ,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  
      REAL ::                        &                                
      FRAC(LAND_PTS,ntiles)       &  ! IN Fractions of surface types.        
      ,CANHT_FT(LAND_PTS,NPFT)    &  ! IN Canopy height (m)                  
      ,TSTAR_TILE(LAND_PTS,NTILES) &   ! IN Surface tile temperatures          
                                     ! actual time of day in secs.
      ,TSOIL_TILE(LAND_PTS,NTILES,SM_LEVELS)   & !
      ,SMCL(LAND_PTS,SM_LEVELS)                & !
      ,SMCL_TILE(LAND_PTS,NTILES,SM_LEVELS)    & !
      ,STHU_TILE(LAND_PTS,NTILES,SM_LEVELS)    &
      ,STHF_TILE(LAND_PTS,NTILES,SM_LEVELS)    &
      ,TSTAR_LAND(ROW_LENGTH,ROWS) &  ! IN Land mean surface temperature (K
      ,TSTAR(ROW_LENGTH,ROWS)       & ! IN GBM surface temperature (K).
      ,BEXP(LAND_PTS)                &      
      ,SATCON(LAND_PTS)              &
                                       !  qrparm.soil.satcon
      ,SATHH(LAND_PTS)               & !  soil water suction
      ,HCAP(LAND_PTS)                &
                                       ! soil/qrparm.soil.hcap
       ,CD(ROW_LENGTH,ROWS)          &  ! OUT Turbulent surface exchange        
                                       ! (bulk transfer) coefficient for   
                                       ! momentum.                         
      ,CH(ROW_LENGTH,ROWS)          &  ! OUT Turbulent surface exchange        
                                       ! (bulk transfer) coefficient for   
                                       ! heat and/or moisture.             
     &,CD_TILE(LAND_PTS,NTILES)     &
                                       ! Drag coefficient
     &,CH_TILE(LAND_PTS,NTILES)     &
                                       ! Transfer coefficient for heat and
                                       ! moisture

      ,FQW_1(ROW_LENGTH,ROWS)       &  ! OUT Moisture flux between layers   
                                       ! (kg per square metre per sec).    
                                       ! FQW(,1) is total water flux       
                                       ! from surface, 'E'.                
      ,FTL_1(ROW_LENGTH,ROWS)       &  ! OUT FTL(,K) contains net turbulent    
                                       ! sensible heat flux into layer K   
                                       ! from below; so FTL(,1) is the     
                                       ! surface sensible heat, H.(W/m2)   
      ,FTL_TILE(LAND_PTS,NTILES)    &  ! OUT Surface FTL for land tiles        
      ,LE_TILE(LAND_PTS,NTILES)     &  ! OUT Surface latent heat flux for      
                                       ! land tiles                        
                                       ! sea-ice (W/m2)                    
      ,RADNET_TILE(LAND_PTS,NTILES) &  ! OUT Surface net radiation 
      , U_S_STD_TILE(LAND_PTS,NTILES) &  ! OUT Surface friction velocity      
 
       ,U_S(ROW_LENGTH,ROWS)        & ! OUT Surface friction velocity (m/s)   
      ,T1_SD(ROW_LENGTH,ROWS)      & ! OUT Standard deviation of turbulent   
                                     !     fluctuations of layer 1 temp;     
                                     !     used in initiating convection.    
      ,Q1_SD(ROW_LENGTH,ROWS)      &  ! OUT Standard deviation of turbulent   
                                     !     flucs of layer 1 humidity;        
                                     !     used in initiating convection.    
        ,RESFS(LAND_PTS,NTILES)      & ! OUT Combined soil, stomatal
                                     !     and aerodynamic resistance
                                     !     factor for fraction (1-FRACA)
                                     !     of snow-free land tiles.
      ,RESFT(LAND_PTS,NTILES)      & ! OUT Total resistance factor.
                                     !     FRACA+(1-FRACA)*RESFS for
                                     !     snow-free L_TILE_PTS, 1 for snow.    
       ,FQW_TILE(LAND_PTS,NTILES)   & ! OUT Surface FQW for land tiles        
      ,FRACA(LAND_PTS,NTILES)      & ! OUT Fraction of surface moisture      
                                     !     flux with only aerodynamic        
                                     !     resistance for snow-free land     
                                     !     tiles.                            
     ,RHOSTAR(ROW_LENGTH,ROWS)    & ! OUT Surface air density               
     !H_BLEND_OROG(ROW_LENGTH,ROWS)&                                 
                                     ! OUT Blending height used as part of   
                                     !     effective roughness scheme        
      ,Z0H(ROW_LENGTH,ROWS)        & ! OUT Roughness length for heat and     

                                     !     moisture (m).                     
      ,Z0H_TILE(LAND_PTS,NTILES)   & ! OUT Tile roughness lengths for heat   
                                     !     and moisture (m).                 
      ,Z0M(ROW_LENGTH,ROWS)        & ! OUT Roughness length for              
                                     !     momentum (m).                     
      ,Z0M_TILE(LAND_PTS,NTILES)   & ! OUT Tile roughness lengths for        
                                     !     momentum.                         
      ,Z0M_EFF(ROW_LENGTH,ROWS)    & ! OUT Effective grid-box roughness      
      ,CANHC_TILE(LAND_PTS,NTILES) & ! IN Areal heat capacity of canopy      
                                     !    for land tiles (J/K/m2).           
      ,WT_EXT_TILE(LAND_PTS,SM_LEVELS,NTILES) &
                                     ! OUT Fraction of evapotranspiration     
                                     !    which is extracted from each       
                                     !    soil layer by each tile.           
      ,FLAKE(LAND_PTS,NTILES)      & ! IN Lake fraction.                  
      ,TILE_FRAC(LAND_PTS,NTILES)  & ! OUT Tile fractions including          
                                     !     snow cover in the ice tile.       
      ,VFRAC_TILE(LAND_PTS,NTILES)  & ! OUT Tile fractions including          
      ,FSMC(LAND_PTS,NPFT)         & ! OUT Moisture availability factor.
      ,SNOW_COND(LAND_PTS,NTILES,3)  &
       ! ,T_SURF(LAND_PTS)             &
      ,T_SURF_TILE(LAND_PTS,NTILES)  & 
      ,SNOW_RHO3L(LAND_PTS,NTILES,3) &                                 
                                       ! Snow density  (3 layer)              
      ,SNOW_RHO1L(LAND_PTS,NTILES)   &                         
                                       ! Mean snow density  (or 1 layer)
      ,SNAGE_TILE(LAND_PTS,NTILES)   &
      ,RTSOIL_TILE(LAND_PTS,NTILES)  &
      ,GFLUX_TILE(LAND_PTS,NTILES)   &
      ,SGFLUX_TILE(LAND_PTS,NTILES)  & 
      ,HCONS(LAND_PTS)                ! OUT Soil thermal conductivity including      

      INTEGER      ::                &
     & TILE_INDEX(LAND_PTS,ntiles)  & ! OUT Index of tile points              
     &,TILE_PTS(ntiles)               ! OUT Number of tile points             

      !In/outs :  
      REAL   ::                      &                                 
      RECIP_L_MO_TILE(LAND_PTS,NTILES) &
      ,EPOT_TILE(LAND_PTS,NTILES)      &
                                     ! Reciprocal of the Monin-Obukhov
                                     ! length for tiles (m^-1).
      ,CH_CAB(LAND_PTS)             &
      ,CD_CAB(LAND_PTS)             &
      ,U_S_CAB(LAND_PTS)            
 
      REAL ::                              &
      SNOW_DEPTH3L(LAND_PTS,NTILES,3)   &  ! Snow depth on tiles (3 layer) (m)
     ,SNOW_MASS3L(LAND_PTS,NTILES,3)    &  ! Snow mass on tiles (3 layer)
     ,SNOW_TMP3L(LAND_PTS,NTILES,3)     &   ! Snow temperature (3 layer)
      !jhan:not used.not aliased
      ,ZH(ROW_LENGTH,ROWS)           ! IN Height above surface of top of     

      Real ::                              &
       FTL_TILE_CAB(LAND_PTS,NTILES)    &
      ,FTL_CAB(LAND_PTS)                &
      ,LE_TILE_CAB(LAND_PTS,NTILES)     &
      ,LE_CAB(LAND_PTS)                 &
      ,TSTAR_TILE_CAB(LAND_PTS,NTILES)  &
      ,TSTAR_CAB(LAND_PTS)              &
      ,SMCL_CAB(LAND_PTS,SM_LEVELS)     &
      ,TSOIL_CAB(LAND_PTS,SM_LEVELS)    &
      ,USTAR_CAB(LAND_PTS)              &
      ,SURF_HTF_CAB(LAND_PTS)           &  
      ,F_ROOT(SM_LEVELS)            

      INTEGER ::                           &                                
      ISNOW_FLG3L(LAND_PTS,NTILES)     &
      , SOIL_TYPE(row_length,rows)      &
      , VEG_TYPE(row_length,rows), ntau   
      
      
      real :: HO2R2_OROG(LAND_PTS)
      logical :: L_CO2_INTERACTIVE   
