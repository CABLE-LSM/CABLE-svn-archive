
         ! parameter list from SF_EXCH
                                                                        
       ! IN atmospheric forcing
       surf_down_sw, cos_zenith_angle, SW_TILE, LW_DOWN,SW_DOWN, &
       sin_theta_latitude,LS_RAIN, LS_SNOW,                 & 
       CON_RAIN, CONV_SNOW, TL_1, QW_1, PSTAR,              &
       L_TILE_PTS,                                          &
       VSHR_LAND,                                           &        
       latitude,longitude,day,time_sec,                     &
       l_cable,                                             &

       ! IN values defining field dimensions and subset to be processed :         
       OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,                  &
       LAND_PTS,LAND_PTS_TRIF,                              &
       NTILES,SM_LEVELS,NPFT_TRIF,                          &          
                                                                       
       ! IN more parameters required from boundary-layer scheme :                      
        Z1_UV,Z1_TQ,                                        &              
                                                                        
       ! IN soil/vegetation/land surface data :                                  
       LAND_INDEX,LAND_MASK,                                & 
       CANOPY_TILE,CATCH,CATCH_SNOW,HCON,HO2R2_OROG,        &            
       FLAND,                                               & 
       !      SNOW_TILE,SIL_OROG_LAND,SMVCCL,SMVCST,SMVCWT,BEXP,   &
       SNOW_TILE,SMVCCL,SMVCST,SMVCWT,BEXP,                 &            
       SATHH,SATCON,HCAP,                                   &
       SMCL,                                                &
       STHF,STHU,Z0_TILE,                                   &            

       ! IN everything not covered so far :                                       
       TIMESTEP,ZH,                                         &            
       CO2_MMR,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE, &
       FRAC,CANHT_FT,LAI_FT,                                &       
       TSOIL,TSTAR,TSTAR_LAND,TSTAR_TILE,                   &         
                                                                        
       ! OUT Diagnostic not requiring STASH flags :                              
       CD,CD_TILE,CH,CH_TILE,FQW_1,FQW_TILE,                &
       FTL_1,FTL_TILE,LE_TILE,RADNET_TILE,                  &  
       U_S_STD_TILE,                                        &
                                                                     
       ! OUT data required elsewhere in UM system :                               
       U_S,T1_SD,Q1_SD,                                     &     
       !                                                                      
       ! OUT data required elsewhere in boundary layer or surface code            
       FRACA,RHOSTAR,RESFS,RESFT,                           & 
       !H_BLEND_OROG,Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z0M_EFF,      &    
       Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z0M_EFF,                   &  
       CANHC_TILE,WT_EXT_TILE,FLAKE,                        &  
       TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC,                  &    
       SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,       &
       SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE, &
       TSOIL_TILE,T_SURF_TILE,HCONS,                        &
       SOIL_TYPE,VEG_TYPE,albsoil,                          &
       ISNOW_FLG3L,ntau,                                    &
       FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,             &
       TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,         &
       USTAR_CAB,SURF_HTF_CAB,F_ROOT                        &
       ,U_S_CAB,CH_CAB,CD_CAB                               &
       ,RECIP_L_MO_TILE,EPOT_TILE                           &
       ,SNAGE_TILE,RTSOIL_TILE                              &
       ,GFLUX_TILE,SGFLUX_TILE,VFRAC_TILE                   &
      
