  �7  �   k820309    9          19.0        %�z^                                                                                                          
       cable_cru_TRENDY.F90 CABLE_CRU              RAIN LWDN SWDN PRES QAIR TMAX TMIN UWIND VWIND PREVTMAX NEXTTMIN ERRSTATUS SECDAY PREF                      @                              
                                                           
       HANDLE_ERR GET_UNIT                      @                              
       LOGN LAND_X LAND_Y EXISTS                   @                               '8                    #WIND    #LWDOWN    #CO2AIR    #PSURF    #SNOWF 	   #AVPRECIP 
   #LAI    #LAI_T    #LAI_M    #LAI_P    #PARAMETERS    #INITIAL    #PATCH    #LAIPATCH                 �                                                               �                                                              �                                                              �                                                              �                               	                               �                               
                               �                                                              �                                                              �                                            	                   �                                    $       
                   �                                    (                          �                                    ,                          �                                    0                          �                                    4             #         @                                                      #STATUS    #MSG                                                                     
                                                   1 #         @                                                      #IUNIT                                                                                                                              @                                                                   &                                                     @                                                                   &                                                                                             8       #INPUT_DETAILS_TYPE    %         @                                                 	          #PATH    #MODE    #NCID    #CHUNKSIZE     #CACHE_SIZE !   #CACHE_NELEMS "   #CACHE_PREEMPTION #   #COMM $   #INFO %             
                                                    1           
                                                                                                              
                                                        
                                 !                     
                                 "                     
                                 #     	                
                                 $                     
                                 %                                                        &                                                       0%         @                                '                           #NCID (   #NAME )   #DIMID *             
                                  (                     
                                )                    1                                            *            %         @                                +                           #NCID ,   #DIMID -   #NAME .   #LEN /             
                                  ,                     
                                  -                                                    .                     1                                            /            %         @                                0                           #NCID 1   #NAME 2   #VARID 3             
                                  1                     
                                2                    1                                            3            %         @                                4                           #NCID 5             
                                  5                             @              @           6     'H                    #METVALS 7              �                               7                              	            &                                                             @               �           8     '�                    #MLAND 9   #NMET :   #XDIMSIZE ;   #YDIMSIZE <   #TDIMSIZE =   #CYEAR >   #METSTART ?   #METEND @   #CTSTEP A   #DTSECS B   #KTAU C   #F_ID D   #V_ID E   #AVG_LWDN F   #CO2VALS G   #DIRECTREAD H   #LEAPYEARS I   #LANDMASK J   #RUN K   #CO2 L   #NDEP M   #FORCING N   #BASEPATH O   #METPATH P   #LANDMASKFILE Q   #VAR_NAME R   #METFILE S   #MET T   #NDEPVALS U   #NDEPF_ID V   #NDEPV_ID W   #NDEP_CTSTEP X                �                               9                                �                               :                               �                               ;                               �                               <                               �                               =                               �                               >                               �                               ?                               �                               @                               �                               A             	                   �                               B     $       
                   �                               C     (                          �                               D     	       ,                   p          p 	           p 	                                      �                               E     	       P                   p          p 	           p 	                                    �                               F            x                 	            &                                                      �                               G            �                 	            &                                                        �                               H                              �                               I                            �                               J                                        &                   &                                                        �                              K            p                          �                              L            �                          �                              M            �                          �                              N            �                          �                              O     �       �                          �                              P     �       �                          �                              Q     �       K             .             �                              R     	                         p          p 	           p 	                                 .             �                              S     	       !      �            p          p 	           p 	                                              �                               T            0      H             #CRU_MET_TYPE 6   p          p            p                                     �                               U            H                	            &                                                        �                               V     �                         �                               W     �                         �                               X     �                                                         Y     �      #CRU_TYPE 8   #         @                                   Z                    #CRU [             D                                 [     �              #CRU_TYPE 8   #         @                                  \                    #CRU ]   #CYEAR ^   #PAR _   #FN `             
                                 ]     �              #CRU_TYPE 8             
                                  ^                     
                                  _                     D @                              `     �                       #         @                                  a                    #CRU b   #CO2AIR c             D                                 b     �              #CRU_TYPE 8             D                                 c     	       #         @                                  d                    #CRU e             
D @                               e     �              #CRU_TYPE 8   #         @                                  f                    #CRU g             
D @                               g     �              #CRU_TYPE 8   #         @                                  h                    #CRU i   #LASTDAYOFYEAR j   #LASTYEAROFMET k             D @                               i     �              #CRU_TYPE 8             
                                  j                     
                                  k           #         @                                   l                   #CRU_GET_SUBDIURNAL_MET%MET_TYPE m   #CRU �   #MET �   #CURYEAR �   #KTAU �   #KEND �   #LASTYEAROFMET �                     @                           m     'H                   #YEAR n   #MOY o   #CA p   #DOY q   #HOD r   #OFSD s   #FLD t   #PRECIP u   #PRECIP_SN v   #TK w   #TVAIR x   #TVRAD y   #PMB z   #UA {   #QV |   #QVAIR }   #DA ~   #DVA    #COSZEN �   #NDEP �   #PDEP �   #FSD �               �                              n                                          &                                                       �                              o            H                             &                                                       �                              p            �                 	            &                                                       �                              q            �                 	            &                                                       �                              r                             	            &                                                       �                              s            h                	            &                                                       �                              t            �                	            &                                                       �                              u            �                	            &                                                       �                              v            @             	   	            &                                                       �                              w            �             
   	            &                                                       �                              x            �                	            &                                                       �                              y                            	            &                                                       �                              z            `                	            &                                                       �                              {            �                	            &                                                       �                              |            �                	            &                                                       �                              }            8                	            &                                                       �                              ~            �                	            &                                                       �                                          �                	            &                                                       �                              �                            	            &                                                       �                              �            X                	            &                                                       �                              �            �                	            &                                                       �                              �            �                	            &                   &                                                     D @                               �     �              #CRU_TYPE 8             D                                 �     H              #CRU_GET_SUBDIURNAL_MET%MET_TYPE m             
  @                               �                     
                                  �                     
                                  �                     
  @                               �              �   '      fn#fn    �   g   b   uapp(CABLE_CRU    .  @   J   NETCDF $   n  T   J  CABLE_COMMON_MODULE %   �  Z   J  CABLE_IO_VARS_MODULE 8     �       INPUT_DETAILS_TYPE+CABLE_IO_VARS_MODULE =     H   a   INPUT_DETAILS_TYPE%WIND+CABLE_IO_VARS_MODULE ?   Z  H   a   INPUT_DETAILS_TYPE%LWDOWN+CABLE_IO_VARS_MODULE ?   �  H   a   INPUT_DETAILS_TYPE%CO2AIR+CABLE_IO_VARS_MODULE >   �  H   a   INPUT_DETAILS_TYPE%PSURF+CABLE_IO_VARS_MODULE >   2  H   a   INPUT_DETAILS_TYPE%SNOWF+CABLE_IO_VARS_MODULE A   z  H   a   INPUT_DETAILS_TYPE%AVPRECIP+CABLE_IO_VARS_MODULE <   �  H   a   INPUT_DETAILS_TYPE%LAI+CABLE_IO_VARS_MODULE >   
  H   a   INPUT_DETAILS_TYPE%LAI_T+CABLE_IO_VARS_MODULE >   R  H   a   INPUT_DETAILS_TYPE%LAI_M+CABLE_IO_VARS_MODULE >   �  H   a   INPUT_DETAILS_TYPE%LAI_P+CABLE_IO_VARS_MODULE C   �  H   a   INPUT_DETAILS_TYPE%PARAMETERS+CABLE_IO_VARS_MODULE @   *  H   a   INPUT_DETAILS_TYPE%INITIAL+CABLE_IO_VARS_MODULE >   r  H   a   INPUT_DETAILS_TYPE%PATCH+CABLE_IO_VARS_MODULE A   �  H   a   INPUT_DETAILS_TYPE%LAIPATCH+CABLE_IO_VARS_MODULE /     ]       HANDLE_ERR+CABLE_COMMON_MODULE 6   _  @   a   HANDLE_ERR%STATUS+CABLE_COMMON_MODULE 3   �  L   a   HANDLE_ERR%MSG+CABLE_COMMON_MODULE -   �  S       GET_UNIT+CABLE_COMMON_MODULE 3   >  @   a   GET_UNIT%IUNIT+CABLE_COMMON_MODULE *   ~  @       LOGN+CABLE_IO_VARS_MODULE ,   �  �       LAND_X+CABLE_IO_VARS_MODULE ,   J	  �       LAND_Y+CABLE_IO_VARS_MODULE ,   �	  X       EXISTS+CABLE_IO_VARS_MODULE !   .
  �       NF90_OPEN+NETCDF &   �
  L   a   NF90_OPEN%PATH+NETCDF &   C  @   a   NF90_OPEN%MODE+NETCDF &   �  @   a   NF90_OPEN%NCID+NETCDF +   �  @   a   NF90_OPEN%CHUNKSIZE+NETCDF ,     @   a   NF90_OPEN%CACHE_SIZE+NETCDF .   C  @   a   NF90_OPEN%CACHE_NELEMS+NETCDF 2   �  @   a   NF90_OPEN%CACHE_PREEMPTION+NETCDF &   �  @   a   NF90_OPEN%COMM+NETCDF &     @   a   NF90_OPEN%INFO+NETCDF $   C  q       NF90_NOWRITE+NETCDF &   �  o       NF90_INQ_DIMID+NETCDF +   #  @   a   NF90_INQ_DIMID%NCID+NETCDF +   c  L   a   NF90_INQ_DIMID%NAME+NETCDF ,   �  @   a   NF90_INQ_DIMID%DIMID+NETCDF .   �  x       NF90_INQUIRE_DIMENSION+NETCDF 3   g  @   a   NF90_INQUIRE_DIMENSION%NCID+NETCDF 4   �  @   a   NF90_INQUIRE_DIMENSION%DIMID+NETCDF 3   �  L   a   NF90_INQUIRE_DIMENSION%NAME+NETCDF 2   3  @   a   NF90_INQUIRE_DIMENSION%LEN+NETCDF &   s  o       NF90_INQ_VARID+NETCDF +   �  @   a   NF90_INQ_VARID%NCID+NETCDF +   "  L   a   NF90_INQ_VARID%NAME+NETCDF ,   n  @   a   NF90_INQ_VARID%VARID+NETCDF "   �  Z       NF90_CLOSE+NETCDF '     @   a   NF90_CLOSE%NCID+NETCDF    H  ]       CRU_MET_TYPE %   �  �   a   CRU_MET_TYPE%METVALS    9  �      CRU_TYPE       H   a   CRU_TYPE%MLAND    h  H   a   CRU_TYPE%NMET "   �  H   a   CRU_TYPE%XDIMSIZE "   �  H   a   CRU_TYPE%YDIMSIZE "   @  H   a   CRU_TYPE%TDIMSIZE    �  H   a   CRU_TYPE%CYEAR "   �  H   a   CRU_TYPE%METSTART       H   a   CRU_TYPE%METEND     `  H   a   CRU_TYPE%CTSTEP     �  H   a   CRU_TYPE%DTSECS    �  H   a   CRU_TYPE%KTAU    8  �   a   CRU_TYPE%F_ID    �  �   a   CRU_TYPE%V_ID "   p  �   a   CRU_TYPE%AVG_LWDN !     �   a   CRU_TYPE%CO2VALS $   �  H   a   CRU_TYPE%DIRECTREAD #   �  H   a   CRU_TYPE%LEAPYEARS "   (  �   a   CRU_TYPE%LANDMASK    �  P   a   CRU_TYPE%RUN    $  P   a   CRU_TYPE%CO2    t  P   a   CRU_TYPE%NDEP !   �  P   a   CRU_TYPE%FORCING "     P   a   CRU_TYPE%BASEPATH !   d  P   a   CRU_TYPE%METPATH &   �  P   a   CRU_TYPE%LANDMASKFILE "     �   a   CRU_TYPE%VAR_NAME !   �  �   a   CRU_TYPE%METFILE    L  �   a   CRU_TYPE%MET "   �  �   a   CRU_TYPE%NDEPVALS "   �   H   a   CRU_TYPE%NDEPF_ID "   �   H   a   CRU_TYPE%NDEPV_ID %   !  H   a   CRU_TYPE%NDEP_CTSTEP    f!  N       CRU    �!  Q       CRU_INIT    "  V   a   CRU_INIT%CRU !   ["  m       CRU_GET_FILENAME %   �"  V   a   CRU_GET_FILENAME%CRU '   #  @   a   CRU_GET_FILENAME%CYEAR %   ^#  @   a   CRU_GET_FILENAME%PAR $   �#  P   a   CRU_GET_FILENAME%FN    �#  ]       GET_CRU_CO2     K$  V   a   GET_CRU_CO2%CRU #   �$  @   a   GET_CRU_CO2%CO2AIR    �$  Q       GET_CRU_NDEP !   2%  V   a   GET_CRU_NDEP%CRU    �%  Q       OPEN_CRU_MET !   �%  V   a   OPEN_CRU_MET%CRU "   /&  w       CRU_GET_DAILY_MET &   �&  V   a   CRU_GET_DAILY_MET%CRU 0   �&  @   a   CRU_GET_DAILY_MET%LASTDAYOFYEAR 0   <'  @   a   CRU_GET_DAILY_MET%LASTYEAROFMET '   |'  �       CRU_GET_SUBDIURNAL_MET D   /(  '     CRU_GET_SUBDIURNAL_MET%MET_TYPE+CABLE_DEF_TYPES_MOD I   V)  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%YEAR+CABLE_DEF_TYPES_MOD H   �)  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%MOY+CABLE_DEF_TYPES_MOD G   ~*  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%CA+CABLE_DEF_TYPES_MOD H   +  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%DOY+CABLE_DEF_TYPES_MOD H   �+  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%HOD+CABLE_DEF_TYPES_MOD I   :,  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%OFSD+CABLE_DEF_TYPES_MOD H   �,  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%FLD+CABLE_DEF_TYPES_MOD K   b-  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%PRECIP+CABLE_DEF_TYPES_MOD N   �-  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%PRECIP_SN+CABLE_DEF_TYPES_MOD G   �.  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%TK+CABLE_DEF_TYPES_MOD J   /  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%TVAIR+CABLE_DEF_TYPES_MOD J   �/  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%TVRAD+CABLE_DEF_TYPES_MOD H   F0  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%PMB+CABLE_DEF_TYPES_MOD G   �0  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%UA+CABLE_DEF_TYPES_MOD G   n1  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%QV+CABLE_DEF_TYPES_MOD J   2  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%QVAIR+CABLE_DEF_TYPES_MOD G   �2  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%DA+CABLE_DEF_TYPES_MOD H   *3  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%DVA+CABLE_DEF_TYPES_MOD K   �3  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%COSZEN+CABLE_DEF_TYPES_MOD I   R4  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%NDEP+CABLE_DEF_TYPES_MOD I   �4  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%PDEP+CABLE_DEF_TYPES_MOD H   z5  �   a   CRU_GET_SUBDIURNAL_MET%MET_TYPE%FSD+CABLE_DEF_TYPES_MOD +   &6  V   a   CRU_GET_SUBDIURNAL_MET%CRU +   |6  m   a   CRU_GET_SUBDIURNAL_MET%MET /   �6  @   a   CRU_GET_SUBDIURNAL_MET%CURYEAR ,   )7  @   a   CRU_GET_SUBDIURNAL_MET%KTAU ,   i7  @   a   CRU_GET_SUBDIURNAL_MET%KEND 5   �7  @   a   CRU_GET_SUBDIURNAL_MET%LASTYEAROFMET 