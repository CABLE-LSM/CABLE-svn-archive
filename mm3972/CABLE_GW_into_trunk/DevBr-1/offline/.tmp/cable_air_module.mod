  Èe  á   k820309    9          19.0        z=c                                                                                                          
       cable_air.F90 CABLE_AIR_MODULE              DEFINE_AIR                                                     
       IAIR_TYPE POINT2CONSTANTS               @                                   
   u #DRIVER_TYPE_PTR    #CBM_TYPE_PTR    #AIR_TYPE_PTR    #ALBEDO_TYPE_PTR    #CANOPY_TYPE_PTR    #CARBON_TYPE_PTR    #RAD_TYPE_PTR    #ROUGH_TYPE_PTR    #SSNOW_TYPE_PTR    #GWHYDRO_TYPE_PTR    #         @     @                                                #C                                                                    #DRIVER_TYPE    #         @     @                                                #C                                                    (               #ICBM_TYPE    #         @     @                                               #C 	                                              	     `               #IAIR_TYPE 
   #         @     @                                                #C                                                                   #IALBEDO_TYPE    #         @     @                                                #C                                                    (              #ICANOPY_TYPE    #         @     @                                                #C                                                                   #ICARBON_TYPE    #         @     @                                                #C                                                                   #IRAD_TYPE    #         @     @                                                #C                                                    `               #IROUGH_TYPE    #         @     @                                                #C                                                    `               #ISSNOW_TYPE    #         @     @                                                #C                                                    P               #IGWHYDRO_TYPE                   @  @                           
     '`                    #TFRZ     #RMAIR !   #RGAS "   #TETENA #   #TETENB $   #TETENC %   #TETENA_ICE &   #TETENB_ICE '   #TETENC_ICE (   #CAPP )   #RMH2O *   #HL +                                                               	                                              !               	                                              "               	                                              #               	                                              $                	                                              %     (          	                                              &     0          	                                              '     8          	                                              (     @       	   	                                              )     H       
   	                                              *     P          	                                              +     X          	                  @  @                               '                     #TFRZ ,   #EMSOIL -   #EMLEAF .   #SBOLTZ /                                              ,                	                                              -               	                                              .               	                                              /               	                  @  @                               '(                    #GRAV 0   #CAPP 1   #EMLEAF 2   #EMSOIL 3   #SBOLTZ 4                                              0                	                                              1               	                                              2               	                                              3               	                                              4                	                  @  @                               '                    #TFRZ 5   #LAI_THRESH 6   #RAD_THRESH 7                                              5                	                                              6               	                                              7               	                  @  @                               '(             %      #TFRZ 8   #RMAIR 9   #RGAS :   #DHEAT ;   #ZETNEG <   #ZETMUL =   #ZETPOS >   #GRAV ?   #UMIN @   #TETENA A   #TETENB B   #TETENC C   #RHOW D   #CTL E   #CSW F   #EMLEAF G   #EMSOIL H   #SBOLTZ I   #PRANDT J   #CAPP K   #RMH2O L   #APOL M   #A33 N   #VONK O   #ZETA0 P   #RGSWC Q   #GAM0 R   #GAM1 S   #GAM2 T   #RGBWC U   #TREFK V   #PI_C W   #LAI_THRESH X   #TETENA_ICE Y   #TETENB_ICE Z   #TETENC_ICE [   #MAXITER \                                              8                	                                              9               	                                              :               	                                              ;               	                                              <                	                                              =     (          	                                              >     0          	                                              ?     8          	                                              @     @       	   	                                              A     H       
   	                                              B     P          	                                              C     X          	                                              D     `          	                                              E     h          	                                              F     p          	                                              G     x          	                                              H               	                                              I               	                                              J               	                                              K               	                                              L                	                                              M     ¨          	                                              N     °          	                                              O     ¸          	                                              P     À          	                                              Q     È          	                                              R     Ð          	                                              S     Ø          	                                              T     à          	                                              U     è          	                                              V     ð          	                                              W     ø           	                                              X            !   	                                              Y           "   	                                              Z           #   	                                              [           $   	                                              \            %                     @  @                               '                    #TFRZ ]                                              ]                	                  @  @                               '              
      #TFRZ ^   #EMSOIL _   #EMLEAF `   #SBOLTZ a   #CAPP b   #LAI_THRESH c   #RAD_THRESH d   #PI180 e   #PI_C f   #GAUSS_W g                                              ^                	                                              _               	                                              `               	                                              a               	                                              b                	                                              c     (          	                                              d     0          	                                              e     8          	                                              f     @       	   	                                             g            H              
   	            &                                                          @  @                               '`                    #CSD h   #CRD i   #CCD j   #CCW_C k   #USUHM l   #VONK m   #A33 n   #CTL o   #ZDLIN p   #CSW q   #GRAV r   #LAI_THRESH s                                              h                	                                              i               	                                              j               	                                              k               	                                              l                	                                              m     (          	                                              n     0          	                                              o     8          	                                              p     @       	   	                                              q     H       
   	                                              r     P          	                                              s     X          	                  @  @                               '`                    #CAPP t   #TFRZ u   #HL v   #HLF w   #HLS x   #DENSITY_LIQ y   #DENSITY_ICE z   #CGSNOW {   #CSWAT |   #CSICE }   #CS_RHO_WAT ~   #CS_RHO_ICE                                               t                	                                              u               	                                              v               	                                              w               	                                              x                	                                              y     (          	                                              z     0          	                                              {     8          	                                              |     @       	   	                                              }     H       
   	                                              ~     P          	                                                   X          	                  @  @                               'P              
      #TFRZ    #HL    #HLF    #HLS    #DENSITY_LIQ    #DENSITY_ICE    #CGSNOW    #CS_RHO_WAT    #CS_RHO_ICE    #PI                                                               	                                                             	                                                             	                                                             	                                                              	                                                   (          	                                                   0          	                                                   8          	                                                   @       	   	                                                   H       
   	                  @  @                              '¼              /      #CAPP    #HL    #HLF    #HLS    #DHEAT    #GRAV    #RGAS    #RMAIR    #RMH2O    #SBOLTZ    #TFRZ    #CGSNOW    #CS_RHO_ICE    #CS_RHO_WAT    #CSICE    #CSWAT    #DENSITY_LIQ    #DENSITY_ICE    #TETENA    #TETENB    #TETENC    #TETENA_ICE     #TETENB_ICE ¡   #TETENC_ICE ¢   #VONK £   #A33 ¤   #CSW ¥   #CTL ¦   #APOL §   #PRANDT ¨   #SCHMID ©   #DIFFWC ª   #RHOW «   #EMLEAF ¬   #EMSOIL ­   #CRD ®   #CSD ¯   #BETA2 °   #CCD ±   #CCW_C ²   #USUHM ³   #ZETMUL ´   #ZETA0 µ   #ZETNEG ¶   #ZETPOS ·   #ZDLIN ¸   #UMIN ¹                                                             	                           Ý                     	                 ö({D            1004.64                                                             	                           Ý                     	                 `¬J            2.5014E6                                                             	                           Ý                     	                  £H            0.334E6                                                             	                           Ý                     	                 à-J            2.8350E6                                                             	                           Ý                     	                 æZ´7            21.5E-6                                                             	                           Ý                     	                 ðA            9.8086                                                             	                           Ý                     	                 _A            8.3143                                                             	                           Ý                     	                 ~Rí<            0.02897                                                            	  	                           Ý                     	                 J<            0.018016                                                    $       
  	                           Ý                     	                 Os3            5.67E-8                                                    (         	                           Ý                     	                 {C            273.16                                                    ,         	                           Ý                     	                   E            2090.0                                                    0         	                           Ý                     	                  ìI            1.9341E6                                                    4         	                           Ý                     	                  ¹J            4.218E6                                                    8         	                           Ý                     	                  @E            2.100E3                                                    <         	                           Ý                     	                  ÐE            4.218E3                                                    @         	                           Ý                     	                   zD            1000.0                                                    D         	                           Ý                     	                  @fD            921.0                                                    H         	                           Ý                     	                 ZdÃ@            6.106                                                    L         	                           Ý                     	                 ö(A            17.27                                                    P         	                           Ý                     	                 ÍLmC            237.3                                                     T         	                           Ý                     	                 sÃ@            6.1078                                               ¡     X         	                           Ý                     	                   ¯A            21.875                                               ¢     \         	                           Ý                     	                  ÀC            265.5                                               £     `         	                           Ý                     	                 ÍÌÌ>            0.40                                               ¤     d         	                           Ý                     	                    ?            1.25                                               ¥     h         	                           Ý                     	                    ?            0.50                                               ¦     l         	                           Ý                     	                 ÍÌÌ>            0.40                                               §     p         	                           Ý                     	                 333?            0.70                                               ¨     t         	                           Ý                     	                 Â5?            0.71                                               ©     x         	                           Ý                     	                 ?            0.60                                               ª     |          	                           Ý                     	                 ÍÌÌ?            1.60                                               «            !  	                           Ý                     	                   zD            1000.0                                               ¬            "  	                           Ý                     	                   ?            1.0                                               ­            #  	                           Ý                     	                   ?            1.0                                               ®            $  	                           Ý                     	                 >            0.3                                               ¯            %  	                           Ý                     	                 ¦D;            0.003                                               °            &  	                           Ý                       	                    ÈB                                                           ±            '  	                           Ý                     	                   pA            15.0                                               ²            (  	                           Ý                     	                    @            2.0                                               ³             )  	                           Ý                     	                 >            0.3                                               ´     ¤       *  	                           Ý                     	                 ÍÌÌ>            0.4                                               µ     ¨       +  	                           Ý                     	                                 0.0                                               ¶     ¬       ,  	                           Ý                       	                    pÁ                                                           ·     °       -  	                           Ý                     	                   ?            1.0                                               ¸     ´       .  	                           Ý                     	                   ?            1.0                                               ¹     ¸       /  	                           Ý                     	                 ÍÌÌ=            0.1    #         @                                   º                  #DEFINE_AIR%MP »   #DEFINE_AIR%AIR_TYPE ¼   #DEFINE_AIR%MET_TYPE Æ   #MET Ý   #AIR Þ                     @                           ¼     '             	      #RHO ½   #VOLM ¾   #RLAM ¿   #QSAT À   #EPSI Á   #VISC Â   #PSYC Ã   #DSATDK Ä   #CMOLAR Å                                             ½                              	            &                                                                                     ¾            H                 	            &                                                                                     ¿                             	            &                                                                                     À            Ø                 	            &                                                                                     Á                             	            &                                                                                     Â            h                	            &                                                                                     Ã            °                	            &                                                                                     Ä            ø                	            &                                                                                     Å            @             	   	            &                                                             @                           Æ     'H                   #YEAR Ç   #MOY È   #CA É   #DOY Ê   #HOD Ë   #OFSD Ì   #FLD Í   #PRECIP Î   #PRECIP_SN Ï   #TK Ð   #TVAIR Ñ   #TVRAD Ò   #PMB Ó   #UA Ô   #QV Õ   #QVAIR Ö   #DA ×   #DVA Ø   #COSZEN Ù   #NDEP Ú   #PDEP Û   #FSD Ü                                             Ç                                          &                                                                                     È            H                             &                                                                                     É                             	            &                                                                                     Ê            Ø                 	            &                                                                                     Ë                             	            &                                                                                     Ì            h                	            &                                                                                     Í            °                	            &                                                                                     Î            ø                	            &                                                                                     Ï            @             	   	            &                                                                                     Ð                         
   	            &                                                                                     Ñ            Ð                	            &                                                                                     Ò                            	            &                                                                                     Ó            `                	            &                                                                                     Ô            ¨                	            &                                                                                     Õ            ð                	            &                                                                                     Ö            8                	            &                                                                                     ×                            	            &                                                                                     Ø            È                	            &                                                                                     Ù                            	            &                                                                                     Ú            X                	            &                                                                                     Û                             	            &                                                                                     Ü            è                	            &                   &                                                                                      »                      
                                  Ý     H             #DEFINE_AIR%MET_TYPE Æ             
D                                 Þ                   #DEFINE_AIR%AIR_TYPE ¼          '      fn#fn &   Ç      b   uapp(CABLE_AIR_MODULE "   â   Z   J  CABLE_DATA_MODULE 6   <        gen@POINT2CONSTANTS+CABLE_DATA_MODULE 2   D  O      DRIVER_TYPE_PTR+CABLE_DATA_MODULE 4     Y   a   DRIVER_TYPE_PTR%C+CABLE_DATA_MODULE /   ì  O      CBM_TYPE_PTR+CABLE_DATA_MODULE 1   ;  W   a   CBM_TYPE_PTR%C+CABLE_DATA_MODULE /     O      AIR_TYPE_PTR+CABLE_DATA_MODULE 1   á  W   a   AIR_TYPE_PTR%C+CABLE_DATA_MODULE 2   8  O      ALBEDO_TYPE_PTR+CABLE_DATA_MODULE 4     Z   a   ALBEDO_TYPE_PTR%C+CABLE_DATA_MODULE 2   á  O      CANOPY_TYPE_PTR+CABLE_DATA_MODULE 4   0  Z   a   CANOPY_TYPE_PTR%C+CABLE_DATA_MODULE 2     O      CARBON_TYPE_PTR+CABLE_DATA_MODULE 4   Ù  Z   a   CARBON_TYPE_PTR%C+CABLE_DATA_MODULE /   3  O      RAD_TYPE_PTR+CABLE_DATA_MODULE 1     W   a   RAD_TYPE_PTR%C+CABLE_DATA_MODULE 1   Ù  O      ROUGH_TYPE_PTR+CABLE_DATA_MODULE 3   (  Y   a   ROUGH_TYPE_PTR%C+CABLE_DATA_MODULE 1     O      SSNOW_TYPE_PTR+CABLE_DATA_MODULE 3   Ð  Y   a   SSNOW_TYPE_PTR%C+CABLE_DATA_MODULE 3   )  O      GWHYDRO_TYPE_PTR+CABLE_DATA_MODULE 5   x  [   a   GWHYDRO_TYPE_PTR%C+CABLE_DATA_MODULE ,   Ó  à      IAIR_TYPE+CABLE_DATA_MODULE 1   ³	  H   a   IAIR_TYPE%TFRZ+CABLE_DATA_MODULE 2   û	  H   a   IAIR_TYPE%RMAIR+CABLE_DATA_MODULE 1   C
  H   a   IAIR_TYPE%RGAS+CABLE_DATA_MODULE 3   
  H   a   IAIR_TYPE%TETENA+CABLE_DATA_MODULE 3   Ó
  H   a   IAIR_TYPE%TETENB+CABLE_DATA_MODULE 3     H   a   IAIR_TYPE%TETENC+CABLE_DATA_MODULE 7   c  H   a   IAIR_TYPE%TETENA_ICE+CABLE_DATA_MODULE 7   «  H   a   IAIR_TYPE%TETENB_ICE+CABLE_DATA_MODULE 7   ó  H   a   IAIR_TYPE%TETENC_ICE+CABLE_DATA_MODULE 1   ;  H   a   IAIR_TYPE%CAPP+CABLE_DATA_MODULE 2     H   a   IAIR_TYPE%RMH2O+CABLE_DATA_MODULE /   Ë  H   a   IAIR_TYPE%HL+CABLE_DATA_MODULE .     ~      DRIVER_TYPE+CABLE_DATA_MODULE 3     H   a   DRIVER_TYPE%TFRZ+CABLE_DATA_MODULE 5   Ù  H   a   DRIVER_TYPE%EMSOIL+CABLE_DATA_MODULE 5   !  H   a   DRIVER_TYPE%EMLEAF+CABLE_DATA_MODULE 5   i  H   a   DRIVER_TYPE%SBOLTZ+CABLE_DATA_MODULE ,   ±        ICBM_TYPE+CABLE_DATA_MODULE 1   9  H   a   ICBM_TYPE%GRAV+CABLE_DATA_MODULE 1     H   a   ICBM_TYPE%CAPP+CABLE_DATA_MODULE 3   É  H   a   ICBM_TYPE%EMLEAF+CABLE_DATA_MODULE 3     H   a   ICBM_TYPE%EMSOIL+CABLE_DATA_MODULE 3   Y  H   a   ICBM_TYPE%SBOLTZ+CABLE_DATA_MODULE /   ¡  z      IALBEDO_TYPE+CABLE_DATA_MODULE 4     H   a   IALBEDO_TYPE%TFRZ+CABLE_DATA_MODULE :   c  H   a   IALBEDO_TYPE%LAI_THRESH+CABLE_DATA_MODULE :   «  H   a   IALBEDO_TYPE%RAD_THRESH+CABLE_DATA_MODULE /   ó  õ     ICANOPY_TYPE+CABLE_DATA_MODULE 4   è  H   a   ICANOPY_TYPE%TFRZ+CABLE_DATA_MODULE 5   0  H   a   ICANOPY_TYPE%RMAIR+CABLE_DATA_MODULE 4   x  H   a   ICANOPY_TYPE%RGAS+CABLE_DATA_MODULE 5   À  H   a   ICANOPY_TYPE%DHEAT+CABLE_DATA_MODULE 6     H   a   ICANOPY_TYPE%ZETNEG+CABLE_DATA_MODULE 6   P  H   a   ICANOPY_TYPE%ZETMUL+CABLE_DATA_MODULE 6     H   a   ICANOPY_TYPE%ZETPOS+CABLE_DATA_MODULE 4   à  H   a   ICANOPY_TYPE%GRAV+CABLE_DATA_MODULE 4   (  H   a   ICANOPY_TYPE%UMIN+CABLE_DATA_MODULE 6   p  H   a   ICANOPY_TYPE%TETENA+CABLE_DATA_MODULE 6   ¸  H   a   ICANOPY_TYPE%TETENB+CABLE_DATA_MODULE 6      H   a   ICANOPY_TYPE%TETENC+CABLE_DATA_MODULE 4   H  H   a   ICANOPY_TYPE%RHOW+CABLE_DATA_MODULE 3     H   a   ICANOPY_TYPE%CTL+CABLE_DATA_MODULE 3   Ø  H   a   ICANOPY_TYPE%CSW+CABLE_DATA_MODULE 6      H   a   ICANOPY_TYPE%EMLEAF+CABLE_DATA_MODULE 6   h  H   a   ICANOPY_TYPE%EMSOIL+CABLE_DATA_MODULE 6   °  H   a   ICANOPY_TYPE%SBOLTZ+CABLE_DATA_MODULE 6   ø  H   a   ICANOPY_TYPE%PRANDT+CABLE_DATA_MODULE 4   @  H   a   ICANOPY_TYPE%CAPP+CABLE_DATA_MODULE 5     H   a   ICANOPY_TYPE%RMH2O+CABLE_DATA_MODULE 4   Ð  H   a   ICANOPY_TYPE%APOL+CABLE_DATA_MODULE 3     H   a   ICANOPY_TYPE%A33+CABLE_DATA_MODULE 4   `  H   a   ICANOPY_TYPE%VONK+CABLE_DATA_MODULE 5   ¨  H   a   ICANOPY_TYPE%ZETA0+CABLE_DATA_MODULE 5   ð  H   a   ICANOPY_TYPE%RGSWC+CABLE_DATA_MODULE 4   8  H   a   ICANOPY_TYPE%GAM0+CABLE_DATA_MODULE 4     H   a   ICANOPY_TYPE%GAM1+CABLE_DATA_MODULE 4   È  H   a   ICANOPY_TYPE%GAM2+CABLE_DATA_MODULE 5     H   a   ICANOPY_TYPE%RGBWC+CABLE_DATA_MODULE 5   X  H   a   ICANOPY_TYPE%TREFK+CABLE_DATA_MODULE 4      H   a   ICANOPY_TYPE%PI_C+CABLE_DATA_MODULE :   è  H   a   ICANOPY_TYPE%LAI_THRESH+CABLE_DATA_MODULE :   0  H   a   ICANOPY_TYPE%TETENA_ICE+CABLE_DATA_MODULE :   x  H   a   ICANOPY_TYPE%TETENB_ICE+CABLE_DATA_MODULE :   À  H   a   ICANOPY_TYPE%TETENC_ICE+CABLE_DATA_MODULE 7     H   a   ICANOPY_TYPE%MAXITER+CABLE_DATA_MODULE /   P  Z      ICARBON_TYPE+CABLE_DATA_MODULE 4   ª  H   a   ICARBON_TYPE%TFRZ+CABLE_DATA_MODULE ,   ò  Ê      IRAD_TYPE+CABLE_DATA_MODULE 1   ¼  H   a   IRAD_TYPE%TFRZ+CABLE_DATA_MODULE 3      H   a   IRAD_TYPE%EMSOIL+CABLE_DATA_MODULE 3   L   H   a   IRAD_TYPE%EMLEAF+CABLE_DATA_MODULE 3      H   a   IRAD_TYPE%SBOLTZ+CABLE_DATA_MODULE 1   Ü   H   a   IRAD_TYPE%CAPP+CABLE_DATA_MODULE 7   $!  H   a   IRAD_TYPE%LAI_THRESH+CABLE_DATA_MODULE 7   l!  H   a   IRAD_TYPE%RAD_THRESH+CABLE_DATA_MODULE 2   ´!  H   a   IRAD_TYPE%PI180+CABLE_DATA_MODULE 1   ü!  H   a   IRAD_TYPE%PI_C+CABLE_DATA_MODULE 4   D"     a   IRAD_TYPE%GAUSS_W+CABLE_DATA_MODULE .   Ø"  Ë      IROUGH_TYPE+CABLE_DATA_MODULE 2   £#  H   a   IROUGH_TYPE%CSD+CABLE_DATA_MODULE 2   ë#  H   a   IROUGH_TYPE%CRD+CABLE_DATA_MODULE 2   3$  H   a   IROUGH_TYPE%CCD+CABLE_DATA_MODULE 4   {$  H   a   IROUGH_TYPE%CCW_C+CABLE_DATA_MODULE 4   Ã$  H   a   IROUGH_TYPE%USUHM+CABLE_DATA_MODULE 3   %  H   a   IROUGH_TYPE%VONK+CABLE_DATA_MODULE 2   S%  H   a   IROUGH_TYPE%A33+CABLE_DATA_MODULE 2   %  H   a   IROUGH_TYPE%CTL+CABLE_DATA_MODULE 4   ã%  H   a   IROUGH_TYPE%ZDLIN+CABLE_DATA_MODULE 2   +&  H   a   IROUGH_TYPE%CSW+CABLE_DATA_MODULE 3   s&  H   a   IROUGH_TYPE%GRAV+CABLE_DATA_MODULE 9   »&  H   a   IROUGH_TYPE%LAI_THRESH+CABLE_DATA_MODULE .   '  â      ISSNOW_TYPE+CABLE_DATA_MODULE 3   å'  H   a   ISSNOW_TYPE%CAPP+CABLE_DATA_MODULE 3   -(  H   a   ISSNOW_TYPE%TFRZ+CABLE_DATA_MODULE 1   u(  H   a   ISSNOW_TYPE%HL+CABLE_DATA_MODULE 2   ½(  H   a   ISSNOW_TYPE%HLF+CABLE_DATA_MODULE 2   )  H   a   ISSNOW_TYPE%HLS+CABLE_DATA_MODULE :   M)  H   a   ISSNOW_TYPE%DENSITY_LIQ+CABLE_DATA_MODULE :   )  H   a   ISSNOW_TYPE%DENSITY_ICE+CABLE_DATA_MODULE 5   Ý)  H   a   ISSNOW_TYPE%CGSNOW+CABLE_DATA_MODULE 4   %*  H   a   ISSNOW_TYPE%CSWAT+CABLE_DATA_MODULE 4   m*  H   a   ISSNOW_TYPE%CSICE+CABLE_DATA_MODULE 9   µ*  H   a   ISSNOW_TYPE%CS_RHO_WAT+CABLE_DATA_MODULE 9   ý*  H   a   ISSNOW_TYPE%CS_RHO_ICE+CABLE_DATA_MODULE 0   E+  Ê      IGWHYDRO_TYPE+CABLE_DATA_MODULE 5   ,  H   a   IGWHYDRO_TYPE%TFRZ+CABLE_DATA_MODULE 3   W,  H   a   IGWHYDRO_TYPE%HL+CABLE_DATA_MODULE 4   ,  H   a   IGWHYDRO_TYPE%HLF+CABLE_DATA_MODULE 4   ç,  H   a   IGWHYDRO_TYPE%HLS+CABLE_DATA_MODULE <   /-  H   a   IGWHYDRO_TYPE%DENSITY_LIQ+CABLE_DATA_MODULE <   w-  H   a   IGWHYDRO_TYPE%DENSITY_ICE+CABLE_DATA_MODULE 7   ¿-  H   a   IGWHYDRO_TYPE%CGSNOW+CABLE_DATA_MODULE ;   .  H   a   IGWHYDRO_TYPE%CS_RHO_WAT+CABLE_DATA_MODULE ;   O.  H   a   IGWHYDRO_TYPE%CS_RHO_ICE+CABLE_DATA_MODULE 3   .  H   a   IGWHYDRO_TYPE%PI+CABLE_DATA_MODULE 5   ß.  l     PHYSICAL_CONSTANTS+CABLE_DATA_MODULE :   K1  «   a   PHYSICAL_CONSTANTS%CAPP+CABLE_DATA_MODULE 8   ö1  ¬   a   PHYSICAL_CONSTANTS%HL+CABLE_DATA_MODULE 9   ¢2  «   a   PHYSICAL_CONSTANTS%HLF+CABLE_DATA_MODULE 9   M3  ¬   a   PHYSICAL_CONSTANTS%HLS+CABLE_DATA_MODULE ;   ù3  «   a   PHYSICAL_CONSTANTS%DHEAT+CABLE_DATA_MODULE :   ¤4  ª   a   PHYSICAL_CONSTANTS%GRAV+CABLE_DATA_MODULE :   N5  ª   a   PHYSICAL_CONSTANTS%RGAS+CABLE_DATA_MODULE ;   ø5  «   a   PHYSICAL_CONSTANTS%RMAIR+CABLE_DATA_MODULE ;   £6  ¬   a   PHYSICAL_CONSTANTS%RMH2O+CABLE_DATA_MODULE <   O7  «   a   PHYSICAL_CONSTANTS%SBOLTZ+CABLE_DATA_MODULE :   ú7  ª   a   PHYSICAL_CONSTANTS%TFRZ+CABLE_DATA_MODULE <   ¤8  ª   a   PHYSICAL_CONSTANTS%CGSNOW+CABLE_DATA_MODULE @   N9  ¬   a   PHYSICAL_CONSTANTS%CS_RHO_ICE+CABLE_DATA_MODULE @   ú9  «   a   PHYSICAL_CONSTANTS%CS_RHO_WAT+CABLE_DATA_MODULE ;   ¥:  «   a   PHYSICAL_CONSTANTS%CSICE+CABLE_DATA_MODULE ;   P;  «   a   PHYSICAL_CONSTANTS%CSWAT+CABLE_DATA_MODULE A   û;  ª   a   PHYSICAL_CONSTANTS%DENSITY_LIQ+CABLE_DATA_MODULE A   ¥<  ©   a   PHYSICAL_CONSTANTS%DENSITY_ICE+CABLE_DATA_MODULE <   N=  ©   a   PHYSICAL_CONSTANTS%TETENA+CABLE_DATA_MODULE <   ÷=  ©   a   PHYSICAL_CONSTANTS%TETENB+CABLE_DATA_MODULE <    >  ©   a   PHYSICAL_CONSTANTS%TETENC+CABLE_DATA_MODULE @   I?  ª   a   PHYSICAL_CONSTANTS%TETENA_ICE+CABLE_DATA_MODULE @   ó?  ª   a   PHYSICAL_CONSTANTS%TETENB_ICE+CABLE_DATA_MODULE @   @  ©   a   PHYSICAL_CONSTANTS%TETENC_ICE+CABLE_DATA_MODULE :   FA  ¨   a   PHYSICAL_CONSTANTS%VONK+CABLE_DATA_MODULE 9   îA  ¨   a   PHYSICAL_CONSTANTS%A33+CABLE_DATA_MODULE 9   B  ¨   a   PHYSICAL_CONSTANTS%CSW+CABLE_DATA_MODULE 9   >C  ¨   a   PHYSICAL_CONSTANTS%CTL+CABLE_DATA_MODULE :   æC  ¨   a   PHYSICAL_CONSTANTS%APOL+CABLE_DATA_MODULE <   D  ¨   a   PHYSICAL_CONSTANTS%PRANDT+CABLE_DATA_MODULE <   6E  ¨   a   PHYSICAL_CONSTANTS%SCHMID+CABLE_DATA_MODULE <   ÞE  ¨   a   PHYSICAL_CONSTANTS%DIFFWC+CABLE_DATA_MODULE :   F  ª   a   PHYSICAL_CONSTANTS%RHOW+CABLE_DATA_MODULE <   0G  §   a   PHYSICAL_CONSTANTS%EMLEAF+CABLE_DATA_MODULE <   ×G  §   a   PHYSICAL_CONSTANTS%EMSOIL+CABLE_DATA_MODULE 9   ~H  §   a   PHYSICAL_CONSTANTS%CRD+CABLE_DATA_MODULE 9   %I  ©   a   PHYSICAL_CONSTANTS%CSD+CABLE_DATA_MODULE ;   ÎI  ¤   a   PHYSICAL_CONSTANTS%BETA2+CABLE_DATA_MODULE 9   rJ  ¨   a   PHYSICAL_CONSTANTS%CCD+CABLE_DATA_MODULE ;   K  §   a   PHYSICAL_CONSTANTS%CCW_C+CABLE_DATA_MODULE ;   ÁK  §   a   PHYSICAL_CONSTANTS%USUHM+CABLE_DATA_MODULE <   hL  §   a   PHYSICAL_CONSTANTS%ZETMUL+CABLE_DATA_MODULE ;   M  §   a   PHYSICAL_CONSTANTS%ZETA0+CABLE_DATA_MODULE <   ¶M  ¤   a   PHYSICAL_CONSTANTS%ZETNEG+CABLE_DATA_MODULE <   ZN  §   a   PHYSICAL_CONSTANTS%ZETPOS+CABLE_DATA_MODULE ;   O  §   a   PHYSICAL_CONSTANTS%ZDLIN+CABLE_DATA_MODULE :   ¨O  §   a   PHYSICAL_CONSTANTS%UMIN+CABLE_DATA_MODULE    OP         DEFINE_AIR 8   îP  ­      DEFINE_AIR%AIR_TYPE+CABLE_DEF_TYPES_MOD <   Q     a   DEFINE_AIR%AIR_TYPE%RHO+CABLE_DEF_TYPES_MOD =   /R     a   DEFINE_AIR%AIR_TYPE%VOLM+CABLE_DEF_TYPES_MOD =   ÃR     a   DEFINE_AIR%AIR_TYPE%RLAM+CABLE_DEF_TYPES_MOD =   WS     a   DEFINE_AIR%AIR_TYPE%QSAT+CABLE_DEF_TYPES_MOD =   ëS     a   DEFINE_AIR%AIR_TYPE%EPSI+CABLE_DEF_TYPES_MOD =   T     a   DEFINE_AIR%AIR_TYPE%VISC+CABLE_DEF_TYPES_MOD =   U     a   DEFINE_AIR%AIR_TYPE%PSYC+CABLE_DEF_TYPES_MOD ?   §U     a   DEFINE_AIR%AIR_TYPE%DSATDK+CABLE_DEF_TYPES_MOD ?   ;V     a   DEFINE_AIR%AIR_TYPE%CMOLAR+CABLE_DEF_TYPES_MOD 8   ÏV  '     DEFINE_AIR%MET_TYPE+CABLE_DEF_TYPES_MOD =   öW     a   DEFINE_AIR%MET_TYPE%YEAR+CABLE_DEF_TYPES_MOD <   X     a   DEFINE_AIR%MET_TYPE%MOY+CABLE_DEF_TYPES_MOD ;   Y     a   DEFINE_AIR%MET_TYPE%CA+CABLE_DEF_TYPES_MOD <   ²Y     a   DEFINE_AIR%MET_TYPE%DOY+CABLE_DEF_TYPES_MOD <   FZ     a   DEFINE_AIR%MET_TYPE%HOD+CABLE_DEF_TYPES_MOD =   ÚZ     a   DEFINE_AIR%MET_TYPE%OFSD+CABLE_DEF_TYPES_MOD <   n[     a   DEFINE_AIR%MET_TYPE%FLD+CABLE_DEF_TYPES_MOD ?   \     a   DEFINE_AIR%MET_TYPE%PRECIP+CABLE_DEF_TYPES_MOD B   \     a   DEFINE_AIR%MET_TYPE%PRECIP_SN+CABLE_DEF_TYPES_MOD ;   *]     a   DEFINE_AIR%MET_TYPE%TK+CABLE_DEF_TYPES_MOD >   ¾]     a   DEFINE_AIR%MET_TYPE%TVAIR+CABLE_DEF_TYPES_MOD >   R^     a   DEFINE_AIR%MET_TYPE%TVRAD+CABLE_DEF_TYPES_MOD <   æ^     a   DEFINE_AIR%MET_TYPE%PMB+CABLE_DEF_TYPES_MOD ;   z_     a   DEFINE_AIR%MET_TYPE%UA+CABLE_DEF_TYPES_MOD ;   `     a   DEFINE_AIR%MET_TYPE%QV+CABLE_DEF_TYPES_MOD >   ¢`     a   DEFINE_AIR%MET_TYPE%QVAIR+CABLE_DEF_TYPES_MOD ;   6a     a   DEFINE_AIR%MET_TYPE%DA+CABLE_DEF_TYPES_MOD <   Êa     a   DEFINE_AIR%MET_TYPE%DVA+CABLE_DEF_TYPES_MOD ?   ^b     a   DEFINE_AIR%MET_TYPE%COSZEN+CABLE_DEF_TYPES_MOD =   òb     a   DEFINE_AIR%MET_TYPE%NDEP+CABLE_DEF_TYPES_MOD =   c     a   DEFINE_AIR%MET_TYPE%PDEP+CABLE_DEF_TYPES_MOD <   d  ¬   a   DEFINE_AIR%MET_TYPE%FSD+CABLE_DEF_TYPES_MOD 2   Æd  @     DEFINE_AIR%MP+CABLE_DEF_TYPES_MOD    e  a   a   DEFINE_AIR%MET    ge  a   a   DEFINE_AIR%AIR 