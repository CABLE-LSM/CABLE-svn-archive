  �  6   k820309    9          19.0        �=c                                                                                                          
       cable_weathergenerator.F90 CABLE_WEATHERGENERATOR              SP PI PIBY2 SECDAY SOLARCONST EPSILON SBOLTZ                   @               @                '�	             &      #NP    #NDTIME    #DELT    #LATDEG    #WINDDAY    #TEMPMINDAY    #TEMPMAXDAY    #TEMPMINDAYNEXT 	   #TEMPMAXDAYPREV 
   #SOLARMJDAY    #DECRAD    #WINDDARK    #WINDLITE    #SOLARNORM    #LATRAD    #DAYLENGTH    #TIMESUNSETPREV    #TIMESUNRISE    #TIMEMAXTEMP    #TIMESUNSET    #TEMPSUNSETPREV    #TEMPSUNSET    #TEMPNIGHTRATE    #TEMPNIGHTRATEPREV    #TEMPRANGEDAY    #TEMPRANGEAFT    #PRECIPDAY    #SNOWDAY    #PMBDAY    #PHISD    #PHILD     #PRECIP !   #SNOW "   #WIND #   #TEMP $   #VAPPMB %   #PMB &   #COSZEN '                �                                                               �                                                              �                                              	              �                                                           
            &                                                      �                                          X                 
            &                                                      �                                          �                 
            &                                                      �                                          �                 
            &                                                      �                              	            0                
            &                                                      �                              
            x             	   
            &                                                      �                                          �             
   
            &                                                        �                                            
              �                                                          
            &                                                      �                                          X                
            &                                                      �                                          �                
            &                                                      �                                          �                
            &                                                      �                                          0                
            &                                                      �                                          x                
            &                                                      �                                          �                
            &                                                      �                                                          
            &                                                      �                                          P                
            &                                                      �                                          �                
            &                                                      �                                          �                
            &                                                      �                                          (                
            &                                                      �                                          p                
            &                                                      �                                          �                
            &                                                      �                                                           
            &                                                      �                                          H                
            &                                                      �                                          �                
            &                                                      �                                          �                
            &                                                      �                                                           
            &                                                      �                                           h                
            &                                                      �                              !            �                 
            &                                                      �                              "            �             !   
            &                                                      �                              #            @             "   
            &                                                      �                              $            �             #   
            &                                                      �                              %            �             $   
            &                                                      �                              &            	             %   
            &                                                      �                              '            `	             &   
            &                                           #         @                                   (                    #WG )   #NP *   #LATITUDE +   #DELS ,             D                                 )     �	              #WEATHER_GENERATOR_TYPE              
                                  *                    
                                  +                    	 #   p          5 � p        r *       5 � p        r *                               
                                  ,     	      #         @                                   -                    #WG .   #NP /   #YEARDAY 0             D                                 .     �	              #WEATHER_GENERATOR_TYPE              
                                  /                     
                                  0           #         @                                   1                    #WG 2   #NP 3   #ITIME 4             D                                 2     �	              #WEATHER_GENERATOR_TYPE              
                                  3                     
  @                               4              �   :      fn#fn ,   �   =   b   uapp(CABLE_WEATHERGENERATOR '     r      WEATHER_GENERATOR_TYPE *   �  H   a   WEATHER_GENERATOR_TYPE%NP .   �  H   a   WEATHER_GENERATOR_TYPE%NDTIME ,     H   a   WEATHER_GENERATOR_TYPE%DELT .   a  �   a   WEATHER_GENERATOR_TYPE%LATDEG /   �  �   a   WEATHER_GENERATOR_TYPE%WINDDAY 2   �  �   a   WEATHER_GENERATOR_TYPE%TEMPMINDAY 2     �   a   WEATHER_GENERATOR_TYPE%TEMPMAXDAY 6   �  �   a   WEATHER_GENERATOR_TYPE%TEMPMINDAYNEXT 6   E  �   a   WEATHER_GENERATOR_TYPE%TEMPMAXDAYPREV 2   �  �   a   WEATHER_GENERATOR_TYPE%SOLARMJDAY .   m  H   a   WEATHER_GENERATOR_TYPE%DECRAD 0   �  �   a   WEATHER_GENERATOR_TYPE%WINDDARK 0   I	  �   a   WEATHER_GENERATOR_TYPE%WINDLITE 1   �	  �   a   WEATHER_GENERATOR_TYPE%SOLARNORM .   q
  �   a   WEATHER_GENERATOR_TYPE%LATRAD 1     �   a   WEATHER_GENERATOR_TYPE%DAYLENGTH 6   �  �   a   WEATHER_GENERATOR_TYPE%TIMESUNSETPREV 3   -  �   a   WEATHER_GENERATOR_TYPE%TIMESUNRISE 3   �  �   a   WEATHER_GENERATOR_TYPE%TIMEMAXTEMP 2   U  �   a   WEATHER_GENERATOR_TYPE%TIMESUNSET 6   �  �   a   WEATHER_GENERATOR_TYPE%TEMPSUNSETPREV 2   }  �   a   WEATHER_GENERATOR_TYPE%TEMPSUNSET 5     �   a   WEATHER_GENERATOR_TYPE%TEMPNIGHTRATE 9   �  �   a   WEATHER_GENERATOR_TYPE%TEMPNIGHTRATEPREV 4   9  �   a   WEATHER_GENERATOR_TYPE%TEMPRANGEDAY 4   �  �   a   WEATHER_GENERATOR_TYPE%TEMPRANGEAFT 1   a  �   a   WEATHER_GENERATOR_TYPE%PRECIPDAY /   �  �   a   WEATHER_GENERATOR_TYPE%SNOWDAY .   �  �   a   WEATHER_GENERATOR_TYPE%PMBDAY -     �   a   WEATHER_GENERATOR_TYPE%PHISD -   �  �   a   WEATHER_GENERATOR_TYPE%PHILD .   E  �   a   WEATHER_GENERATOR_TYPE%PRECIP ,   �  �   a   WEATHER_GENERATOR_TYPE%SNOW ,   m  �   a   WEATHER_GENERATOR_TYPE%WIND ,     �   a   WEATHER_GENERATOR_TYPE%TEMP .   �  �   a   WEATHER_GENERATOR_TYPE%VAPPMB +   )  �   a   WEATHER_GENERATOR_TYPE%PMB .   �  �   a   WEATHER_GENERATOR_TYPE%COSZEN    Q  p       WGEN_INIT    �  d   a   WGEN_INIT%WG    %  @   a   WGEN_INIT%NP #   e  �   a   WGEN_INIT%LATITUDE      @   a   WGEN_INIT%DELS %   Y  e       WGEN_DAILY_CONSTANTS (   �  d   a   WGEN_DAILY_CONSTANTS%WG (   "  @   a   WGEN_DAILY_CONSTANTS%NP -   b  @   a   WGEN_DAILY_CONSTANTS%YEARDAY $   �  c       WGEN_SUBDIURNAL_MET '     d   a   WGEN_SUBDIURNAL_MET%WG '   i  @   a   WGEN_SUBDIURNAL_MET%NP *   �  @   a   WGEN_SUBDIURNAL_MET%ITIME 