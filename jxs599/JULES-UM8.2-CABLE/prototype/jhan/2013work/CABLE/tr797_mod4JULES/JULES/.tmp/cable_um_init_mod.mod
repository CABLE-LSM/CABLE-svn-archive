	  �E  O   k820309    �          13.0        �S                                                                                                           
       cable_jules_init.F90 CABLE_UM_INIT_MOD #         @                                                   2   #EXPLICIT_CALL_INITIALIZATION%SUM    #EXPLICIT_CALL_INITIALIZATION%SIN    #EXPLICIT_CALL_INITIALIZATION%ALLOCATED    #ROW_LENGTH    #ROWS    #LAND_PTS    #NTILES    #NPFT 	   #SM_LEVELS 
   #TIMESTEP    #LATITUDE    #LONGITUDE    #LAND_INDEX    #TILE_FRAC    #TILE_PTS    #TILE_INDEX    #DZSOIL    #BEXP    #HCON    #SATCON    #SATHH    #SMVCST    #SMVCWT    #SMVCCL    #ALBSOIL    #CANHT_FT    #LAI_FT    #SW_DOWN    #LW_DOWN    #LS_RAIN    #LS_SNOW     #TL_1 !   #QW_1 "   #VSHR_LAND #   #PSTAR $   #Z1_TQ %   #Z1_UV &   #CANOPY_TILE '   #FLAND (   #CO2_MMR )   #COS_ZENITH_ANGLE *   #SNOW_TILE +   #SNAGE_TILE ,   #SNOW_RHO1L -   #ISNOW_FLG3L .   #SNOW_RHO3L /   #SNOW_DEPTH3L 0   #SNOW_TMP3L 1   #SNOW_MASS3L 2   #SNOW_COND 3   #SMCL_TILE 4   #STHF_TILE 5   #TSOIL_TILE 6                                                   SUM                                                 SIN                                                 ALLOCATED           
  @                                                    
  @                                                    
  @                                                    
  @                                                    
  @                               	                     
  @                               
                     
  @                                    	               
  @                                                   	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                   	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                       p          5 � p        r        5 � p        r                               
  @                                                   	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                       p          5 � p        r        5 � p        r                               
  @                                                         p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                   	    p          5 � p        r 
       5 � p        r 
                              
  @                                                   	    p          5 � p        r        5 � p        r                               
  @                                                   	    p          5 � p        r        5 � p        r                               
  @                                                   	    p          5 � p        r        5 � p        r                               
  @                                                   	 	   p          5 � p        r        5 � p        r                               
  @                                                   	 
   p          5 � p        r        5 � p        r                               
  @                                                   	    p          5 � p        r        5 � p        r                               
  @                                                   	    p          5 � p        r        5 � p        r                               
  @                                                   	    p          5 � p        r        5 � p        r                               
  @                                                   	       p        5 � p        r    p          5 � p        r      5 � p        r 	       5 � p        r      5 � p        r 	                              
  @                                                   	 !     p        5 � p        r    p          5 � p        r      5 � p        r 	       5 � p        r      5 � p        r 	                              
D @                                                   	       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                   	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                   	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                                                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               !                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               "                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               #                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               $                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               %                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               &                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               '                    	 "     p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               (                    	    p          5 � p        r        5 � p        r                                
  @                               )     	               
D @                               *                    	       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
D @                               +                    	       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               ,                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               -                    	      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               .                          p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               
  @                               /                    	 $       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      p            5 � p        r      5 � p        r      p                                   
  @                               0                    	 %       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      p            5 � p        r      5 � p        r      p                                   
  @                               1                    	 '       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      p            5 � p        r      5 � p        r      p                                   
  @                               2                    	 &       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      p            5 � p        r      5 � p        r      p                                   
D @                               3                    	 #        p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      p            5 � p        r      5 � p        r      p                                   
  @                               4                    	 )       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      5 � p        r 
       5 � p        r      5 � p        r      5 � p        r 
                              
  @                               5                    	 (       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      5 � p        r 
       5 � p        r      5 � p        r      5 � p        r 
                              
  @                               6                    	 *       p        5 � p        r    p        5 � p        r    p          5 � p        r      5 � p        r      5 � p        r 
       5 � p        r      5 � p        r      5 � p        r 
                     #         @                                  7                   #ASSIGN_UM_BASICS_TO_UM1%INT 8   #ROW_LENGTH 9   #ROWS :   #LAND_PTS ;   #NTILES <   #NPFT =   #SM_LEVELS >   #TIMESTEP ?   #LATITUDE @   #LONGITUDE A   #LAND_INDEX B   #TILE_FRAC C   #TILE_PTS D   #TILE_INDEX E                                              8     INT           
                                  9                     
                                  :                     
                                  ;                     
                                  <                     
                                  =                     
                                  >                     
  @                               ?     	               
                                  @                    	 8     p        5 � p        r 9   p          5 � p        r 9     5 � p        r :       5 � p        r 9     5 � p        r :                              
                                  A                    	 9     p        5 � p        r 9   p          5 � p        r 9     5 � p        r :       5 � p        r 9     5 � p        r :                              
                                  B                     :   p          5 � p        r ;       5 � p        r ;                              
                                  C                    	 =     p        5 � p        r ;   p          5 � p        r ;     5 � p        r <       5 � p        r ;     5 � p        r <                              
                                  D                     ;   p          5 � p        r <       5 � p        r <                              
                                  E                     <     p        5 � p        r ;   p          5 � p        r ;     5 � p        r <       5 � p        r ;     5 � p        r <                     #         @                                   F                    #ROW_LENGTH G   #ROWS H   #LS_RAIN I   #LS_SNOW J   #CONV_RAIN K   #CONV_SNOW L   #DTL_1 M   #DQW_1 N             
                                  G                     
                                  H                    
                                  I                    	 .     p        5 � p        r G   p          5 � p        r G     5 � p        r H       5 � p        r G     5 � p        r H                              
                                  J                    	 /     p        5 � p        r G   p          5 � p        r G     5 � p        r H       5 � p        r G     5 � p        r H                                                               K                    	 0      p        5 � p        r G   p          5 � p        r G     5 � p        r H       5 � p        r G     5 � p        r H                                                               L                    	 1      p        5 � p        r G   p          5 � p        r G     5 � p        r H       5 � p        r G     5 � p        r H                               @                               M                    	 2      p        5 � p        r G   p          5 � p        r G     5 � p        r H       5 � p        r G     5 � p        r H                               @                               N                    	 3      p        5 � p        r G   p          5 � p        r G     5 � p        r H       5 � p        r G     5 � p        r H                        �   /      fn#fn -   �   n      EXPLICIT_CALL_INITIALIZATION 1   =  <      EXPLICIT_CALL_INITIALIZATION%SUM 1   y  <      EXPLICIT_CALL_INITIALIZATION%SIN 7   �  B      EXPLICIT_CALL_INITIALIZATION%ALLOCATED 8   �  @   a   EXPLICIT_CALL_INITIALIZATION%ROW_LENGTH 2   7  @   a   EXPLICIT_CALL_INITIALIZATION%ROWS 6   w  @   a   EXPLICIT_CALL_INITIALIZATION%LAND_PTS 4   �  @   a   EXPLICIT_CALL_INITIALIZATION%NTILES 2   �  @   a   EXPLICIT_CALL_INITIALIZATION%NPFT 7   7  @   a   EXPLICIT_CALL_INITIALIZATION%SM_LEVELS 6   w  @   a   EXPLICIT_CALL_INITIALIZATION%TIMESTEP 6   �  $  a   EXPLICIT_CALL_INITIALIZATION%LATITUDE 7   �  $  a   EXPLICIT_CALL_INITIALIZATION%LONGITUDE 8   �  �   a   EXPLICIT_CALL_INITIALIZATION%LAND_INDEX 7   �	  $  a   EXPLICIT_CALL_INITIALIZATION%TILE_FRAC 6   �
  �   a   EXPLICIT_CALL_INITIALIZATION%TILE_PTS 8   �  $  a   EXPLICIT_CALL_INITIALIZATION%TILE_INDEX 4   �  �   a   EXPLICIT_CALL_INITIALIZATION%DZSOIL 2   c  �   a   EXPLICIT_CALL_INITIALIZATION%BEXP 2     �   a   EXPLICIT_CALL_INITIALIZATION%HCON 4   �  �   a   EXPLICIT_CALL_INITIALIZATION%SATCON 3     �   a   EXPLICIT_CALL_INITIALIZATION%SATHH 4   3  �   a   EXPLICIT_CALL_INITIALIZATION%SMVCST 4   �  �   a   EXPLICIT_CALL_INITIALIZATION%SMVCWT 4   �  �   a   EXPLICIT_CALL_INITIALIZATION%SMVCCL 5   O  �   a   EXPLICIT_CALL_INITIALIZATION%ALBSOIL 6     $  a   EXPLICIT_CALL_INITIALIZATION%CANHT_FT 4   '  $  a   EXPLICIT_CALL_INITIALIZATION%LAI_FT 5   K  $  a   EXPLICIT_CALL_INITIALIZATION%SW_DOWN 5   o  $  a   EXPLICIT_CALL_INITIALIZATION%LW_DOWN 5   �  $  a   EXPLICIT_CALL_INITIALIZATION%LS_RAIN 5   �  $  a   EXPLICIT_CALL_INITIALIZATION%LS_SNOW 2   �  $  a   EXPLICIT_CALL_INITIALIZATION%TL_1 2   �  $  a   EXPLICIT_CALL_INITIALIZATION%QW_1 7   #  $  a   EXPLICIT_CALL_INITIALIZATION%VSHR_LAND 3   G  $  a   EXPLICIT_CALL_INITIALIZATION%PSTAR 3   k  $  a   EXPLICIT_CALL_INITIALIZATION%Z1_TQ 3   �  $  a   EXPLICIT_CALL_INITIALIZATION%Z1_UV 9   �   $  a   EXPLICIT_CALL_INITIALIZATION%CANOPY_TILE 3   �!  �   a   EXPLICIT_CALL_INITIALIZATION%FLAND 5   �"  @   a   EXPLICIT_CALL_INITIALIZATION%CO2_MMR >   �"  $  a   EXPLICIT_CALL_INITIALIZATION%COS_ZENITH_ANGLE 7   �#  $  a   EXPLICIT_CALL_INITIALIZATION%SNOW_TILE 8   %  $  a   EXPLICIT_CALL_INITIALIZATION%SNAGE_TILE 8   7&  $  a   EXPLICIT_CALL_INITIALIZATION%SNOW_RHO1L 9   ['  $  a   EXPLICIT_CALL_INITIALIZATION%ISNOW_FLG3L 8   (  t  a   EXPLICIT_CALL_INITIALIZATION%SNOW_RHO3L :   �)  t  a   EXPLICIT_CALL_INITIALIZATION%SNOW_DEPTH3L 8   g+  t  a   EXPLICIT_CALL_INITIALIZATION%SNOW_TMP3L 9   �,  t  a   EXPLICIT_CALL_INITIALIZATION%SNOW_MASS3L 7   O.  t  a   EXPLICIT_CALL_INITIALIZATION%SNOW_COND 7   �/  �  a   EXPLICIT_CALL_INITIALIZATION%SMCL_TILE 7   W1  �  a   EXPLICIT_CALL_INITIALIZATION%STHF_TILE 8   �2  �  a   EXPLICIT_CALL_INITIALIZATION%TSOIL_TILE (   4        ASSIGN_UM_BASICS_TO_UM1 ,   �5  <      ASSIGN_UM_BASICS_TO_UM1%INT 3   �5  @   a   ASSIGN_UM_BASICS_TO_UM1%ROW_LENGTH -   6  @   a   ASSIGN_UM_BASICS_TO_UM1%ROWS 1   Y6  @   a   ASSIGN_UM_BASICS_TO_UM1%LAND_PTS /   �6  @   a   ASSIGN_UM_BASICS_TO_UM1%NTILES -   �6  @   a   ASSIGN_UM_BASICS_TO_UM1%NPFT 2   7  @   a   ASSIGN_UM_BASICS_TO_UM1%SM_LEVELS 1   Y7  @   a   ASSIGN_UM_BASICS_TO_UM1%TIMESTEP 1   �7  $  a   ASSIGN_UM_BASICS_TO_UM1%LATITUDE 2   �8  $  a   ASSIGN_UM_BASICS_TO_UM1%LONGITUDE 3   �9  �   a   ASSIGN_UM_BASICS_TO_UM1%LAND_INDEX 2   �:  $  a   ASSIGN_UM_BASICS_TO_UM1%TILE_FRAC 1   �;  �   a   ASSIGN_UM_BASICS_TO_UM1%TILE_PTS 3   m<  $  a   ASSIGN_UM_BASICS_TO_UM1%TILE_INDEX -   �=  �       IMPLICIT_CALL_INITIALIZATION 8   A>  @   a   IMPLICIT_CALL_INITIALIZATION%ROW_LENGTH 2   �>  @   a   IMPLICIT_CALL_INITIALIZATION%ROWS 5   �>  $  a   IMPLICIT_CALL_INITIALIZATION%LS_RAIN 5   �?  $  a   IMPLICIT_CALL_INITIALIZATION%LS_SNOW 7   	A  $  a   IMPLICIT_CALL_INITIALIZATION%CONV_RAIN 7   -B  $  a   IMPLICIT_CALL_INITIALIZATION%CONV_SNOW 3   QC  $  a   IMPLICIT_CALL_INITIALIZATION%DTL_1 3   uD  $  a   IMPLICIT_CALL_INITIALIZATION%DQW_1 