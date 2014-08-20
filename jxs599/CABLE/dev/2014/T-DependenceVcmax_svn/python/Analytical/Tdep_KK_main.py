#!/usr/bin/python 

# main program governing iT-dependence Kattage & Knorr (KK) 

# usage: called from main

#import python modules
import numpy as np 
import sys
import matplotlib.pyplot as plt

#import local, application specific modules
from Tdep_KK_funcs import Tdep_KK_func, Tdep_KK_func_pl
from Biomes_KK import Set_biomes
from Speices_KK import Set_plants
from DataSet import Read_dataset
from main_data import PlantDepParams_KK, BiomeDepParams_KK
from main_data import pl, bi, T0C_degK, R_gas, nT_L
from main_data import nPlants_max, nBiomes_max

def Tdep_KK_main( T_leaf, ifile ):
   
   # arrays to be returned to main for plotting
   KK_pl    = np.zeros( (15, nT_L) ) 
   KK_bi    = np.zeros( (nT_L) ) 
   
   # insert a break from the CLI for reading output
   print "\n"
   
   # f^n, model parameter: Vcmax_25 reference, normalization Temperature
   T_ref = T0C_degK + 25. # 25 deg-C #hardwired 0 deg-C (K)

#########################################################################

   # Declare mutable object => lines read from config file
   nlines = []
   data = Read_dataset( nlines, ifile )

#########################################################################
 
   # Fix nBiomes[i], nPlants[i] from dataSet
   nBiomes     = [] 
   nPlants     = [] 

   # instantiate classes: for plants
   pl = PlantDepParams_KK()

   # Set Plant dependent Parameters
   
   # hardwired 0 ~ bc 0th model,dataset here
   nPlants.append( Set_plants( nlines[0], nBiomes, pl, data ) )

#########################################################################

   # instantiate classes: for biomes
   for i in range( nBiomes[0] ):
      bi.append( BiomeDepParams_KK() )
   
   # Set Biome dependent Parameters
   
   Set_biomes( nBiomes[0], bi, nPlants[0], pl ) 
   
#########################################################################
   ################################
   # Temperature depence function #
   ################################

   # Vcmax & Jmax described by KK(2007) Eq.(1)
   # Kmax = K_25 * f(T_leaf) [ described by fn_T_L here ]

   # eveluate f(T_leaf) per biome 
   fn_T_L = np.zeros((nBiomes[0], len( T_leaf ) ))

   # Call function defining Temperature dependence 
   # with average paramater values for all plants, biomes 
   for j in range( nBiomes[0] ):
      for i in range( len( T_leaf ) ):
         fn_T_L[j][i] = Tdep_KK_func( j, T_leaf[i], T_ref, R_gas, bi, pl )
  
#########################################################################
   
   ################################
   # Models/plots from KK dataset #
   ################################

   # import library to handle writing plots to pdfs
   from matplotlib.backends.backend_pdf import PdfPages
   
#########################################################################
   
   ################################
   # Models/plots from KK dataset #
   ################################

   # Reproduce general fn(T_leaf) wout T acclim, as per Figure 4c in KK 2007 
   # Average Params H_a, DeltaS, Vcmax_25 over biome and then normalize
   # by Vcmax_25
   
   # quick switch
   show_plot = False
   show_plot = True

   # this will endure over following plots = leaf temperature
   x = T_leaf - T0C_degK

   Vcmax_KK_Bi = np.zeros( ( nBiomes[0], len(T_leaf) ) )
   # Vcmax & Jmax described by KK(2007) Eq.(1)
   for j in range( nBiomes[0] ):
      for i in range( len( T_leaf ) ):
         Vcmax_KK_Bi[j][i]= fn_T_L[j][i]
         KK_bi[i] = fn_T_L[0][i]

   Ksum = KK_bi 

   # This is: Vcmax_25_Bi[0][:], j=0th biome plotted only as all are normaized 
   # by corresponding Vcmax_25 (bi[0].Vcmax_25[j] anyway) and so lie on top of 
   # each other.

   y1= Vcmax_KK_Bi[0]
   
   if show_plot is True:   
      plt.title('KK - per biome')
      plt.ylabel('Vcmax/Vcmax_25')
      plot1 = plt.plot(x, y1, 'r-', linewidth=1 )
     
      #pp = PdfPages('KK.pdf')
      plt.savefig("KK_reproTheirFig4c.pdf")
      plt.close() 

#########################################################################

   # Average Params H_a, DeltaS, Vcmax_25 over biome BUT dont normalize
   # by Vcmax_25
   
   # quick switch
   show_plot = False
   show_plot = True

   # Vcmax & Jmax described by KK(2007) Eq.(1)
   for j in range( nBiomes[0] ):
      for i in range( len( T_leaf ) ):
         Vcmax_KK_Bi[j][i]= bi[0].Vcmax_25[j] * fn_T_L[j][i]
   
   if show_plot is True:   
      plt.title('KK Temp. dependence per biome')
      plt.ylabel('Vcmax')
      
      # Biome 1
      y1= Vcmax_KK_Bi[0]
      plot2 = plt.plot(x, y1, 'r-', linewidth=1 )

      # Biome 2
      y1= Vcmax_KK_Bi[1]
      plot3 = plt.plot(x, y1, 'b-', linewidth=1 )

      # Biome 3
      y1= Vcmax_KK_Bi[2]
      plot4 = plt.plot(x, y1, 'g-', linewidth=1 )

   plt.savefig("KK_TheirFig4c_NOnorm.pdf")
   plt.close() 
   #plt.show()
  
#########################################################################

   # Per species params H_a, DeltaS, Vcmax_25 for each biome AND is normalized
   # by Vcmax_25
   
   # quick switch
   show_plot = False
   show_plot = True

   # eveluate f(T_leaf) per plant 
   gn_T_L = np.zeros((nPlants[0], len( T_leaf ) ))

   for j in range( nPlants[0] ):
      for i in range( len( T_leaf ) ):
         gn_T_L[j][i] = Tdep_KK_func_pl( j, T_leaf[i], T_ref, R_gas, bi,pl )
  
   Vcmax_KK_pl = np.zeros( ( nPlants[0], len(T_leaf) ) )
   # Vcmax per plant 
   for j in range( nPlants[0] ):
      for i in range( len( T_leaf ) ):
         Vcmax_KK_pl[j][i]= gn_T_L[j][i]

   if show_plot is True:   
      # plot Vcmax per plant in first Biome group 
      plt.title('KK Temp-dep per plant - Biome 1')
      plt.ylabel('Vcmax/Vcmax_25')
      
      for j in range ( 0,14 ):
         y3= Vcmax_KK_pl[j]
         plt.plot(x, y3, 'g-', linewidth=1 )
         
      # generic desc of T-dependece
      Vcmax_KK_Bi[0]= fn_T_L[0]
      y1= Vcmax_KK_Bi[0]
      plt.plot(x, y1, 'r-', linewidth=2 )
      
      plt.savefig("KK_Norm_perSpecies1.pdf")
      plt.close() 
      #plt.show()
 
##################

      plt.title('KK Temp-dep per plant - Biome 2')
      # plot Vcmax per plant in 2nd Biome group 
      for j in range ( 15,23 ):
         y3= Vcmax_KK_pl[j]
         plt.plot(x, y3, 'g-', linewidth=1 )
         
      # generic desc of T-dependece
      Vcmax_KK_Bi[1]= fn_T_L[1]
      y1= Vcmax_KK_Bi[1]
      plt.plot(x, y1, 'r-', linewidth=2 )
      
      plt.savefig("KK_Norm_perSpecies2.pdf")
      plt.close() 
      #plt.show()
 
##################

      plt.title('KK Temp-dep per plant - Biome 3')
      # plot Vcmax per plant in third Biome group 
      for j in range ( 24,47 ):
         y3= Vcmax_KK_pl[j]
         plt.plot(x, y3, 'g-', linewidth=1 )
         
      # generic desc of T-dependece
      Vcmax_KK_Bi[2]= fn_T_L[2]
      y1= Vcmax_KK_Bi[2]
      plt.plot(x, y1, 'r-', linewidth=2 )
      
      plt.savefig("KK_Norm_perSpecies3.pdf")
      plt.close() 
      #plt.show()
   
##################

#########################################################################
 
   # Per species params H_a, DeltaS, Vcmax_25 for each biome BUT dont normalize
   # by Vcmax_25
   
   # quick switch
   show_plot = False
   show_plot = True

   # plot Un-Normed Vcmax per Biome group 
   Vcmax_KK_Bi[0]= bi[0].Vcmax_25[0] * fn_T_L[0]
   Vcmax_KK_Bi[1]= bi[0].Vcmax_25[1] * fn_T_L[1]
   Vcmax_KK_Bi[2]= bi[0].Vcmax_25[2] * fn_T_L[2]

   # eveluate f(T_leaf) per plant 
   gn_T_L = np.zeros((nPlants[0], len( T_leaf ) ))

   for j in range( nPlants[0] ):
      for i in range( len( T_leaf ) ):
         gn_T_L[j][i] = Tdep_KK_func_pl( j, T_leaf[i], T_ref, R_gas, bi,pl )
  
   Vcmax_KK_pl = np.zeros( ( nPlants[0], len(T_leaf) ) )
   # Vcmax per plant 
   for j in range( nPlants[0] ):
      for i in range( len( T_leaf ) ):
         Vcmax_KK_pl[j][i]= pl.Vcmax_25[j] * gn_T_L[j][i]

   if show_plot is True:   
      plt.title('KK Temp-dep per plant - Biome 1')
      plt.ylabel('Vcmax')
      # plot Vcmax per plant in first Biome group 
      for j in range ( 0,14 ):
         y3= Vcmax_KK_pl[j]
         plt.plot(x, y3, 'g-', linewidth=1 )
         
      # generic desc of T-dependece
      y1= Vcmax_KK_Bi[0]
      plt.plot(x, y1, 'r-', linewidth=2 )
      
      plt.savefig("KK_perSpecies1.pdf")
      plt.close() 
 
##################

      plt.title('KK Temp-dep per plant - Biome 2')
      # plot Vcmax per plant in 2nd Biome group 
      for j in range ( 15,23 ):
         y3= Vcmax_KK_pl[j]
         plt.plot(x, y3, 'g-', linewidth=1 )
         
      # generic desc of T-dependece
      y1= Vcmax_KK_Bi[1]
      plt.plot(x, y1, 'r-', linewidth=2 )
      
      plt.savefig("KK_perSpecies2.pdf")
      plt.close() 
      
##################

      plt.title('KK Temp-dep per plant - Biome 3')
      # plot Vcmax per plant in third Biome group 
      for j in range ( 24,47 ):
         y3= Vcmax_KK_pl[j]
         plt.plot(x, y3, 'g-', linewidth=1 )
         
      # generic desc of T-dependece
      y1= Vcmax_KK_Bi[2]
      plt.plot(x, y1, 'r-', linewidth=2 )
      
      plt.savefig("KK_perSpecies3.pdf")
      plt.close() 
   
##################

   # plot Vcmax per plant in first Biome group 
   for j in range ( 0,14 ):
      KK_pl[j] = Vcmax_KK_pl[j] / pl.Vcmax_25[j] 
   
   #print "KK_bi ", KK_pl 
   return( KK_pl, Ksum ) 
   #return( Ksum  ) 
 
##################


   
    
   # Vcmax & Jmax described by KK(2007) Eq.(1)
   #for j in range( nBiomes[0] ):
   #j=99
   #for i in range( len( T_leaf ) ):
   #   Vcmax[j][i]= bi[0].Vcmax_25[j] * fn_T_L[j][i]
   #   # normalized by Vcmax_25 per Biome
 
   #for i in range( len( T_leaf ) ):
   #   print "Vcmax_25 ", Vcmax[2][i] -  Vcmax[1][i]

   #for i in range( nBiomes):
   #   for j in range( len( T_leaf ) ):
   #      Vcmax = Vcmax_25[i] * fn_T_L[j]


################################################################################



























   # Vcmax & Jmax described by KK(2007) Eq.(1)
   for j in range( nBiomes[0] ):
      for i in range( len( T_leaf ) ):
         Vcmax_KK_Bi[j][i]= bi[0].Vcmax_25[j] * fn_T_L[j][i]
         # normalized by Vcmax_25 per Biome
         #Vcmax_KK_Bi[j][i]= fn_T_L[j][i]
  
   # eveluate f(T_leaf) per biome 
   gn_T_L = np.zeros((nPlants[0], len( T_leaf ) ))

   ## Call function 
   for j in range( nPlants[0] ):
      for i in range( len( T_leaf ) ):
         gn_T_L[j][i] = Tdep_KK_func_pl( j, T_leaf[i], T_ref, R_gas, bi,pl )
  
   # Vcmax & Jmax described by KK(2007) Eq.(1)
   for j in range( nPlants[0] ):
      for i in range( len( T_leaf ) ):
         Vcmax_KK_pl[j][i]= pl.Vcmax_25[j] * gn_T_L[j][i]

   # Vcmax & Jmax described by KK(2007) Eq.(1)
   #for j in range( nBiomes[0] ):
   #j=99
   #for i in range( len( T_leaf ) ):
   #   Vcmax[j][i]= bi[0].Vcmax_25[j] * fn_T_L[j][i]
   #   # normalized by Vcmax_25 per Biome
 
   #for i in range( len( T_leaf ) ):
   #   print "Vcmax_25 ", Vcmax[2][i] -  Vcmax[1][i]

   #for i in range( nBiomes):
   #   for j in range( len( T_leaf ) ):
   #      Vcmax = Vcmax_25[i] * fn_T_L[j]


################################################################################



























