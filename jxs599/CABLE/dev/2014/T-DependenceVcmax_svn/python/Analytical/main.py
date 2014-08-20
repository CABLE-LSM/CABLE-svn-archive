#!/usr/bin/python #-i
# use -i here to drop into innteractive mode t end of script
# and preserve namespce. have to comment out below as well

__author__ = 'Jhan Srbinovsky'

# main program governing project to do .... 

# usage:

#import python modules
import sys
import matplotlib.pyplot as plt
import numpy as np

# End python header
#########################################################################

# commented to run interactively
#def main(argv):

# GLOBAL decs so can plot models on same domain 
from main_data import T_leaf#, KK_pl, KK_bi, CABLE_c3 
from main_data import T0C_degK, nT_L

#########################################################################

#def test_fn( K ):
#   print K
#   K = [1,2,3,4,5]
#   print K
#   return K
#########################################################################

# insert a break from the CLI for reading output
print "\n"

#KK_pl    = []
#KK_bi    = []
CABLE_c3 = []
#K = test_fn( KK_pl ) 
#print K
#sys.exit()
#CABLE_c3 = np.zeros( (nT_L) ) 

KK_pl    = np.zeros( (15, nT_L) ) 
KK_bi    = np.zeros( (nT_L) ) 
# call model(s)

#import local, application specific modules
from Tdep_KK_main import Tdep_KK_main

# Dataset to read 
ifile = "KK_dataset.txt"

# KK Temperature dependence model
KK_pl, KK_bi = Tdep_KK_main( T_leaf, ifile )
#KK_bi = Tdep_KK_main( T_leaf, ifile )

#########################################################################

from Tdep_CABLE import T_dep_CABLE_C3 
CABLE_c3 = T_dep_CABLE_C3( T_leaf, ifile, CABLE_c3 ) 

#########################################################################

# plot CABLE against KK

show_plot = False
show_plot = True

x = T_leaf - T0C_degK 

if show_plot is True:   
   plt.title('Temp-dep per plant - KK Biome 1 vs CABLE c3')
   plt.ylabel('Vcmax')

   # plot Vcmax per plant in first Biome group 
   for j in range ( 0,14 ):
      y3 = KK_pl[j]
      plt.plot(x, y3, 'g-', linewidth=1 )
   
   # generic desc of T-dependece
   y1= KK_bi
   plt.plot(x, KK_bi, 'r-', linewidth=2 )
   y1 = CABLE_c3 
   plt.plot(x, y1, 'b-', linewidth=2 )
   #plt.show()    
   #
   plt.savefig("KKvCABLE.pdf")
   plt.close() 
 
##################














#
## plot model(s)
#
## plot generically x against y
##plot_generic( model[0].x, model[0].y )
#
## plot KK instance for T-dependence - averaged params other than PFT dep. Vcmax 
#x= model[0].x - T0C_degK
#
#y1= model[1].y[0]      
##y1= np.log(model[1].y[0] )     
#plt.subplot(3, 1, 1)
#plt.plot(x, y1, 'g-', linewidth=1 )
#plt.title('KK - per biome')
#plt.ylabel('Vcmax/Vcmax_25')
#
## broadleaf
#y2= model[0].y[0]     
#plt.subplot(3, 1, 2)
#plt.plot(x, y2, 'g-', linewidth=1 )
#
##coniferous
#y2= model[0].y[1]     
#plt.subplot(3, 1, 2)
#plt.plot(x, y2, 'r-', linewidth=1 )
#
##herbaceous
#y2= model[0].y[2]     
#plt.subplot(3, 1, 2)
#plt.plot(x, y2, 'b-', linewidth=1 )
#
#y3= model[2].y[0]     
#plt.subplot(3, 1, 3)
#plt.plot(x, y3, 'b-', linewidth=1 )
#
#y3= model[2].y[1]     
#plt.subplot(3, 1, 3)
#plt.plot(x, y3, 'b-', linewidth=1 )
#
#y3= model[2].y[2]     
#plt.subplot(3, 1, 3)
#plt.plot(x, y3, 'b-', linewidth=1 )
#
##for i in range ( nPlants[0] ):
#for i in range ( 17 ):
#   y3= model[2].y[i]     
#   plt.subplot(3, 1, 3)
#   plt.plot(x, y3, 'b-', linewidth=1 )
#
#y2= model[1].y[0]  * bi[0].Vcmax_25[0]    
#plt.subplot(3, 1, 3)
#plt.plot(x, y2, 'r-', linewidth=1 )
#
###for i in range ( nPlants[0] ):
##j =99
##y3= model[2].y[j]     
##plt.subplot(3, 1, 3)
##plt.plot(x, y3, 'y-', linewidth=1 )
##
#plt.xlabel('Leaf Temperature (deg C)')
#plt.ylabel('Vcmax')
#
#plt.show()
#
#########################################################################


#plt.plot( x, y, 'bd' )
#
#
##this is first dim(biome)=0
#y= model[0].y[0]     
#plt.plot( x, y, '-', linewidth=1 )
#
##sys.exit()
## plot Second Biome, T-dependence KK 
#y= model[0].y[1]     #this is first dim(biome)=1
#plt.plot( x, y, 'r-', linewidth=1 )
#
###plt.plot( x, y, 'rs' )
#
## plot Third Biome, T-dependence KK 
#y= model[0].y[2]     #this is first dim(biome)=2
#plt.plot( x, y, 'g-', linewidth=1 )
###plt.plot( x, y, 'go' )
#
#
#plt.show()
#######################

#if running purely as a script from the command line
################################################################################
#if __name__ == "__main__":
#   main(sys.argv[1:])

################################################################################

