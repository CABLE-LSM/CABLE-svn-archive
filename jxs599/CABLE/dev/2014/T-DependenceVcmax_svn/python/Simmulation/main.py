#!/usr/bin/python #-i
# use -i here to drop into innteractive mode t end of script
# and preserve namespce. have to comment out below as well

__author__ = 'Jhan Srbinovsky'

# main program governing project to do .... 

# usage:
from main_data import CABLEraw, CABLEdiurnal, recl
#import python modules
import sys
import matplotlib.pyplot as plt
import numpy as np

# End python header
#########################################################################

# commented to run interactively
#def main(argv):

# GLOBAL decs so can plot models on same domain 
#from main_data import T_leaf#, KK_pl, KK_bi, CABLE_c3 

#########################################################################

# insert a break from the CLI for reading output
print "\n"

#import local, application specific modules
from DataSet import ReadBinary


# Dataset to read 
ifile = "Vcmax_sunlit000.bin"

# KK Temperature dependence model
ReadBinary( ifile )

##########################################################################
show_plot = False
show_plot = True

x = np.zeros( len(CABLEraw) ) 

if show_plot is True:   
   plt.title('Vcmax out of CABLE c3 for Tumbarumba')
   plt.ylabel('Vcmax')

   # plot Vcmax per plant in first Biome group 
   for j in range (  len(CABLEraw) -1 ):
      x[j] = j
   plt.plot(x, CABLEraw, 'g-', linewidth=1 )
   
   plt.show()    
   plt.close() 

   x = np.zeros( len(CABLEraw) / 24 ) 
   j=0
   for i in range ( 0, len(CABLEraw)-1, 24 ):
      x[j] =j 
      CABLEdiurnal[j] = ( np.sum( CABLEraw[i:i+23] ) /24)
      print x[j],  CABLEdiurnal[j] 
      j=j+1
   plt.plot(x, CABLEdiurnal, 'r-', linewidth=1 )

#   # generic desc of T-dependece
#   y1= KK_bi
#   plt.plot(x, KK_bi, 'r-', linewidth=2 )
#   y1 = CABLE_c3 
#   plt.plot(x, y1, 'b-', linewidth=2 )
   plt.show()    
   #
#   plt.savefig("KKvCABLE.pdf")
   plt.close() 
 
## plot CABLE against KK
#
#show_plot = False
#show_plot = True
#
#x = T_leaf - T0C_degK 
#
#if show_plot is True:   
#   plt.title('Temp-dep per plant - KK Biome 1 vs CABLE c3')
#   plt.ylabel('Vcmax')
#
#   # plot Vcmax per plant in first Biome group 
#   for j in range ( 0,14 ):
#      y3 = KK_pl[j]
#      plt.plot(x, y3, 'g-', linewidth=1 )
#   
#   # generic desc of T-dependece
#   y1= KK_bi
#   plt.plot(x, KK_bi, 'r-', linewidth=2 )
#   y1 = CABLE_c3 
#   plt.plot(x, y1, 'b-', linewidth=2 )
#   #plt.show()    
#   #
#   plt.savefig("KKvCABLE.pdf")
#   plt.close() 
 
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

