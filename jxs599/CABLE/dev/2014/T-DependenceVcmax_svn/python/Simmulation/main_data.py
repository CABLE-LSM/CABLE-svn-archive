#!/usr/bin/python

import numpy as np

#data = []

# length of data - # timesteps
Tumba = 17520 #hourly, 2 years

# declare numpy array to hold real data 
#data= np.zeros(recl/3)

recl=[]
CABLEraw= np.zeros(Tumba)
CABLEdiurnal= np.zeros(Tumba/24)
CABLEseas= np.zeros(Tumba/4)
   

# 0 deg-C in K
T0C_degK = 273.15
Tmin0_degC = 0.; Tmax0_degC = 40.
Tmin0 = Tmin0_degC + T0C_degK 
Tmax0 = Tmax0_degC + T0C_degK 

# "square" dimensions over models

# x-dimension: number of values of T_leaf
nT_L = 100 #; n_max = 200

# Set domain of Leaf Temperature (Kelvin) with nx grads from Tmin:Tmax
T_leaf = []
T_leaf = np.linspace(Tmin0, Tmax0, nT_L)

# Constants
R_gas = 8.3144621      # Universal Gas Constant 


#########################################################################
# Put KK dataset classes here so we can use same data in other models

# these maximums are more than sufficient
nPlants_max = 100; nBiomes_max = 5 

# KK(2007) actually gives values per plant
class PlantDepParams_KK(object):
   def __init__(self):
      self.PlNumber = np.zeros( nPlants_max ) 
      self.BiNumber = np.zeros( nPlants_max ) 
      self.Vcmax_25 = np.zeros( nPlants_max ) 
      self.Jmax_25  = np.zeros( nPlants_max ) 
      self.H_a      = np.zeros( nPlants_max )    # activation energy
      self.H_d      = np.zeros( nPlants_max )    # de-activation energy
      self.DeltaS   = np.zeros( nPlants_max )    # Entropy factor
 
class BiomeDepParams_KK(object):
   def __init__(self):
      self.BiNumber = np.zeros( nPlants_max ) 
      self.Vcmax_25 = np.zeros( nBiomes_max ) 
      self.Jmax_25  = np.zeros( nBiomes_max ) 
      self.H_a      = np.zeros( nBiomes_max )    # activation energy
      self.H_d      = np.zeros( nBiomes_max )    # de-activation energy
      self.DeltaS   = np.zeros( nBiomes_max )    # Entropy factor
      self.T_opt    = np.zeros( nBiomes_max )    # Optimum temp 
    
# instantiate classes: for plants
pl = []
bi = []
#########################################################################

#CABLE_c3 = []
#KK_pl    = []
#KK_bi    = []
 
## class containing x & y data to plot    
#class Models(object):
#   def __init__(self):
#      self.x = np.zeros(nT_L) 
#      self.y = np.zeros( (n_max, nT_L) ) # assumes1:1 mapping b/n x&y
# End global class declarations

