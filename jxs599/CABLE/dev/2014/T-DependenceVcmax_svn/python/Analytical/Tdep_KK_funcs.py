import sys
import math as m
import numpy as np

###############################################################################

def Tdep_KK_func( i, T_L, T_ref, R_gas, bi, pl ):
   
   # Activation moderated by deactivation
   
   # Activation 
   ####################
   DeltaT = T_L - T_ref
   
   Eff_activ  = float( bi[0].H_a[i] * DeltaT )
   norm_activ = float( T_ref * T_L * R_gas )
   
   Activation = Tdep_KK_funcExpFuncs( Eff_activ, norm_activ )
   
   # Deactivation: ratio b/n ref & leaf
   ####################
   
   # _ref Temperature
   Entropy_ref =  ( T_ref  * bi[0].DeltaS[i] ) - bi[0].H_d[i]
   norm_ref = R_gas * T_ref
   Deactiv_ref = 1. + Tdep_KK_funcExpFuncs( Entropy_ref, norm_ref )
   
   # leaf Temperature
   Entropy_leaf = ( T_L * bi[0].DeltaS[i] ) - bi[0].H_d[i]
   norm_leaf = R_gas * T_L
   Deactiv_leaf = 1. + Tdep_KK_funcExpFuncs( Entropy_leaf, norm_leaf )
  
   fn_T_L =  Activation * Deactiv_ref / Deactiv_leaf
   
   return fn_T_L

###############################################################################

def Tdep_KK_func_pl( i, T_L, T_ref, R_gas, bi, pl ):
   
   # Activation moderated by deactivation
   
   # Activation 
   ####################
   DeltaT = T_L - T_ref
   
   Eff_activ  = float( pl.H_a[i] * DeltaT )
   norm_activ = float( T_ref * T_L * R_gas )
   
   Activation = Tdep_KK_funcExpFuncs( Eff_activ, norm_activ )
   
   # Deactivation: ratio b/n ref & leaf
   ####################
   
   # _ref Temperature
   Entropy_ref =  ( T_ref  * pl.DeltaS[i] ) - pl.H_d[i]
   norm_ref = R_gas * T_ref
   Deactiv_ref = 1. + Tdep_KK_funcExpFuncs( Entropy_ref, norm_ref )
   
   # leaf Temperature
   Entropy_leaf = ( T_L * pl.DeltaS[i] ) - pl.H_d[i]
   norm_leaf = R_gas * T_L
   Deactiv_leaf = 1. + Tdep_KK_funcExpFuncs( Entropy_leaf, norm_leaf )
  
   fn_T_L =  Activation * Deactiv_ref / Deactiv_leaf
   
   return fn_T_L


###############################################################################

def Tdep_KK_funcExpFuncs( numer, denom ):
   term = np.exp( numer/denom )
   return term 

###############################################################################


