from ctypes import *

mylib = CDLL('cable_python_driver.so')
import numpy as np
nparams=c_int(18)
ntimesteps=c_int(87648)
nfluxes=c_int(2)

param_vec = c_double * 18

parms=param_vec(7e-4,5.0,2.0,1.5,17.0, 0.01 ,   0.055,0.062,0.302,0.010,2.0,0.943,20.0,0.00004,3.0,2.0,1.0,0.1)


time_vec = c_double * 87648



#(7e-4,5.0,2.0,1.5,17.0, 0.01 ,   0.055,0.062,0.302,0.010,2.0,0.943,20.0,0.00004,3.0,2.0,1.0,0.1)

#parms[0:3] = (7e-4,5.0,2.0,1.5)
#parms[4:]  = ( 17.0, 0.01 ,   0.055,0.062,0.302,0.010,2.0,0.943,20.0,0.00004,3.0,2.0,1.0,0.1)

#byref(parms)


#QHfluxes=time_vec(0.0)
#QLEfluxes=time_vec(0.0)

QHfluxes=time_vec()
QLEfluxes=time_vec()

pntr_params=pointer(parms)
pntr_qh = pointer(QHfluxes)
pntr_qle = pointer(QLEfluxes)

mylib.c_cable.restype = None
#mylib.c_cable(nparams,ntimesteps,nfluxes,pntr_params,pntr_qh,pntr_qle)
mylib.c_cable(byref(nparams),byref(ntimesteps),byref(nfluxes),byref(parms),byref(QHfluxes),byref(QLEfluxes))

