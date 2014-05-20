import cdms2, sys
import numpy as np

cfile = cdms2.open('MMDC_%syrs.nc' % sys.argv[1])

# For a start date other than Jan

tas   = cfile['tas']     # (18,TS,12)
lh    = cfile['hfls']
sh    = cfile['hfss']
pr    = cfile['field5226']
rnet  = cfile['field3333']
tas2  = cfile['variable_14'] # (21,TS,12,20)
lh2   = cfile['variable_18']
sh2   = cfile['variable_22']
pr2   = cfile['variable_26']
rnet2 = cfile['variable_30']

tas3  = np.roll(tas  ,-3,axis=0)
lh3   = np.roll(lh   ,-3,axis=0)
sh3   = np.roll(sh   ,-3,axis=0)
pr3   = np.roll(pr   ,-3,axis=0)
rnet3 = np.roll(rnet ,-3,axis=0)
tas4  = np.roll(tas2 ,-3,axis=1)
lh4   = np.roll(lh2  ,-3,axis=1)
sh4   = np.roll(sh2  ,-3,axis=1)
pr4   = np.roll(pr2  ,-3,axis=1)
rnet4 = np.roll(rnet2,-3,axis=1)

cfout = cdms2.createDataset('MeanMnthDailyCycles_%syrs.nc' % sys.argv[1])
cfout.write(tas3 ,id=tas.id ,axes=[tas.getAxisList()[0] ,tas.getAxisList()[1] ,tas.getAxisList()[2]])
cfout.write(lh3  ,id=lh.id  ,axes=[lh.getAxisList()[0]  ,lh.getAxisList()[1]  ,lh.getAxisList()[2]])
cfout.write(sh3  ,id=sh.id  ,axes=[sh.getAxisList()[0]  ,sh.getAxisList()[1]  ,sh.getAxisList()[2]])
cfout.write(pr3  ,id=pr.id  ,axes=[pr.getAxisList()[0]  ,pr.getAxisList()[1]  ,pr.getAxisList()[2]])
cfout.write(rnet3,id=rnet.id,axes=[rnet.getAxisList()[0],rnet.getAxisList()[1],rnet.getAxisList()[2]])

cfout.write(tas4 ,id=tas2.id ,axes=[tas2.getAxisList()[0] ,tas2.getAxisList()[1] ,tas2.getAxisList()[2] ,tas2.getAxisList()[3]])
cfout.write(lh4  ,id=lh2.id  ,axes=[lh2.getAxisList()[0]  ,lh2.getAxisList()[1]  ,lh2.getAxisList()[2]  ,lh2.getAxisList()[3]])
cfout.write(sh4  ,id=sh2.id  ,axes=[sh2.getAxisList()[0]  ,sh2.getAxisList()[1]  ,sh2.getAxisList()[2]  ,sh2.getAxisList()[3]])
cfout.write(pr4  ,id=pr2.id  ,axes=[pr2.getAxisList()[0]  ,pr2.getAxisList()[1]  ,pr2.getAxisList()[2]  ,pr2.getAxisList()[3]])
cfout.write(rnet4,id=rnet2.id,axes=[rnet2.getAxisList()[0],rnet2.getAxisList()[1],rnet2.getAxisList()[2],rnet2.getAxisList()[3]])
cfout.sync()

cfile.close()

