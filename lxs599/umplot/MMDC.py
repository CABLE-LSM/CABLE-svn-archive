import cdms2, MV2, os, sys

#Sflag = cdms2.getNetcdfShuffleFlag()
#Dflag = cdms2.getNetcdfDeflateFlag()
#Lflag = cdms2.getNetcdfDeflateLevelFlag()
#print Sflag,Dflag,Lflag
#if hasattr(cdms2, 'setNetcdfDeflateFlag'):
cdms2.setNetcdfShuffleFlag(0) ## where value is either 0 or 1
cdms2.setNetcdfDeflateFlag(0) ## where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(0) ## where value is a integer bw 0 and 9 inc
#0 1 1
#/apps/python/2.7.2/lib/python2.7/site-packages/cdms2/dataset.py:1660: Warning: Since CDAT Version 5.2 File are now written with compression and shuffling
#You can query different values of compression using the functions:
#cdms2.getNetcdfShuffleFlag() returning 1 if shuffling is enabled, 0 otherwise
#cdms2.getNetcdfDeflateFlag() returning 1 if deflate is used, 0 otherwise
#cdms2.getNetcdfDeflateLevelFlag() returning the level of compression for the deflate method
#If you want to turn that off or set different values of compression use the functions:
#cdms2.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
#cdms2.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
#cdms2.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included
#Turning all values to 0 will produce NetCDF3 Classic files


year  = os.getenv('YR')
cfile = cdms2.open('Timeseries_%syrs.nc' % sys.argv[1])
ta    = cfile['tas']
qh    = cfile['hfss']
qle   = cfile['hfls']
apr   = cfile['field5226']
rnt   = cfile['field3333']

tstep = (ta.shape[0]/(12*30*int(year)))

tas   = MV2.reshape(ta,(int(year),12,30,tstep,ta.shape[1]))
sh    = MV2.reshape(qh,(int(year),12,30,tstep,qh.shape[1]))
lh    = MV2.reshape(qle,(int(year),12,30,tstep,qle.shape[1]))
ap    = MV2.reshape(apr,(int(year),12,30,tstep,apr.shape[1]))
rn    = MV2.reshape(rnt,(int(year),12,30,tstep,rnt.shape[1]))
T1    = MV2.average(tas[:,:,:,:,:],2)
L1    = MV2.average(lh[:,:,:,:,:],2)
S1    = MV2.average(sh[:,:,:,:,:],2)
P1    = MV2.average(ap[:,:,:,:,:],2)
R1    = MV2.average(rn[:,:,:,:,:],2)
T2    = MV2.average(T1[:,:,:,:],0)
L2    = MV2.average(L1[:,:,:,:],0)
S2    = MV2.average(S1[:,:,:,:],0)
P2    = MV2.average(P1[:,:,:,:],0)
R2    = MV2.average(R1[:,:,:,:],0)

cfout = cdms2.createDataset('MMDC_%syrs.nc' % sys.argv[1])
cfout.write(T2,id=ta.id)
cfout.write(L2,id=qle.id)
cfout.write(S2,id=qh.id)
cfout.write(P2,id=apr.id)
cfout.write(R2,id=rnt.id)
cfout.write(T1)
cfout.write(L1)
cfout.write(S1)
cfout.write(P1)
cfout.write(R1)
cfout.sync()

cfile.close()
