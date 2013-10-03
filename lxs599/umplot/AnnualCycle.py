import cdms2, MV2, os, sys

cdms2.setNetcdfShuffleFlag(0) ## where value is either 0 or 1
cdms2.setNetcdfDeflateFlag(0) ## where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(0) ## where value is a integer bw 0 and 9 inc

year  = os.getenv('YR')
cfile = cdms2.open('Timeseries_%syrs.nc' % sys.argv[1])
ta    = cfile['tas']
qh    = cfile['hfss']
qle   = cfile['hfls']
apr   = cfile['field5226']
rnt   = cfile['field3333']

tstep = (ta.shape[0]/(12*30*int(year)))

tas = MV2.reshape(ta,(int(year),12,30,tstep,ta.shape[1]))
sh  = MV2.reshape(qh,(int(year),12,30,tstep,qh.shape[1]))
lh  = MV2.reshape(qle,(int(year),12,30,tstep,qle.shape[1]))
ap  = MV2.reshape(apr,(int(year),12,30,tstep,apr.shape[1]))
rn  = MV2.reshape(rnt,(int(year),12,30,tstep,rnt.shape[1]))
T1  = MV2.average(tas[:,:,:,:,:],3)
L1  = MV2.average(lh[:,:,:,:,:],3)
S1  = MV2.average(sh[:,:,:,:,:],3)
P1  = MV2.average(ap[:,:,:,:,:],3)
R1  = MV2.average(rn[:,:,:,:,:],3)
T2  = MV2.average(T1[:,:,:,:],2)
L2  = MV2.average(L1[:,:,:,:],2)
S2  = MV2.average(S1[:,:,:,:],2)
P2  = MV2.average(P1[:,:,:,:],2)
R2  = MV2.average(R1[:,:,:,:],2)
T3  = MV2.average(T2[:,:,:],0)
L3  = MV2.average(L2[:,:,:],0)
S3  = MV2.average(S2[:,:,:],0)
P3  = MV2.average(P2[:,:,:],0)
R3  = MV2.average(R2[:,:,:],0)

cfout = cdms2.createDataset('AnnCycle_%syrs.nc' % sys.argv[1])
cfout.write(T3,id=ta.id)
cfout.write(L3,id=qle.id)
cfout.write(S3,id=qh.id)
cfout.write(P3,id=apr.id)
cfout.write(R3,id=rnt.id)
cfout.write(T2,id='tair_%syrs' % sys.argv[1])
cfout.write(L2,id='hfls_%syrs' % sys.argv[1])
cfout.write(S2,id='hfss_%syrs' % sys.argv[1])
cfout.write(P2,id='ppt_%syrs' % sys.argv[1])
cfout.write(R2,id='rnet_%syrs' % sys.argv[1])
cfout.sync()

cfile.close()
