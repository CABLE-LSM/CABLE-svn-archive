import cdms2, sys
import numpy as np

cdms2.setNetcdfShuffleFlag(0)      ## where value is either 0 or 1
cdms2.setNetcdfDeflateFlag(0)      ## where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(0) ## where value is a integer bw 0 and 9 inc

cfile=cdms2.open('MeanMnthDailyCycles_%syrs.nc' % sys.argv[1])
tas = cfile['tas']
qsh = cfile['hfss']
qle = cfile['hfls']
apr = cfile['field5226']
rnt = cfile['field3333']
#swr = cfile['solar']
#lwr = cfile['longwave']

tshift = [-7,20,19,2,11,4,2,2,20,2,-9,15,-10,-8,2,-7,-10,2]
sites  = ['manaus','hay','dalywate','hapexbat','india','hyytiala','africaN','africaS','tumbarum','tharandt','bondvill','dinghush','litwashi','walkerbr','loobos','harvard','nsabor','vielsalm']
nsite  = tshift.__len__()
nsite2 = sites.__len__()

if nsite2 != nsite:
    print("ERROR in number of sites",nsite,nsite2)

# ------------------------------------------------------------------------------

for ns in range(0,nsite):
    #print(ns,tshift[ns],sites[ns])

    tas_rl=np.roll(tas[:,:,ns],tshift[ns],axis=1)
    qsh_rl=np.roll(qsh[:,:,ns],tshift[ns],axis=1)
    qle_rl=np.roll(qle[:,:,ns],tshift[ns],axis=1)
    apr_rl=np.roll(apr[:,:,ns],tshift[ns],axis=1)
    rnt_rl=np.roll(rnt[:,:,ns],tshift[ns],axis=1)
    #swr_rl=swr[:,:,ns] #np.roll(swr[:,:,ns],tshift[ns],axis=1)
    #lwr_rl=lwr[:,:,ns] #np.roll(lwr[:,:,ns],tshift[ns],axis=1)

    cfout=cdms2.createDataset('mmdc_%(first)s_roll_%(last)syrs.nc' % {'first': sites[ns], 'last': sys.argv[1]})
    cfout.write(tas_rl,id=tas.id)
    cfout.write(qsh_rl,id=qsh.id)
    cfout.write(qle_rl,id=qle.id)
    cfout.write(apr_rl,id=apr.id)
    cfout.write(rnt_rl,id=rnt.id)
    #cfout.write(swr_rl,id=swr.id)
    #cfout.write(lwr_rl,id=lwr.id)
    cfout.sync()

cfile.close()

