import cdms2, sys
import numpy as np
import os

cdms2.setNetcdfShuffleFlag(0)      ## where value is either 0 or 1
cdms2.setNetcdfDeflateFlag(0)      ## where value is either 0 or 1
cdms2.setNetcdfDeflateLevelFlag(0) ## where value is a integer bw 0 and 9 inc

cfile=cdms2.open('MeanMnthDailyCycles_%syrs.nc' % sys.argv[1])
tas = cfile['temp']       # tas
qsh = cfile['sh']         # hfss
qle = cfile['lh']         # hfls
apr = cfile['tot_precip'] # field5226
rnt = cfile['field202']   # field3333
#swr = cfile['solar']
#lwr = cfile['longwave']

tstep = int(os.getenv('TSTEP'))
#pals_shift_hrs = [,,,,,2,,,10,1,-6 ,,,,1,-5 ,-6 ,]
#pals_shift     = [,,,,,4,,,20,2,-12,,,,2,-10,-12,]

#real_shift_hrs = [-4,10,9.5,1,5.5,2,1,1,10,1,-6,8,-6,-5,1,-5,-6,1]
#tshift_new = [-8,20,19,2,11,4,2,2,20,2,-12,16,-12,-10,2,-10,-12,2]
#tshift_old = [-7,20,19,2,11,4,2,2,20,2,-9,15,-10,-8,2,-7,-10,2]

if tstep == 48:
    tshift = [-7,20,19,2,11,4,2,2,20,2,-9,15,-10,-8,2,-7,-10,2]
if tstep == 72:
    tshift = [-7*1.5,20*1.5,19*1.5,2*1.5,11*1.5,4*1.5,2*1.5,2*1.5,20*1.5,2*1.5,-9*1.5,15*1.5,-10*1.5,-8*1.5,2*1.5,-7*1.5,-10*1.5,2*1.5]
if tstep == 288:
    tshift = [-7*6,20*6,19*6,2*6,11*6,4*6,2*6,2*6,20*6,2*6,-9*6,15*6,-10*6,-8*6,2*6,-7*6,-10*6,2*6]
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

