import cdms2, MV2, cdutil, cdtime, calendar
import numpy as np
import os, sys

cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(0)
cdms2.setNetcdfDeflateLevelFlag(0)

year = os.getenv('YR')
tfile = cdms2.open('Timeseries_%syrs.nc' % sys.argv[1])
ystrt = int(sys.argv[2])
yend = int(sys.argv[3])

tas = tfile['temp']       # tas
qh  = tfile['sh']         # hfss
qle = tfile['lh']         # hfls
apr = tfile['tot_precip'] # field5226
rnt = tfile['field202']   # field3333

tstep = int(os.getenv('TSTEP'))
#tstep = 48 # hardwired, (tas.shape[0]/(12*30*int(year)))
nsites = tas.shape[1]

a = MV2.zeros((int(year),12,tstep,tas.shape[1]))
b = MV2.zeros((int(year),12,tstep,tas.shape[1]))
c = MV2.zeros((int(year),12,tstep,tas.shape[1]))
d = MV2.zeros((int(year),12,tstep,tas.shape[1]))
i = MV2.zeros((int(year),12,tstep,tas.shape[1]))
l = [31,28,31,30,31,30,31,31,30,31,30,31] 
L = [31,29,31,30,31,30,31,31,30,31,30,31]
y = 0

#if ystrt==yend:
#    yend=yend+1

for yr in range(ystrt,yend+1):
 
    for mnth in range(1,13):

        #year = yr
        #m = mnth
        #if mnth == 12:
        #    m = 0
        #    year = yr+1

        if calendar.isleap(yr):

            #t = MV2.average(tas(time=(cdtime.comptime(yr,mnth,1,0,30,0),cdtime.comptime(yr,m+1,1,0,0,0))).reshape(L[mnth-1],tstep,1),0) 
            t = MV2.average(tas(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],tstep,tas.shape[1]),0) 
            a[y, mnth-1,:,:] = t[:,:] 

            q = MV2.average(qh(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],tstep,tas.shape[1]),0)
            b[y, mnth-1,:,:] = q[:,:] 

            s = MV2.average(qle(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],tstep,tas.shape[1]),0)
            c[y, mnth-1,:,:] = s[:,:] 

            h = MV2.average(apr(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],tstep,tas.shape[1]),0)
            d[y, mnth-1,:,:] = h[:,:]

            j = MV2.average(rnt(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],tstep,tas.shape[1]),0)
            i[y, mnth-1,:,:] = j[:,:]

        else:
            
            t = MV2.average(tas(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],tstep,tas.shape[1]),0)
            a[y, mnth-1,:,:] = t[:,:]

            q = MV2.average(qh(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],tstep,tas.shape[1]),0)
            b[y, mnth-1,:,:] = q[:,:]
            
            s = MV2.average(qle(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],tstep,tas.shape[1]),0)
            c[y, mnth-1,:,:] = s[:,:]

            h = MV2.average(apr(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],tstep,tas.shape[1]),0)
            d[y, mnth-1,:,:] = h[:,:]

            j = MV2.average(rnt(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],tstep,tas.shape[1]),0)
            i[y, mnth-1,:,:] = j[:,:]

    y = y + 1


A  = MV2.average(a,2)
B  = MV2.average(b,2)
C  = MV2.average(c,2)
D  = MV2.average(d,2)
I  = MV2.average(i,2)
aa = MV2.average(A,0)
bb = MV2.average(B,0) 
cc = MV2.average(C,0) 
dd = MV2.average(D,0) 
ii = MV2.average(I,0) 

cfout = cdms2.createDataset('AnnCycle_%syrs.nc' % sys.argv[1])
cfout.write(aa,id=tas.id)
cfout.write(bb,id=qh.id)
cfout.write(cc,id=qle.id)
cfout.write(dd,id=apr.id)
cfout.write(ii,id=rnt.id)
cfout.write(A,id='tas_yrmean')
cfout.write(B,id='qh_yrmean')
cfout.write(C,id='qle_yrmean')
cfout.write(D,id='apr_yrmean')
cfout.write(I,id='rnt_yrmean')
cfout.sync()

tfile.close()
