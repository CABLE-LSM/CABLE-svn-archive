import cdms2, MV2, cdutil, cdtime, calendar
import numpy as np
import os, sys

cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(0)
cdms2.setNetcdfDeflateLevelFlag(0)

#set shwtm=`cdo showtime Timeseries_27yrs.nc`
#echo $shwtm[$i+48]

year = os.getenv('YR')
tfile = cdms2.open('Timeseries_%syrs.nc' % sys.argv[1])
#t2file = cdms2.open('Timeseries_%syrs_noswlw.nc' % sys.argv[1])
#t3file = cdms2.open('Timeseries_%syrs_swlw.nc' % sys.argv[1])
ystrt = int(sys.argv[2])
yend = int(sys.argv[3])

tas = tfile['tas']
qh  = tfile['hfss']
qle = tfile['hfls']
apr = tfile['field5226']
rnt = tfile['field3333']
#tas = t2file['temp']
#qh  = t2file['lh']
#qle = t2file['sh']
##apr = t2file['precip']
#apr = t2file['tot_precip']
#rnt = t2file['field202']
#u = t2file['u']
#v = t2file['v']
#mflx = t2file['field184']
#u1 = t2file['u_1']
#v1 = t2file['v_1']
#wnd = t2file['wind']
#sevp = t2file['field1526']
#cevp = t2file['field1527']
#ievp = t2file['field1528']
#tmht = t2file['field1534']
#friv = t2file['field1696']
#tsl = t2file['soiltemp']
#clt = t2file['field30']
#ts = t2file['temp_1']
#u2 = t2file['u_2']
#v2 = t2file['v_2']
#sndep = t2file['snowdepth']
#ts2 = t2file['temp_2']
#blht = t2file['blht']
#sw = t3file['solar']
#swd = t3file['field203'] 
#lw = t3file['longwave'] 
#lwd = t3file['ilr'] 

tstep = 48 # hardwired, (tas.shape[0]/(12*30*int(year)))
#tstep = (tas.shape[0]/(365*int(year)))
#tstep = (tas.shape[0]/(365.25*int(year)))
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

            # Mn = L[mnth-1]
            #else
            # Mn = l[mnth-1]
            #print(yr,mnth,L[mnth-1],tas.shape[1])    

            #t = MV2.average(tas(time=(cdtime.comptime(yr,mnth,1,0,30,0),cdtime.comptime(yr,m+1,1,0,0,0))).reshape(L[mnth-1],48,1),0)
            t = MV2.average(tas(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],48,tas.shape[1]),0) 
            a[y, mnth-1,:,:] = t[:,:] 

            q = MV2.average(qh(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],48,tas.shape[1]),0)
            b[y, mnth-1,:,:] = q[:,:] 

            s = MV2.average(qle(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],48,tas.shape[1]),0)
            c[y, mnth-1,:,:] = s[:,:] 

            h = MV2.average(apr(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],48,tas.shape[1]),0)
            d[y, mnth-1,:,:] = h[:,:]

            j = MV2.average(rnt(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,L[mnth-1],23,59,59))).reshape(L[mnth-1],48,tas.shape[1]),0)
            i[y, mnth-1,:,:] = j[:,:]

        else:
            #print(yr,mnth,l[mnth-1],tas.shape[1])    
            t = MV2.average(tas(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],48,tas.shape[1]),0)
            a[y, mnth-1,:,:] = t[:,:]

            q = MV2.average(qh(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],48,tas.shape[1]),0)
            b[y, mnth-1,:,:] = q[:,:]
            
            s = MV2.average(qle(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],48,tas.shape[1]),0)
            c[y, mnth-1,:,:] = s[:,:]

            h = MV2.average(apr(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],48,tas.shape[1]),0)
            d[y, mnth-1,:,:] = h[:,:]

            j = MV2.average(rnt(time=(cdtime.comptime(yr,mnth,1,0,0,0),cdtime.comptime(yr,mnth,l[mnth-1],23,59,59))).reshape(l[mnth-1],48,tas.shape[1]),0)
            i[y, mnth-1,:,:] = j[:,:]

    y = y + 1


A = MV2.average(a,0)
B = MV2.average(b,0)
C = MV2.average(c,0)
D = MV2.average(d,0)
I = MV2.average(i,0)

cfout = cdms2.createDataset('MMDC_%syrs.nc' % sys.argv[1])
cfout.write(A,id=tas.id)
cfout.write(B,id=qh.id)
cfout.write(C,id=qle.id)
cfout.write(D,id=apr.id)
cfout.write(I,id=rnt.id)
cfout.write(a)#,id='tas_yrmean')
cfout.write(b)#,id='qh_yrmean')
cfout.write(c)#,id='qle_yrmean')
cfout.write(d)#,id='apr_yrmean')
cfout.write(i)#,id='rnt_yrmean')
cfout.sync()

tfile.close()
