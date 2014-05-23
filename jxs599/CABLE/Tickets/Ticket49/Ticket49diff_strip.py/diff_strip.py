#!/usr/bin/python

#r1755 | jxs599 | 2014-04-29 21:59:22 +1000 (Tue, 29 Apr 2014) | 1 line
#export SLI over the top and ci changes to existing files

#to use system calls
#import os

#import numpy as N
#ddir = N.array( [ [0,lrdir],[0,lftot] ])
#print ddir.shape

# root of branches to diff
SLIdir = 'Ticket49_tr287_SLI'
trdir =  'Ticket49_tr287'

# files to diff are across model dirs
chemdir = 'core/biogeochem'
physdir = 'core/biogeophys'
offline = 'offline'

#files to process in the biogeochem dir
fchem = [                              
          'casa_cable.F90',             
          'casa_cnp.F90',                
          'casa_inout.F90',              
          'casa_variable.F90'            
        ]   

#files to process in the biogeophys dir
fphys = [                              
          'cable_air.F90',
          'cable_albedo.F90',
          'cable_canopy.F90',
          'cable_carbon.F90',
          'cable_cbm.F90',
          'cable_common.F90',
          'cable_data.F90',
          'cable_define_types.F90',
          'cable_radiation.F90',
          'cable_roughness.F90',
          'cable_soilsnow.F90',
        ]

#files to process in the offline dir
foffl = [                              
          'Makefile_offline',
          'build.ksh',
          'cable_abort.F90',
          'cable_checks.F90',
          'cable_driver.F90',
          'cable_initialise.F90',
          'cable_input.F90',
          'cable_iovars.F90',
          'cable_output.F90',
          'cable_parameters.F90',
          'cable_read.F90',
          'cable_write.F90'
        ] 

rootdir = [
            SLIdir,  
            trdir
          ]
            
modeldir = [
             chemdir,
             physdir,
             offline
           ]     
# lengths for future convenience
lrdir = len(rootdir)
lmdir = len(modeldir)
lfch = len(fchem)
lfph = len(fphys)
lfof = len(foffl)
lftot = lfch + lfph + lfof

ddir = [[ 0 for x in xrange(lftot)] for x in xrange(lrdir)] 

i=0
for rdir in rootdir:
    j=0
    for mdir in modeldir: 
        if mdir == chemdir: 
            k=0
            for fch in fchem:
                k=k+1
                dirq = rdir + '/' + mdir + '/' + fch
                #print dirq     
                ddir[i][j] = dirq
                j=j+1
                #print mdir
                #print fch
        
        if mdir == physdir: 
            for fph in fphys:
                dirq = rdir + '/' + mdir + '/' + fph 
                ddir[i][j] = dirq
                j=j+1
    
        if mdir == offline: 
            for fof in foffl:
                dirq = rdir + '/' + mdir + '/' + fof 
                ddir[i][j] = dirq
                j=j+1
    i=i+1

for i in range(lrdir):
    for j in range(lftot):
        print ddir[i][j]

#cmd = 'diff ' + ddir[1][1] + ' ' + dir[2][1] + ' > junk'
#cmd = 'touch junk'
#os.system(cmd)

import subprocess

#execut system comand
cmd = "/usr/bin/diff ../tralb/cable_canopy.F90 ../SLIalb/cable_canopy.F90 -i -w --strip-trailing-cr"
#print cmd

#subprocess.call(cmd)
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
(output, err) = p.communicate()
print "Today is", output
#p.wait()


#
##	diff tralb/cable_albedo.F90 SLIalb//cable_albedo.F90 -i -w --strip-trailing-cr > junk
#
##python testing
##hw = "Hello World!"
##print hw
##for i in range( len(fchem) ): 
