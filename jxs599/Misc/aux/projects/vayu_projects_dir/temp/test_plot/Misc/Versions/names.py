#!/usr/bin/python

import os
filename = "Output/namelist.dat"
namelist = ('myGamma.001\n', 'myGamma.005\n','myGamma.010\n','myGamma.040\n','myGamma.060\n','myGamma.080\n','myGamma.100\n') 
for value in namelist:
#xvalue="Output/myGamma.001\n"
	FILE = open(filename,"w")
	FILE.write(value)
	FILE.close()
	os.system("./match")
   
#    FILE.writelines(namelist)

