#!/usr/bin/python
import os

###############################################################################
TestUM = True 
TestMPI = True
TestSerial = True
###############################################################################

# define dirs
user = 'jxs599'
institute = '599'  # CSIRO
project = 'p66'    # ACCESS-CSIRO-general project

# branch from the CABLE repository branch to test
TicketNum = "52"
SVNURL = "https://trac.nci.org.au/svn/cable/branches/Share/Tickets2015/" + TicketNum    
# i.e. the root of your CABLE code containing UM,core, offline sub-directories
SVNURLROOT = "/" + TicketNum

# root of TestSuite Directory
root = os.getcwd()
uHome = '/home/' + institute + '/'+ user 
# root of TestSuiteApps Directory, etc,etc
root_app = root + '/TestSuiteApps'
src = root_app + '/src' 
bin = root_app + '/bin'
run = root_app + '/Run'
cfgs = root + '/TestSuiteConfigs/'
trunk = src + SVNURLROOT
UM = trunk + '/UM'
offline = trunk + '/offline'

#Store for built executables

# Single Sites
binSerial = bin + '/Serial'

# MPI-parallel - generally gswp2
binParallel= bin + '/Parallel'

# UM Build and Run dirs - this is config dependent on the particular UM run to be tested
# jaaaD is ~ACCESS-1.4 with carbon-cycle switched off 
#UMcfgs = '/g/data1/p66/CABLE/TestSuiteConfigs/JAAAF'
UMsrc = 'vaheb'
UMrun = UMsrc + '-*' 
UMcfgs = cfgs + UMsrc 

###############################################################################





