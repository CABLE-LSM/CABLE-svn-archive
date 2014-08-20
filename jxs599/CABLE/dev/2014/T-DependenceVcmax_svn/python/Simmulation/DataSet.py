import numpy as np
import struct
import sys
from main_data import CABLEraw, CABLEdiurnal, recl

def ReadBinary( ifile ):
   
   # open the file for reading,binary=rb
   f = open(ifile, 'rb')
   
   # read everytnig in the file
   raw = f.read( )
   
   # how long is the read data in bytes
   recl = len(raw) // 12
   print recl/3 

   # unpack binary data as a string
   # Fortran formats bin data as recl,rec,recl
   # where recl = record length. therefore use format structure "ifi"
   dat = struct.unpack("ifi" * recl , raw)
   
   j=0 
   # fill array with real data (NOT "recl" bytes)
   for i in range (0,(len(dat) -1),3):
      CABLEraw[j] = dat[i+1]
      j=j+1
   
   f.close() 
   
   return CABLEraw
###############################################################################

