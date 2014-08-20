import numpy as np

def Read_dataset( nlines, ifile ):

   # Declare mutable object => lines read from config file
   #lines = []

   # Reads all the lines in the config file specified at CLI
   lines = tuple( open(ifile, 'r' ))

   nlines.append( int( len(lines) ) )

   return lines
###############################################################################

