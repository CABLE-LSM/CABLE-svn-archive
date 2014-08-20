import numpy as np
import sys 

def Set_plants( nlines, nBiomes, pl, lines ): 
   
   # Loop over all of these read lines
   ll = 0
   for l in range( nlines ):

      # Create a string from each line
      lstr = str(lines[l])
      # nullify case sensitivity in that string
      lstr = lstr.lower()

      # if it is meant to be interpreted - interpret based on what
      # it starts with. See configfile example format
      # check if the string is a comment
      cfield = lines[l].strip().split()
      if not lstr.startswith('#'):
         if int(cfield[5]) > -1:
            pl.PlNumber[ll] = (cfield[0])
            pl.BiNumber[ll] = (cfield[1])
            pl.Vcmax_25[ll] = (cfield[4])
            pl.H_a[ll] = (cfield[5])
            pl.H_d[ll] = (200000.)
            pl.DeltaS[ll] = (cfield[6])
            ll = ll + 1
         #pl.Jmax_25.append(cfield[0])
   
   nBiomes.append( int( pl.BiNumber[ll-1] ) )
   return ll 
         ## Mode
         #if lstr.startswith( "mode" ):
         #   cfield = lines[ll].strip().split()
         #   cfg.mode.append(cfield[1])

         ## Name
         #if lstr.startswith( "name" ):
         #   cfield = lines[ll].strip().split()
         #   cstr = ""
         #   for istr in range( 1,len(cfield) ):
         #      cstr = cstr + cfield[istr]
         #   cfg.name.append(cstr)

         ## Path
         #if lstr.startswith( "path" ):
         #   cfield = lines[ll].strip().split()
         #   cfg.path.append(cfield[1])

      # End IF not comment line

   # End Loop over read lines[]





