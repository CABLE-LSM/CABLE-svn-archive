import os
import subprocess

def driver( fields,
            pars,
            config,
            flags ):

    dir_cwd = os.getcwd()
    execu2 = dir_cwd + '/' + pars.path[2]
    execu3 = dir_cwd + '/' + pars.path[3]
    dir_map = dir_cwd + '/' + pars.path[0]
    dir_cat = dir_cwd + '/' + pars.path[1]
    dir_ncdf = dir_cwd + '/' + pars.path[4]
    
    # IF mapping is to be processed
    if config.map[len(config.map)-1] is True:

        os.chdir(dir_map)

        for mm in range( len(pars.mapfield) ):
            #execute system commands
            catcmd = ( execu2 + " " + pars.mapfield[mm] + " " + config.nodes[0] )

            print '\nExecuting...\n', catcmd
            #subprocess.call(catcmd)
            p2 = subprocess.Popen(catcmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = p2.communicate()
            print 'Done.'

            #subprocess.call(cpcmd)
            cpcmd = ( "/bin/cp " + pars.mapfield[mm] + ".* " + dir_cat )
            print '\nCopying...\n', cpcmd
            cpcat = subprocess.Popen(cpcmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = cpcat.communicate()

    # IF there are fields to be processed
    if( len(fields.name) > 0 ):
        
        for ff in range( len(fields.name) ):
        
            dir_var = dir_cwd + '/' + fields.path[ff]
            os.chdir(dir_var)
            #print '\n\npwd: \n\n ', os.getcwd()
            funit = open("input.dat",'w')
            funit.write(fields.name[ff]+'\n')
            funit.write(fields.name[ff] +'.nc\n')
            funit.close()

            #execute system commands
            catcmd = ( execu2 + " " + fields.name[ff] + " " + config.nodes[0] )

            print '\nExecuting...\n', catcmd
            #subprocess.call(catcmd)
            p2 = subprocess.Popen(catcmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = p2.communicate()
            print 'Done.'
  
            #subprocess.call(cpcmd)
            cpcmd = ( "/bin/cp " + fields.name[ff] + ".* " + dir_cat )
            print '\nCopying...\n', cpcmd
            cpcat = subprocess.Popen(cpcmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = cpcat.communicate()
            
            # call fortran executable to generate mapped netcdf file
            nccmd = ( execu3 + " " + dir_cat + "/" )

            print '\nExecuting...\n', nccmd
            #subprocess.call(nccmd)
            p3 = subprocess.Popen(nccmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = p3.communicate()
            print 'Done.'
            
            #   system("mv $basename.nc $dir_ncdf" );
            #subprocess.call(mvcmd)
            mvcmd = ( "/bin/mv " + fields.name[ff] + ".nc " + dir_ncdf )
            print '\nCopying...\n', mvcmd
            cpnc = subprocess.Popen(mvcmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = cpnc.communicate()


