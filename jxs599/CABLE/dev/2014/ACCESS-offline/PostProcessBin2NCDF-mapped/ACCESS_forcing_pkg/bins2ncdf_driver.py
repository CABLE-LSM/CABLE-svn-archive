import os
import subprocess

def driver( fields,
            pars,
            config,
            flags ):

    dir_cwd = os.getcwd()
    execu = dir_cwd + '/' + pars.path[2]
    dir_map = dir_cwd + '/' + pars.path[0]
    
    if config.map[0] is True:

        os.chdir(dir_map)

        for mm in range( len(pars.mapfield) ):
            #execute system commands
            cmd = ( execu + " " + pars.mapfield[mm] + " " + config.nodes[0] )
            print cmd

    #subprocess.call(cmd)
    #p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    #(output, err) = p.communicate()

    #for ff in range( len(fields.name) ):
    #    print fields.name[ff]
