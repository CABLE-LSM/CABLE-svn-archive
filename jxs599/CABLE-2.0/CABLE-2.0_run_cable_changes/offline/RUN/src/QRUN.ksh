
#PBS -l walltime=800
#PBS -l vmem=3072MB
#PBS -l ncpus=1
#PBS -l jobfs=2GB
#PBS -q express 
#PBS -j oe
#PBS -p 10
#PBS -N CABLE 

#if qsub-ed go to the dir where job was launched from 
cd $PBS_O_WORKDIR

RUN_CABLE run plot
