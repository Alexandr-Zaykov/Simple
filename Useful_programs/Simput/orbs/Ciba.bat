#!/bin/bash

##SBATCH --job-name
#SBATCH --nodes 1
#SBATCH --output job.out
#SBATCH --ntasks-per-node 4
#SBATCH --qos=preemptable
##SBATCH --partition="bnode0112,bnode0113"

###ENVIRONMENT###

module load slurm/blanca

module load gcc

module unload openmpi

export NBOEXE=/projects/$LOGNAME/quantum/nbo6/bin/nbo6.i4.exe

export PATH=/projects/$LOGNAME/local/openmpi-2.1.1/bin:/projects/$LOGNAME/quantum/orca_4_0_0_2_linux_x86-64:${PATH}

export LD_LIBRARY_PATH=/projects/$LOGNAME/local/openmpi-2.1.1/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

###JOB INPUT###

export JOB=$1
export WorkDir=/rc_scratch/$LOGNAME/$JOB
mkdir -p $WorkDir
export HomeDir=$PWD
cp -rp $HomeDir/* $WorkDir
cd $WorkDir

###JOB START###

/projects/$LOGNAME/quantum/orca_4_0_0_2_linux_x86-64/orca "$JOB".inp > "$HomeDir"/"$JOB".out 2> "$HomeDir"/"$JOB".err

###JOB END###

rm $WorkDir/*.tmp.* > /dev/null 2>&1
rm $WorkDir/*.tmp > /dev/null 2>&1
rm $WorkDir/core.* > /dev/null 2>&1
#rm $WorkDir/FILE.* > /dev/null 2>&1
cp $WorkDir/* $HomeDir/
cd $HomeDir
#rm -r $WorkDir

