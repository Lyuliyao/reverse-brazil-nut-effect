#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
#casePath="$PWD"
echo "casePath: $casePath"
logpath=$casePath
headerText="run_parallel_cfdemSolverPiso_periodicChannel_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs="${CFDEM_NPROCS:-8}"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
separateDEM="false"
reconstuctCase="true"
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="false"
postproc="false"
cleanCase="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode $separateDEM $reconstuctCase
python3 $casePath/check_force_balance.py "$casePath" --latest --balance-tol 1e-8 --warn-only
forcecheck_status=$?
if [ $forcecheck_status -ne 0 ]
    then
        exit $forcecheck_status
fi


if [ $runOctave == "true" ]
    then
        #------------------------------#
        #  octave

        #- change path
        cd octave

        #- rmove old graph
        rm *.png

        #- run octave
        matlab --no-gui checkVolFlow.m

        #- show plot 
        eog volflow.png

        #- copy log file to test harness
        cp ../../$logfileName $testHarnessPath
        cp *.png $testHarnessPath
fi

if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    echo "simulation finished? ...press enter to proceed"
    read

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py dump*.liggghts_restart

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK
    #source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh                       #- include functions
    #pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    paraview

    #- keep terminal open (if started in new terminal)
    echo "...press enter to clean up case"
    echo "press Ctr+C to keep data"
    read

fi

#- clean up case
if [ $cleanCase == "true" ]
  then
    #- clean up case
    keepDEMrestart="false"
    cleanCFDEMcase $casePath/CFD $keepDEMrestart
fi
