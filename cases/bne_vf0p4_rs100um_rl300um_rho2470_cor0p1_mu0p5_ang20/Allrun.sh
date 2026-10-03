#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run periodicChannel
# Christoph Goniva - Sept. 2010
#===================================================================#

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=^openib
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

mkdir -p $casePath/DEM/post/restart $casePath/DEM/post/reports

cd $casePath/CFD
blockMesh

. $casePath/parDEMrun.sh

#- run parallel CFD-DEM in new terminal
. $casePath/parCFDDEMrun.sh
