#!/bin/bash

#Variables section
SALOME_INSTALL_DIR=/home/alexshtil/SALOME-9.6.0-DB08
# SALOME_INSTALL_DIR=/home/itam-alex/SALOME-9.6.0
SOLUTION_DIR=$(pwd)
SALOME_SCRIPT=${SOLUTION_DIR}/src/makeLattice.GEOMETRY.py

initialSalomePort=2880
salomePort=${initialSalomePort}

scaleForFolderNames=6
scaleForCalculations=20

caseName=case2

function makeCase(){
    mkdir -p ${caseName}

    cp -r src/0.orig    ${caseName}
    cp -r src/constant  ${caseName}
    cp -r src/system    ${caseName}
}

function makeLattice(){
    #Execution SALOME's script
    cd ${SALOME_INSTALL_DIR}
        ./salome kill $1
        ./salome start -t -b --port=$1 ${SALOME_SCRIPT} args:${SOLUTION_DIR}/${caseName}
        ./salome kill $1
    cd ${SOLUTION_DIR}
}

makeCase

makeLattice ${salomePort}

echo -e "\nGeneration done!\n"

exit 0
