#!/bin/bash

#Variables section
# SALOME_INSTALL_DIR=/home/alexshtil/SALOME-9.6.0
# SALOME_INSTALL_DIR=/home/alexgubkin/SALOME-9.6.0
# SALOME_INSTALL_DIR=/home/alex/SALOME-9.7.0-native-UB20.04
SALOME_INSTALL_DIR=/home/nikita/SALOME-9.12.0-native-UB22.04-SRC
SOLUTION_DIR=$(pwd)
SALOME_SCRIPT=${SOLUTION_DIR}/src/splittedChannel.ellipse.GEOMETRY.py
# SALOME_SCRIPT=${SOLUTION_DIR}/src/splittedChannel.circle.GEOMETRY.py

initialSalomePort=2810
salomePort=${initialSalomePort}

scaleForFolderNames=6
scaleForCalculations=20

caseName=theta135

function makeCase(){
    mkdir -p ${caseName}

    cp -r src/0   ${caseName}
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

# mv /home/nikita/SALOME-9.12.0-native-UB22.04-SRC/forTopoSetBox  /home/nikita/Hele-ShawCell_KND/staticMesh/2D/${caseName}/system/

echo -e "\nGeneration done!\n"

exit 0
