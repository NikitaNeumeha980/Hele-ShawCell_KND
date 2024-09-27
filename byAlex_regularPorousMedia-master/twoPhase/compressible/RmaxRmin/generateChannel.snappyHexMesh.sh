#!/bin/bash

#Variables section
SALOME_INSTALL_DIR=/home/alexshtil/SALOME-8.5.0-DB08
SOLUTION_DIR=$(pwd)
SALOME_SCRIPT=${SOLUTION_DIR}/src/makeChannel.snappyHexMesh.py

numThreads=34
threadNum=0
initialSalomePort=2810
salomePort=${initialSalomePort}

scaleForFolderNames=6
scaleForCalculations=20

channelLength0=$(echo "scale=${scaleForFolderNames}; 20000.0/1.0;" | bc -l)

Rmin0=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)
Rmax0=$(echo "scale=${scaleForFolderNames}; 2000.0/1.0;" | bc -l)
RMinStep=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)
RMaxStep=$(echo "scale=${scaleForFolderNames}; 20000.0/1.0;" | bc -l)

# numStructures=$(echo "((${Rmax0} - ${Rmin0}) / ${Rstep} + 1) * ((${Rmax0} - ${Rmin0}) / ${Rstep} + 2) / 2;" | bc)
structureNum=0

function make_channel(){
    local CASE_DIR=Rmax$(echo "scale=${scaleForFolderNames}; $2;" | bc -l)Rmin$(echo "scale=${scaleForFolderNames}; $3;" | bc -l)

    #Preparation case directories
    mkdir -p ${CASE_DIR}

    #Execution SALOME's script
    cd ${SALOME_INSTALL_DIR}
        ./salome kill $1
        ./salome start -t -b --port=$1 ${SALOME_SCRIPT} args:${SOLUTION_DIR}/${CASE_DIR},$2,$3,${channelLength0}
        ./salome kill $1
    cd ${SOLUTION_DIR}
}

Rmax=${Rmax0}

while [[ $(echo "scale=${scaleForFolderNames}; ${Rmax} >= ${Rmin0};" | bc -l) -ne 0 ]]
do
    Rmin=${Rmin0}

    while [[ $(echo "scale=${scaleForFolderNames}; ${Rmin} <= ${Rmax};" | bc -l) -ne 0 ]]
    do
        make_channel ${salomePort} ${Rmax} ${Rmin} &

        let "salomePort += 1"
        let "threadNum += 1"
        let "structureNum += 1"

        if [[ $(echo "${threadNum} % ${numThreads};" | bc) -eq 0 ]]
        then
            salomePort=${initialSalomePort}
            threadNum=0
#             echo -e "\nGeneration progress: " ${structureNum} / ${numStructures} "\n"
            wait
        fi

        Rmin=$(echo "${Rmin} + ${RMinStep}" | bc -l)
    done

    Rmax=$(echo "${Rmax} - ${RMaxStep}" | bc -l)
done

wait

# echo -e "\nGenerated: " ${structureNum} / ${numStructures}
echo -e "\nGeneration done!\n"

exit 0

 
