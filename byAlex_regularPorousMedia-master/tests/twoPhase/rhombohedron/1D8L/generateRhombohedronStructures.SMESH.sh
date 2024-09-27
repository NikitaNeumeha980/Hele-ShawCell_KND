#!/bin/bash

#Variables section
SALOME_INSTALL_DIR=/home/tb-itam-alexshtil/SALOME-8.5.0-DB08
SOLUTION_DIR=$(pwd)
# SALOME_SCRIPT=${SOLUTION_DIR}/src/makePorousCell.snappyHexMesh.py
SALOME_SCRIPT=${SOLUTION_DIR}/src/makePorousCell.SMESH.py

numThreads=34
threadNum=0
initialSalomePort=2810
salomePort=${initialSalomePort}

scaleForFolderNames=6
scaleForCalculations=20

minTheta=$(echo "scale=${scaleForFolderNames}; 90.0/1.0;" | bc -l)
maxTheta=$(echo "scale=${scaleForFolderNames}; 90.0/1.0;" | bc -l)
thetaStep=$(echo "scale=${scaleForFolderNames}; 5.0/1.0;" | bc -l)

minIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)
maxIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.9 - c(0.5*${maxTheta}/180*4*a(1));" | bc -l)
intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)

# scaleFactor=$(echo "scale=${scaleForCalculations}; 0.001/1.0;" | bc -l)

increasingScaleFactor=$(echo "scale=${scaleForCalculations}; 100.0/1.0;" | bc -l)

UTranslationDirectionNbTimes=8
VTranslationDirectionNbTimes=1
WTranslationDirectionNbTimes=1

numStructures=$(echo "((${maxTheta} - ${minTheta}) / ${thetaStep} + 1) * ((${maxIntersectionParameter} - ${minIntersectionParameter}) / ${intersectionParameterStep} + 1);" | bc)
structureNum=0

function make_porous_cell(){
    local CASE_DIR=theta$(echo "scale=${scaleForFolderNames}; $2;" | bc -l)IP$(echo "scale=${scaleForFolderNames}; $3;" | bc -l)

    #Preparation case directories
    mkdir -p ${CASE_DIR}

    #Execution SALOME's script
    cd ${SALOME_INSTALL_DIR}
        ./salome kill $1
        ./salome start -t -b --port=$1 ${SALOME_SCRIPT} args:${SOLUTION_DIR}/${CASE_DIR},$2,$3,${increasingScaleFactor},${UTranslationDirectionNbTimes},${VTranslationDirectionNbTimes},${WTranslationDirectionNbTimes}
        ./salome kill $1
    cd ${SOLUTION_DIR}
}

theta=${minTheta}

while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
do
    intersectionParameter=${minIntersectionParameter}

    while [[ $(echo "scale=${scaleForFolderNames}; ${intersectionParameter} <= 0.9 - c(0.5*${theta}/180*4*a(1));" | bc -l) -ne 0 ]]
    do
        make_porous_cell ${salomePort} ${theta} ${intersectionParameter} &

        let "salomePort += 1"
        let "threadNum += 1"
        let "structureNum += 1"

        if [[ $(echo "${threadNum} % ${numThreads};" | bc) -eq 0 ]]
        then
            salomePort=${initialSalomePort}
            threadNum=0
            echo -e "\nGeneration progress: " ${structureNum} / ${numStructures} "\n"
            wait
        fi

        intersectionParameter=$(echo "scale=${scaleForCalculations}; ${intersectionParameter} + ${intersectionParameterStep};" | bc -l)
    done

    theta=$(echo "scale=${scaleForCalculations}; ${theta} + ${thetaStep};" | bc -l)
done

wait

echo -e "\nGenerated: " ${structureNum} / ${numStructures}
echo -e "\nGeneration done!\n"

exit 0

