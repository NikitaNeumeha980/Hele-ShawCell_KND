#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

rm  permeabilities\
    volumetricFlows\
    massFlows

#Variables section
SALOME_INSTALL_DIR=/home/alexshtil/SALOME-8.4.0-DB08
#SALOME_INSTALL_DIR=/home/igoshin/SALOME-8.3.0-DB08
SOLUTION_DIR=$(pwd)
#SALOME_SCRIPT=${SOLUTION_DIR}/src/makeChannel.SMESH.py
SALOME_SCRIPT=${SOLUTION_DIR}/src/makeChannel.snappyHexMesh.py

scaleForFolderNames=6
scaleForCalculations=20

channelLength0=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)

Rmin0=$(echo "scale=${scaleForFolderNames}; 1.0/1.0;" | bc -l)
Rmax0=$(echo "scale=${scaleForFolderNames}; 10.0/1.0;" | bc -l)
Rstep=$(echo "scale=${scaleForFolderNames}; 0.2/1.0;" | bc -l)

scaleSize=$(echo "scale=${scaleForCalculations}; 0.000001/1.0;" | bc -l)

Rmax=${Rmin0}

touch   permeabilities\
        volumetricFlows\
        massFlows

while [[ $(echo "scale=${scaleForFolderNames}; ${Rmax} <= ${Rmax0};" | bc -l) -ne 0 ]]
do
    Rmin=${Rmin0}

    mapFieldsSwitcher=0

    while [[ $(echo "scale=${scaleForFolderNames}; ${Rmin} <= ${Rmax};" | bc -l) -ne 0 ]]
    do
        echo -e "\n###################################"
        echo "Rmax="${Rmax}
        echo "Rmin="${Rmin}
        echo -e "###################################\n"

        CASE_DIR=Rmax$(echo "scale=${scaleForFolderNames}; ${Rmax};" | bc -l)Rmin$(echo "scale=${scaleForFolderNames}; ${Rmin};" | bc -l)

        if ! [ -d ${CASE_DIR} ]
        then
            #Preparation case directories
            mkdir ${CASE_DIR}

            #Execution SALOME's script
            cd ${SALOME_INSTALL_DIR}
            ./salome kill 2811
            ./salome start -t -b --port=2811 ${SALOME_SCRIPT} args:${SOLUTION_DIR}/${CASE_DIR},${Rmin},${Rmax},${channelLength0}
            ./salome kill 2811
            cd ${SOLUTION_DIR}

            #Copy sources to directory
            cp -r src/0/ ${CASE_DIR}
            cp -r src/constant/ ${CASE_DIR}
            cp -r src/system/ ${CASE_DIR}

            #Preparation mesh for snappyHexMesh
            cd ${CASE_DIR}
            mv *.stl constant/triSurface/
            #runApplication surfaceOrient constant/triSurface/planeChannel3D.stl constant/triSurface/planeChannel3D.stl "(1e10 1e10 1e10)"
            runApplication surfaceCheck constant/triSurface/planeChannel3D.stl
            runApplication ideasUnvToFoam backgroundMesh.unv
            runApplication surfaceFeatureExtract
            runApplication surfaceFeatureConvert constant/triSurface/*.eMesh constant/triSurface/edges.vtk
            cp system/decomposeParDict.snappyHexMesh system/decomposeParDict
            mv locationInMesh system/
            runApplication decomposePar -copyZero
            runParallel snappyHexMesh -overwrite
            runApplication reconstructParMesh -constant
            rm -r processor* log.decomposePar
            runApplication extrudeMesh
            runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
            runApplication checkMesh -allGeometry

            #Preparation mesh for SMESH
#             cd ${CASE_DIR}
#             runApplication ideasUnvToFoam channel.unv
#             runApplication checkMesh -allGeometry
#             runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"

            sed -r -i "s/.*channelLength.*/channelLength   channelLength   [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${scaleSize}*${channelLength0}/1.0;" | bc -l);/" constant/transportProperties
#             sed -r -i '/.*defaultFaces.*/{N;N;N;N;N;d}' constant/polyMesh/boundary
#             sed -r -i '/wedge_1/{N;N;s/patch/wedge/}' constant/polyMesh/boundary
#             sed -r -i '/wedge_2/{N;N;s/patch/wedge/}' constant/polyMesh/boundary
#             sed -r -i '/channelWall/{N;N;s/patch/wall/}' constant/polyMesh/boundary

            cp system/decomposeParDict.simpleFoam system/decomposeParDict

            #First approximation
            #potentialFoam
#             if [[ ${mapFieldsSwitcher} -eq 0 ]]
#             then
# #                 runApplication decomposePar
# #                 runParallel potentialFoam
#                 runApplication potentialFoam
# #                runApplication reconstructPar -withZero
# #                rm -r processor* log.decomposePar log.reconstructPar
#                 sed -r -i '/inlet/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 sed -r -i '/outlet/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 runApplication decomposePar
# 
#                 let "mapFieldsSwitcher += 1"
#             else
#                 #mapFields
#                 sed -r -i '/inlet/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 sed -r -i '/outlet/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 runApplication decomposePar
#                 runApplication mapFields ../${PREVIOUS_CASE_DIR} -parallelSource -parallelTarget -sourceTime latestTime
#             fi

            #simpleFoam calculation
            sed -r -i '/inlet/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
            sed -r -i '/outlet/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
            runApplication decomposePar
            runParallel simpleFoam

            #Post-processing
            echo "${Rmax} ${Rmin} $(sed '$!d' permeabilityLUFO)" >> ../permeabilities
            echo "${Rmax} ${Rmin} $(tail -n 1 postProcessing/flowRatePatch\(name\=inlet\)/0/surfaceFieldValue.dat | sed -r 's/\t//' | sed -r -e 's/[ ]+/ /g')" >> ../volumetricFlows
#             echo "${Rmax} ${Rmin} $(tail -n 1 inputOutputMassFlowLUFO | sed -r -n '/inlet.*outlet.*relativeDiff/p' | sed -r 's/inlet.*outlet.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../massFlows

            cd ..

        else
            #First direction post-processing
            cd ${CASE_DIR}

            echo "${Rmax} ${Rmin} $(sed '$!d' permeabilityLUFO)" >> ../permeabilities
            echo "${Rmax} ${Rmin} $(tail -n 1 postProcessing/flowRatePatch\(name\=inlet\)/0/surfaceFieldValue.dat | sed -r 's/\t//' | sed -r -e 's/[ ]+/ /g')" >> ../volumetricFlows
#             echo "${Rmax} ${Rmin} $(tail -n 1 inputOutputMassFlowLUFO | sed -r -n '/inlet.*outlet.*relativeDiff/p' | sed -r 's/inlet.*outlet.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../massFlows

            cd ..
        fi

        PREVIOUS_CASE_DIR=${CASE_DIR}

        Rmin=$(echo "${Rmin} + ${Rstep}" | bc -l)
    done

    Rmax=$(echo "${Rmax} + ${Rstep}" | bc -l)
done

exit 0
