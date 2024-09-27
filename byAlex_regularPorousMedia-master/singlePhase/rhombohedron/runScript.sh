#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
SALOME_INSTALL_DIR=/home/alexshtil/SALOME-8.3.0-DB08
# SALOME_INSTALL_DIR=/home/igoshin/SALOME-8.3.0-DB08
SOLUTION_DIR=$(pwd)
SALOME_SCRIPT=${SOLUTION_DIR}/src/makePorousCell.snappyHexMesh.py
# SALOME_SCRIPT=${SOLUTION_DIR}/src/makePorousCell.SMESH.py

scaleForFolderNames=6
scaleForCalculations=20

minTheta=$(echo "scale=${scaleForFolderNames}; 60.0/1.0;" | bc -l)
maxTheta=$(echo "scale=${scaleForFolderNames}; 90.0/1.0;" | bc -l)
# thetaStep=$(echo "scale=${scaleForFolderNames}; 0.5/1.0;" | bc -l)
thetaStep=$(echo "scale=${scaleForFolderNames}; 10.0/1.0;" | bc -l)

minIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.01/1.0;" | bc -l)
maxIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.99 - c(0.5*${maxTheta}/180*4*a(1));" | bc -l)
# intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.005/1.0;" | bc -l)
intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)

UTranslationDirectionNbTimes=1
VTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}
WTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}
scaleSize=$(echo "scale=${scaleForCalculations}; 0.00001/1.0;" | bc -l)
unitCellSize=$(echo "scale=${scaleForCalculations}; ${scaleSize};" | bc -l)
# cellSize=$(echo "scale=${scaleForCalculations}; ${scaleSize}*${UTranslationDirectionNbTimes};" | bc -l)

theta=${minTheta}

#Preparation of text files
rm      permeabilityTensor\
        I1\
        I2\
        I3\
        permeabilitiesFirstDirection\
        permeabilitiesSecondDirection\
        permeabilitiesThirdDirection\
        massFlowsFirstDirection\
        massFlowsSecondDirection\
        massFlowsThirdDirection

touch   permeabilityTensor\
        I1\
        I2\
        I3\
        permeabilitiesFirstDirection\
        permeabilitiesSecondDirection\
        permeabilitiesThirdDirection\
        massFlowsFirstDirection\
        massFlowsSecondDirection\
        massFlowsThirdDirection

#Compile some functions
g++ -O3 src/permeabilityTensorCalculation.cpp -o permeabilityTensorCalculation

cp -r src/functionObjects/ functionObjects/
./functionObjects/Allwmake

while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
do
    intersectionParameter=${minIntersectionParameter}

    mapFieldsForFirstDirectionSwitcher=0
    mapFieldsForSecondDirectionSwitcher=0
    mapFieldsForThirdDirectionSwitcher=0

    while [[ $(echo "scale=${scaleForFolderNames}; ${intersectionParameter} <= 0.99 - c(0.5*${theta}/180*4*a(1));" | bc -l) -ne 0 ]]
    do
        echo -e "\n###################################"
        echo "theta="${theta}
        echo "intersectionParameter="${intersectionParameter}
        echo -e "###################################\n"

        CASE_DIR=theta$(echo "scale=${scaleForFolderNames}; ${theta};" | bc -l)IP$(echo "scale=${scaleForFolderNames}; ${intersectionParameter};" | bc -l)

        if ! [ -d ${CASE_DIR} ]
        then
            #Preparation case directories
            mkdir -p ${CASE_DIR}/firstDirection ${CASE_DIR}/secondDirection ${CASE_DIR}/thirdDirection

            #Execution SALOME's script
            cd ${SALOME_INSTALL_DIR}
            ./salome kill 2812
            ./salome start -t -b --port=2812 ${SALOME_SCRIPT} args:${SOLUTION_DIR}/${CASE_DIR}/firstDirection,${theta},${intersectionParameter},${UTranslationDirectionNbTimes},${VTranslationDirectionNbTimes},${WTranslationDirectionNbTimes}
            ./salome kill 2812
            cd ${SOLUTION_DIR}

            ###############################
            ##First direction calculation##
            ###############################

            #Copy sources to directory for first direction
            cp -r src/firstDirection/0/ ${CASE_DIR}/firstDirection/
            cp -r src/firstDirection/constant/ ${CASE_DIR}/firstDirection/
            cp -r src/firstDirection/system/ ${CASE_DIR}/firstDirection/

            #Preparation STL-files for snappyHexMesh
            cd ${CASE_DIR}/firstDirection
            mv *.stl constant/triSurface/
            mv centerOfMass system/
            cp system/snappyHexMeshSrc/* system/

            sed -i '1 s/^.*$/solid backgroundMeshBox/' constant/triSurface/backgroundMeshBox.stl
#             runApplication surfaceClean constant/triSurface/backgroundMeshBox.stl constant/triSurface/backgroundMeshBox.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid skeletonWall/' constant/triSurface/skeletonWall.stl
#             runApplication surfaceClean constant/triSurface/skeletonWall.stl constant/triSurface/skeletonWall.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid leftSide/' constant/triSurface/leftSide.stl
#             runApplication surfaceClean constant/triSurface/leftSide.stl constant/triSurface/leftSide.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid rightSide/' constant/triSurface/rightSide.stl
#             runApplication surfaceClean constant/triSurface/rightSide.stl constant/triSurface/rightSide.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid frontSide/' constant/triSurface/frontSide.stl
#             runApplication surfaceClean constant/triSurface/frontSide.stl constant/triSurface/frontSide.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid oppositeSide/' constant/triSurface/oppositeSide.stl
#             runApplication surfaceClean constant/triSurface/oppositeSide.stl constant/triSurface/oppositeSide.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid topSide/' constant/triSurface/topSide.stl
#             runApplication surfaceClean constant/triSurface/topSide.stl constant/triSurface/topSide.stl 1e-6 1e-6
#             rm log.surfaceClean

            sed -i '1 s/^.*$/solid bottomSide/' constant/triSurface/bottomSide.stl
#             runApplication surfaceClean constant/triSurface/bottomSide.stl constant/triSurface/bottomSide.stl 1e-6 1e-6
#             rm log.surfaceClean

            touch constant/triSurface/porousCell.stl

            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/leftSide.stl constant/triSurface/porousCell.stl
            rm log.surfaceAdd
            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/rightSide.stl constant/triSurface/porousCell.stl
            rm log.surfaceAdd
            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/frontSide.stl constant/triSurface/porousCell.stl
            rm log.surfaceAdd
            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/oppositeSide.stl constant/triSurface/porousCell.stl
            rm log.surfaceAdd
            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/topSide.stl constant/triSurface/porousCell.stl
            rm log.surfaceAdd
            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/bottomSide.stl constant/triSurface/porousCell.stl
            rm log.surfaceAdd
            runApplication surfaceAdd constant/triSurface/porousCell.stl constant/triSurface/skeletonWall.stl constant/triSurface/porousCell.stl

#             cat \
#                 constant/triSurface/leftSide.stl \
#                 constant/triSurface/rightSide.stl \
#                 constant/triSurface/topSide.stl \
#                 constant/triSurface/bottomSide.stl \
#                 constant/triSurface/oppositeSide.stl \
#                 constant/triSurface/frontSide.stl \
#                 constant/triSurface/skeletonWall.stl > \
#                 constant/triSurface/porousCell.stl

#             runApplication surfaceOrient constant/triSurface/porousCell.stl constant/triSurface/porousCell.stl '(0 0 0)'
#             runApplication surfaceClean constant/triSurface/porousCell.stl constant/triSurface/porousCell.stl 1e-6 1e-6
#             runApplication surfaceCheck constant/triSurface/porousCell.stl
            runApplication surfaceFeatureExtract
            runApplication surfaceFeatureConvert constant/triSurface/*.eMesh constant/triSurface/edges.vtk

            #Preparation background mesh for snappyHexMesh
            runApplication ideasUnvToFoam backgroundMesh.unv

            #Run snappyHexMesh
            runApplication decomposePar -copyZero
            runParallel snappyHexMesh -overwrite
            runParallel collapseEdges -overwrite
#             rm log.collapseEdges
#             runParallel collapseEdges -overwrite -collapseFaces
            runParallel checkMesh -allTopology -allGeometry
            runParallel transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
            runApplication reconstructParMesh -constant
#             runApplication checkMesh -allTopology -allGeometry
            rm -r processor* log.decomposePar

            #Preparation mesh for SMESH
#             runApplication ideasUnvToFoam porousCell.unv
#             runApplication checkMesh -allGeometry
#             runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
            
            #Copy mesh for another directions calculation
            cp -r constant/ ../secondDirection/constant/
            cp -r constant/ ../thirdDirection/constant/

            #Change the case properties
            sed -r -i "s/.*thetaAngle.*/thetaAngle                      thetaAngle                      [0 0 0 0 0 0 0]     ${theta};/" constant/transportProperties
            sed -r -i "s/.*UTranslationDirectionNbTimes.*/UTranslationDirectionNbTimes    UTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${UTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*VTranslationDirectionNbTimes.*/VTranslationDirectionNbTimes    VTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${VTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*WTranslationDirectionNbTimes.*/WTranslationDirectionNbTimes    WTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${WTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*unitCellSize.*/unitCellSize                    unitCellSize                    [0 0 0 0 0 0 0]     ${unitCellSize};/" constant/transportProperties
            sed -r -i "s/.*L.*/L                               L                               [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${unitCellSize}*${UTranslationDirectionNbTimes};" | bc -l);/" constant/transportProperties

            #Assign the translation vectors for cyclicAMI
            #Assign the translation vector for Front/Opposite side
            Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${VTranslationDirectionNbTimes}*c(${theta}/180*4*a(1));" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${VTranslationDirectionNbTimes}*s(${theta}/180*4*a(1));" | bc -l)
            Tz=0

            sed -r -i "/cyclicOppositeSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicFrontSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes}*c(${theta}/180*4*a(1));" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes}*s(${theta}/180*4*a(1));" | bc -l)
            Tz=0

            sed -r -i "/cyclicFrontSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicOppositeSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            #Assign the translation vector for Bottom/Top side
            x=$(echo "scale=${scaleForCalculations}; c(${theta}/180*4*a(1))/c(0.5*${theta}/180*4*a(1));" | bc -l)
            alpha=$(echo "scale=${scaleForCalculations}; 2*a(sqrt((1-${x})/(1+${x})));" | bc -l)
            Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*c(0.5*${theta}/180*4*a(1));" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*s(0.5*${theta}/180*4*a(1));" | bc -l)
            Tz=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*s(${alpha});" | bc -l)

            sed -r -i "/cyclicTopSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicBottomSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*c(0.5*${theta}/180*4*a(1));" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*s(0.5*${theta}/180*4*a(1));" | bc -l)
            Tz=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*s(${alpha});" | bc -l)

            sed -r -i "/cyclicBottomSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicTopSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            #Run the createPatch utility for creation cyclicAMI
            runApplication createPatch -overwrite
            sed -r -i '/skeletonWall/{N;N;s/patch/wall/}' constant/polyMesh/boundary
            runApplication topoSet
            cp system/decomposeParDict.simpleFoam system/decomposeParDict

            #First approximation
            #potentialFoam
#             if [[ ${mapFieldsForFirstDirectionSwitcher} -eq 0 ]]
#             then
#                 runApplication decomposePar
#                 runParallel potentialFoam
#                 runApplication reconstructPar -withZero
#                 rm -r processor* log.decomposePar log.reconstructPar
#                 sed -r -i '/leftSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 sed -r -i '/rightSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 runApplication decomposePar
# 
#                 let "mapFieldsForFirstDirectionSwitcher += 1"
#             else
#                 #mapFields
#                 sed -r -i '/leftSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 sed -r -i '/rightSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 runApplication decomposePar
#                 runApplication mapFields ../../${PREVIOUS_CASE_DIR}/firstDirection -parallelSource -parallelTarget -sourceTime latestTime
#             fi

            #simpleFoam calculation
            sed -r -i '/leftSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
            sed -r -i '/rightSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
            runApplication decomposePar
            runParallel simpleFoam

            #Post-processing for LUFO-functions
#             echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesFirstDirection
#             echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsFirstDirection
#             echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAveragedGradPDivNu
#             echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../averagedU

            #Post-processing for functionObjects
            echo \
                "${theta} ${intersectionParameter}" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=leftSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1)" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=rightSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicTopSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicBottomSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicOppositeSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicFrontSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                >> ../../massFlowsFirstDirection
            echo "$(sed '$!d' postProcessing/averagedU1/0/averagedU.dat | sed -r 's/\t/\n/g' | tail -3)" >> ../averagedU
            echo "$(sed '$!d' postProcessing/mAveragedGradPDivNu1/0/mAveragedGradPDivNu.dat | sed -r 's/\t/\n/g' | tail -3)" >> ../mAveragedGradPDivNu
            echo "${theta} ${intersectionParameter} $(sed '$!d' postProcessing/rhombohedronStructurePermeability1/0/rhombohedronStructurePermeability.dat | sed -r 's/\t/\n/g' | tail -1)" >> ../../permeabilitiesFirstDirection

            cd ..

            ################################
            ##Second direction calculation##
            ################################

            #Copy sources to directory for second direction
            cp -r ../src/secondDirection/0/ secondDirection/
            cp -r ../src/secondDirection/system/ secondDirection/

            #Preparation mesh
            cd secondDirection

            #Change the case properties
            sed -r -i "s/.*thetaAngle.*/thetaAngle                      thetaAngle                      [0 0 0 0 0 0 0]     ${theta};/" constant/transportProperties
            sed -r -i "s/.*UTranslationDirectionNbTimes.*/UTranslationDirectionNbTimes    UTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${UTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*VTranslationDirectionNbTimes.*/VTranslationDirectionNbTimes    VTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${VTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*WTranslationDirectionNbTimes.*/WTranslationDirectionNbTimes    WTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${WTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*unitCellSize.*/unitCellSize                    unitCellSize                    [0 0 0 0 0 0 0]     ${unitCellSize};/" constant/transportProperties
            sed -r -i "s/.*L.*/L                               L                               [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes};" | bc -l);/" constant/transportProperties

            #Assign the translation vectors for cyclicAMI
            #Assign the translation vector for Left/Right side
            Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
            Ty=0
            Tz=0

            sed -r -i "/cyclicLeftSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicRightSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
            Ty=0
            Tz=0

            sed -r -i "/cyclicRightSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicLeftSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            #Assign the translation vector for Bottom/Top side
            x=$(echo "scale=${scaleForCalculations}; c(${theta}/180*4*a(1))/c(0.5*${theta}/180*4*a(1))" | bc -l)
            alpha=$(echo "scale=${scaleForCalculations}; 2*a(sqrt((1-${x})/(1+${x})))" | bc -l)
            Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*c(0.5*${theta}/180*4*a(1))" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*s(0.5*${theta}/180*4*a(1))" | bc -l)
            Tz=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*s(${alpha})" | bc -l)

            sed -r -i "/cyclicTopSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicBottomSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*c(0.5*${theta}/180*4*a(1))" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*s(0.5*${theta}/180*4*a(1))" | bc -l)
            Tz=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*s(${alpha})" | bc -l)

            sed -r -i "/cyclicBottomSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicTopSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            runApplication createPatch -overwrite
            sed -r -i '/skeletonWall/{N;N;s/patch/wall/}' constant/polyMesh/boundary
            runApplication topoSet
            cp system/decomposeParDict.simpleFoam system/decomposeParDict

            #First approximation
            #potentialFoam
#         if [[ $(echo "${mapFieldsForSecondDirectionSwitcher} == 0" | bc -l) -ne 0 ]]
#         then
#             runApplication decomposePar
#             runParallel potentialFoam
#             runApplication reconstructPar -withZero
#             rm -r processor* log.decomposePar log.reconstructPar
#             sed -r -i '/oppositeSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#             sed -r -i '/frontSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#             runApplication decomposePar
# 
#             let "mapFieldsForSecondDirectionSwitcher += 1"
#         else
#             #mapFields
#             sed -r -i '/oppositeSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#             sed -r -i '/frontSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#             runApplication decomposePar
#             runApplication mapFields ../../${PREVIOUS_CASE_DIR}/secondDirection -parallelSource -parallelTarget -sourceTime latestTime
#         fi
# 
            #simpleFoam calculation
            runApplication decomposePar
            runParallel simpleFoam

            #Post-processing for LUFO-functions
#             echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesSecondDirection
#             echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsSecondDirection
#             echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
#             echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

            #Post-processing for functionObjects
            echo \
                "${theta} ${intersectionParameter}" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicLeftSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1)" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicRightSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicTopSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicBottomSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=oppositeSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=frontSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                >> ../../massFlowsSecondDirection
            echo "$(sed '$!d' postProcessing/averagedU1/0/averagedU.dat | sed -r 's/\t/\n/g' | tail -3)" >> ../averagedU
            echo "$(sed '$!d' postProcessing/mAveragedGradPDivNu1/0/mAveragedGradPDivNu.dat | sed -r 's/\t/\n/g' | tail -3)" >> ../mAveragedGradPDivNu
            echo "${theta} ${intersectionParameter} $(sed '$!d' postProcessing/rhombohedronStructurePermeability1/0/rhombohedronStructurePermeability.dat | sed -r 's/\t/\n/g' | tail -1)" >> ../../permeabilitiesSecondDirection

            cd ..

            ###############################
            ##Third direction calculation##
            ###############################

            #Copy sources to directory for third direction
            cp -r ../src/thirdDirection/0/ thirdDirection/
            cp -r ../src/thirdDirection/system/ thirdDirection/

            #Preparation mesh
            cd thirdDirection

            #Change the case properties
            sed -r -i "s/.*thetaAngle.*/thetaAngle                      thetaAngle                      [0 0 0 0 0 0 0]     ${theta};/" constant/transportProperties
            sed -r -i "s/.*UTranslationDirectionNbTimes.*/UTranslationDirectionNbTimes    UTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${UTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*VTranslationDirectionNbTimes.*/VTranslationDirectionNbTimes    VTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${VTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*WTranslationDirectionNbTimes.*/WTranslationDirectionNbTimes    WTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${WTranslationDirectionNbTimes};/" constant/transportProperties
            sed -r -i "s/.*unitCellSize.*/unitCellSize                    unitCellSize                    [0 0 0 0 0 0 0]     ${unitCellSize};/" constant/transportProperties
            sed -r -i "s/.*L.*/L                               L                               [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes};" | bc -l);/" constant/transportProperties

            #Assign the translation vectors for cyclicAMI
            #Assign the translation vector for Left/Right side
            Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
            Ty=0
            Tz=0

            sed -r -i "/cyclicLeftSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicRightSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
            Ty=0
            Tz=0

            sed -r -i "/cyclicRightSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicLeftSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            #Assign the translation vector for Front/Opposite side
            Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${VTranslationDirectionNbTimes}*c(${theta}/180*4*a(1));" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${VTranslationDirectionNbTimes}*s(${theta}/180*4*a(1));" | bc -l)
            Tz=0

            sed -r -i "/cyclicOppositeSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicFrontSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes}*c(${theta}/180*4*a(1));" | bc -l)
            Ty=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes}*s(${theta}/180*4*a(1));" | bc -l)
            Tz=0

            sed -r -i "/cyclicFrontSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
            neighbourPatch cyclicOppositeSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
            }" system/createPatchDict

            runApplication createPatch -overwrite
            sed -r -i '/skeletonWall/{N;N;s/patch/wall/}' constant/polyMesh/boundary
            runApplication topoSet
            cp system/decomposeParDict.simpleFoam system/decomposeParDict

            #First approximation
            #potentialFoam
#             if [[ $(echo "${mapFieldsForThirdDirectionSwitcher} == 0" | bc -l) -ne 0 ]]
#             then
#                 runApplication decomposePar
#                 runParallel potentialFoam
#                 runApplication reconstructPar -withZero
#                 rm -r processor* log.decomposePar log.reconstructPar
#                 sed -r -i '/topSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 sed -r -i '/bottomSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 runApplication decomposePar
# 
#                 let "mapFieldsForThirdDirectionSwitcher += 1"
#             else
#                 #mapFields
#                 sed -r -i '/topSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 sed -r -i '/bottomSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
#                 runApplication decomposePar
#                 runApplication mapFields ../../${PREVIOUS_CASE_DIR}/thirdDirection -parallelSource -parallelTarget -sourceTime latestTime
#             fi

            #simpleFoam calculation
            runApplication decomposePar
            runParallel simpleFoam

            #Post-processing for LUFO-functions
#             echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesThirdDirection
#             echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsThirdDirection
#             echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
#             echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

            #Post-processing for functionObjects
            echo \
                "${theta} ${intersectionParameter}" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicLeftSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1)" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicRightSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=topSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=bottomSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicOppositeSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                "$(sed '$!d' postProcessing/flowRatePatch\(name=cyclicFrontSide\)/0/surfaceFieldValue.dat | sed -r 's/\t/\n/g' | tail -1 | sed 's/^//')" \
                >> ../../massFlowsThirdDirection
            echo "$(sed '$!d' postProcessing/averagedU1/0/averagedU.dat | sed -r 's/\t/\n/g' | tail -3)" >> ../averagedU
            echo "$(sed '$!d' postProcessing/mAveragedGradPDivNu1/0/mAveragedGradPDivNu.dat | sed -r 's/\t/\n/g' | tail -3)" >> ../mAveragedGradPDivNu
            echo "${theta} ${intersectionParameter} $(sed '$!d' postProcessing/rhombohedronStructurePermeability1/0/rhombohedronStructurePermeability.dat | sed -r 's/\t/\n/g' | tail -1)" >> ../../permeabilitiesThirdDirection

            cd ../../

            ./permeabilityTensorCalculation ${CASE_DIR}

            echo "${theta} ${intersectionParameter}" >> permeabilityTensor
            echo "$(sed -n '1,3p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> permeabilityTensor
            echo "${theta} ${intersectionParameter} $(sed -n '4p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I1
            echo "${theta} ${intersectionParameter} $(sed -n '5p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I2
            echo "${theta} ${intersectionParameter} $(sed -n '6p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I3
        else
            #First direction post-processing
            cd ${CASE_DIR}/firstDirection

            echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesFirstDirection
            echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsFirstDirection
            echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
            echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

            cd ..

            #Second direction post-processing
            cd secondDirection

            echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesSecondDirection
            echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsSecondDirection
            echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
            echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

            cd ..

            #Third direction post-processing
            cd thirdDirection

            echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesThirdDirection
            echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsThirdDirection
            echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
            echo "$(tail -n 4 mAGradPDivNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

            cd ../../

            ./permeabilityTensorCalculation ${CASE_DIR}

            echo "${theta} ${intersectionParameter}" >> permeabilityTensor
            echo "$(sed -n '1,3p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> permeabilityTensor
            echo "${theta} ${intersectionParameter} $(sed -n '4p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I1
            echo "${theta} ${intersectionParameter} $(sed -n '5p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I2
            echo "${theta} ${intersectionParameter} $(sed -n '6p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I3
        fi

        PREVIOUS_CASE_DIR=${CASE_DIR}

        intersectionParameter=$(echo "scale=${scaleForCalculations}; ${intersectionParameter} + ${intersectionParameterStep};" | bc -l)
    done

    theta=$(echo "scale=${scaleForCalculations}; ${theta} + ${thetaStep};" | bc -l)
done

exit 0
