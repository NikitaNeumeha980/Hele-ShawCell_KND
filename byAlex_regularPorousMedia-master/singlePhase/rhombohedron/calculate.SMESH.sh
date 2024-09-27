 
#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
SALOME_INSTALL_DIR=/home/alexshtil/SALOME-8.5.0-DB08
# SALOME_INSTALL_DIR=/home/igoshin/SALOME-8.3.0-DB08
SOLUTION_DIR=$(pwd)

scaleForFolderNames=6
scaleForCalculations=20

minTheta=$(echo "scale=${scaleForFolderNames}; 60.0/1.0;" | bc -l)
maxTheta=$(echo "scale=${scaleForFolderNames}; 90.0/1.0;" | bc -l)
thetaStep=$(echo "scale=${scaleForFolderNames}; 5.0/1.0;" | bc -l)

minIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.01/1.0;" | bc -l)
maxIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.99 - c(0.5*${maxTheta}/180*4*a(1));" | bc -l)
intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)

UTranslationDirectionNbTimes=1
VTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}
WTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}
scaleSize=$(echo "scale=${scaleForCalculations}; 0.0000001/1.0;" | bc -l)
unitCellSize=$(echo "scale=${scaleForCalculations}; 100.0 * ${scaleSize};" | bc -l)
# cellSize=$(echo "scale=${scaleForCalculations}; ${scaleSize}*${UTranslationDirectionNbTimes};" | bc -l)

work=($(find . -type d -name "theta*"))
workNum=0

#Compile some functions
cp -r src/functionObjects/ functionObjects/
./functionObjects/Allwmake

theta=${minTheta}

while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
do
    mapFieldsForFirstDirectionSwitcher=0
    mapFieldsForSecondDirectionSwitcher=0
    mapFieldsForThirdDirectionSwitcher=0

    intersectionParameter=${minIntersectionParameter}

    while [[ $(echo "scale=${scaleForFolderNames}; ${intersectionParameter} <= 0.99 - c(0.5*${theta}/180*4*a(1));" | bc -l) -ne 0 ]]
    do
        startTime=$(date +%s.%N)

        let "workNum += 1"

        CASE_DIR=theta$(echo "scale=${scaleForFolderNames}; ${theta};" | bc -l)IP$(echo "scale=${scaleForFolderNames}; ${intersectionParameter};" | bc -l)

        echo -e "\n###################################"
        echo -e "Case: "${CASE_DIR}
        echo -e "Calculation progress: " ${workNum} / ${#work[*]}
        echo -e "###################################\n"

        cd ${CASE_DIR}
        
        #Preparation case directories
        mkdir -p firstDirection secondDirection thirdDirection

        ###############################
        ##First direction calculation##
        ###############################

        #Copy sources to directory for first direction
        cp -r ../src/firstDirection/constant/ firstDirection/
        cp -r ../src/firstDirection/system/ firstDirection/

        #Preparation STL-files for snappyHexMesh
        cd firstDirection
        mv ../*.unv .
#         mv ../centreOfMass system/

        #Preparation mesh for SMESH
        runApplication ideasUnvToFoam porousCell.unv
        runApplication checkMesh -allTopology -allGeometry
        mv log.checkMesh log.checkMesh.originMesh

        #Run the polyDualMesh utility for creation dual mesh
        runApplication polyDualMesh -overwrite 180
        runApplication checkMesh -allTopology -allGeometry
        mv log.checkMesh log.checkMesh.polyDualMesh

        runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"

        #Copy mesh for another directions calculation
        cp -r constant/ ../secondDirection/constant/
        cp -r constant/ ../thirdDirection/constant/

        #Copy 0 to directory for first direction
        cp -r ../../src/firstDirection/0/ .

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

        #Run the topoSet utility for creation AMI-groups
        runApplication topoSet

        cp system/decomposeParDict.simpleFoam system/decomposeParDict

        #First approximation
        #potentialFoam
    #     if [[ ${mapFieldsForFirstDirectionSwitcher} -eq 0 ]]
    #     then
    #         runApplication decomposePar
    #         runParallel potentialFoam
    #         runApplication reconstructPar -withZero
    #         rm -r processor* log.decomposePar log.reconstructPar
    #         sed -r -i '/leftSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    #         sed -r -i '/rightSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    #         runApplication decomposePar
    # 
    #         let "mapFieldsForFirstDirectionSwitcher += 1"
    #     else
    #         #mapFields
    #         sed -r -i '/leftSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    #         sed -r -i '/rightSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    #         runApplication decomposePar
    #         runApplication mapFields ../../${PREVIOUS_CASE_DIR}/firstDirection -parallelSource -parallelTarget -sourceTime latestTime
    #     fi

        #simpleFoam calculation
        sed -r -i '/leftSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
        sed -r -i '/rightSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
        runApplication decomposePar
        runParallel simpleFoam

    #     cd ..
    # 
    #     ################################
    #     ##Second direction calculation##
    #     ################################
    # 
    #     #Copy sources to directory for second direction
    #     cp -r ../src/secondDirection/0/ secondDirection/
    #     cp -r ../src/secondDirection/system/ secondDirection/
    # 
    #     #Preparation mesh
    #     cd secondDirection
    # 
    #     #Change the case properties
    #     sed -r -i "s/.*thetaAngle.*/thetaAngle                      thetaAngle                      [0 0 0 0 0 0 0]     ${theta};/" constant/transportProperties
    #     sed -r -i "s/.*UTranslationDirectionNbTimes.*/UTranslationDirectionNbTimes    UTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${UTranslationDirectionNbTimes};/" constant/transportProperties
    #     sed -r -i "s/.*VTranslationDirectionNbTimes.*/VTranslationDirectionNbTimes    VTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${VTranslationDirectionNbTimes};/" constant/transportProperties
    #     sed -r -i "s/.*WTranslationDirectionNbTimes.*/WTranslationDirectionNbTimes    WTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${WTranslationDirectionNbTimes};/" constant/transportProperties
    #     sed -r -i "s/.*unitCellSize.*/unitCellSize                    unitCellSize                    [0 0 0 0 0 0 0]     ${unitCellSize};/" constant/transportProperties
    #     sed -r -i "s/.*L.*/L                               L                               [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes};" | bc -l);/" constant/transportProperties
    # 
    #     #Assign the translation vectors for cyclicAMI
    #     #Assign the translation vector for Left/Right side
    #     Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
    #     Ty=0
    #     Tz=0
    # 
    #     sed -r -i "/cyclicLeftSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicRightSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
    #     Ty=0
    #     Tz=0
    # 
    #     sed -r -i "/cyclicRightSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicLeftSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     #Assign the translation vector for Bottom/Top side
    #     x=$(echo "scale=${scaleForCalculations}; c(${theta}/180*4*a(1))/c(0.5*${theta}/180*4*a(1))" | bc -l)
    #     alpha=$(echo "scale=${scaleForCalculations}; 2*a(sqrt((1-${x})/(1+${x})))" | bc -l)
    #     Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*c(0.5*${theta}/180*4*a(1))" | bc -l)
    #     Ty=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*s(0.5*${theta}/180*4*a(1))" | bc -l)
    #     Tz=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${WTranslationDirectionNbTimes}*s(${alpha})" | bc -l)
    # 
    #     sed -r -i "/cyclicTopSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicBottomSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*c(0.5*${theta}/180*4*a(1))" | bc -l)
    #     Ty=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*c(${alpha})*s(0.5*${theta}/180*4*a(1))" | bc -l)
    #     Tz=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes}*s(${alpha})" | bc -l)
    # 
    #     sed -r -i "/cyclicBottomSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicTopSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     runApplication createPatch -overwrite
    #     sed -r -i '/skeletonWall/{N;N;s/patch/wall/}' constant/polyMesh/boundary
    #     runApplication topoSet
    #     cp system/decomposeParDict.simpleFoam system/decomposeParDict
    # 
    #     #First approximation
    #     #potentialFoam
    # #     if [[ $(echo "${mapFieldsForSecondDirectionSwitcher} == 0" | bc -l) -ne 0 ]]
    # #     then
    # #         runApplication decomposePar
    # #         runParallel potentialFoam
    # #         runApplication reconstructPar -withZero
    # #         rm -r processor* log.decomposePar log.reconstructPar
    # #         sed -r -i '/oppositeSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         sed -r -i '/frontSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         runApplication decomposePar
    # # 
    # #         let "mapFieldsForSecondDirectionSwitcher += 1"
    # #     else
    # #         #mapFields
    # #         sed -r -i '/oppositeSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         sed -r -i '/frontSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         runApplication decomposePar
    # #         runApplication mapFields ../../${PREVIOUS_CASE_DIR}/secondDirection -parallelSource -parallelTarget -sourceTime latestTime
    # #     fi
    # 
    #     #simpleFoam calculation
    #     runApplication decomposePar
    #     runParallel simpleFoam
    # 
    #     cd ..
    # 
    #     ###############################
    #     ##Third direction calculation##
    #     ###############################
    # 
    #     #Copy sources to directory for third direction
    #     cp -r ../src/thirdDirection/0/ thirdDirection/
    #     cp -r ../src/thirdDirection/system/ thirdDirection/
    # 
    #     #Preparation mesh
    #     cd thirdDirection
    # 
    #     #Change the case properties
    #     sed -r -i "s/.*thetaAngle.*/thetaAngle                      thetaAngle                      [0 0 0 0 0 0 0]     ${theta};/" constant/transportProperties
    #     sed -r -i "s/.*UTranslationDirectionNbTimes.*/UTranslationDirectionNbTimes    UTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${UTranslationDirectionNbTimes};/" constant/transportProperties
    #     sed -r -i "s/.*VTranslationDirectionNbTimes.*/VTranslationDirectionNbTimes    VTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${VTranslationDirectionNbTimes};/" constant/transportProperties
    #     sed -r -i "s/.*WTranslationDirectionNbTimes.*/WTranslationDirectionNbTimes    WTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${WTranslationDirectionNbTimes};/" constant/transportProperties
    #     sed -r -i "s/.*unitCellSize.*/unitCellSize                    unitCellSize                    [0 0 0 0 0 0 0]     ${unitCellSize};/" constant/transportProperties
    #     sed -r -i "s/.*L.*/L                               L                               [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${unitCellSize}*${WTranslationDirectionNbTimes};" | bc -l);/" constant/transportProperties
    # 
    #     #Assign the translation vectors for cyclicAMI
    #     #Assign the translation vector for Left/Right side
    #     Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
    #     Ty=0
    #     Tz=0
    # 
    #     sed -r -i "/cyclicLeftSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicRightSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${UTranslationDirectionNbTimes}" | bc -l)
    #     Ty=0
    #     Tz=0
    # 
    #     sed -r -i "/cyclicRightSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicLeftSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     #Assign the translation vector for Front/Opposite side
    #     Tx=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${VTranslationDirectionNbTimes}*c(${theta}/180*4*a(1));" | bc -l)
    #     Ty=$(echo "scale=${scaleForCalculations}; -${unitCellSize}*${VTranslationDirectionNbTimes}*s(${theta}/180*4*a(1));" | bc -l)
    #     Tz=0
    # 
    #     sed -r -i "/cyclicOppositeSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicFrontSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     Tx=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes}*c(${theta}/180*4*a(1));" | bc -l)
    #     Ty=$(echo "scale=${scaleForCalculations}; ${unitCellSize}*${VTranslationDirectionNbTimes}*s(${theta}/180*4*a(1));" | bc -l)
    #     Tz=0
    # 
    #     sed -r -i "/cyclicFrontSide/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
    #     neighbourPatch cyclicOppositeSide;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
    #     }" system/createPatchDict
    # 
    #     runApplication createPatch -overwrite
    #     sed -r -i '/skeletonWall/{N;N;s/patch/wall/}' constant/polyMesh/boundary
    #     runApplication topoSet
    #     cp system/decomposeParDict.simpleFoam system/decomposeParDict
    # 
    #     #First approximation
    #     #potentialFoam
    # #     if [[ $(echo "${mapFieldsForThirdDirectionSwitcher} == 0" | bc -l) -ne 0 ]]
    # #     then
    # #         runApplication decomposePar
    # #         runParallel potentialFoam
    # #         runApplication reconstructPar -withZero
    # #         rm -r processor* log.decomposePar log.reconstructPar
    # #         sed -r -i '/topSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         sed -r -i '/bottomSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         runApplication decomposePar
    # # 
    # #         let "mapFieldsForThirdDirectionSwitcher += 1"
    # #     else
    # #         #mapFields
    # #         sed -r -i '/topSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         sed -r -i '/bottomSide/{N;N;N;s/fixedValue/zeroGradient/;s/value.*$//}' 0/U
    # #         runApplication decomposePar
    # #         runApplication mapFields ../../${PREVIOUS_CASE_DIR}/thirdDirection -parallelSource -parallelTarget -sourceTime latestTime
    # #     fi
    # 
    #     #simpleFoam calculation
    #     runApplication decomposePar
    #     runParallel simpleFoam

        cd ../../

    #     PREVIOUS_CASE_DIR=${CASE_DIR}

        caseExecutionTime=$(echo "$(date +%s.%N) - ${startTime}" | bc)

        echo -e "\nCase execution time:" ${caseExecutionTime} "seconds\n"

        intersectionParameter=$(echo "scale=${scaleForCalculations}; ${intersectionParameter} + ${intersectionParameterStep};" | bc -l)
    done

    theta=$(echo "scale=${scaleForCalculations}; ${theta} + ${thetaStep};" | bc -l)
done

exit 0
