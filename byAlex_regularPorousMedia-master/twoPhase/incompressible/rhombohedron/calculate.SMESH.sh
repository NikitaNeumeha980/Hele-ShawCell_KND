#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
# SALOME_INSTALL_DIR=/home/alexshtil/SALOME-8.5.0-DB08
SOLUTION_DIR=$(pwd)

scaleForFolderNames=6
scaleForCalculations=20

pi=$(echo "scale=${scaleForCalculations}; 4*a(1);" | bc -l)

minTheta=$(echo "scale=${scaleForFolderNames}; 78.0/1.0;" | bc -l)
maxTheta=$(echo "scale=${scaleForFolderNames}; 78.0/1.0;" | bc -l)
thetaStep=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)

# minIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)
# intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)

# minPorosity=$(echo "scale=${scaleForFolderNames}; 0.2.0;" | bc -l)
# maxPorosity=$(echo "scale=${scaleForFolderNames}; 0.2/1.0;" | bc -l)
# porosityStep=$(echo "scale=${scaleForFolderNames}; 0.05/1.0;" | bc -l)

scaleFactor=$(echo "scale=${scaleForCalculations}; 0.001/1.0;" | bc -l)

increasingScaleFactor=$(echo "scale=${scaleForCalculations}; 100.0/1.0;" | bc -l)

scaledScaleFactor=$(echo "scale=${scaleForCalculations}; ${scaleFactor}/${increasingScaleFactor};" | bc -l)

# scaleSize=$(echo "scale=${scaleForCalculations}; ${scaleFactor}/${increasingScaleFactor};" | bc -l)

unitCellSize=$(echo "scale=${scaleForCalculations}; ${scaleFactor};" | bc -l)

# cellSize=$(echo "scale=${scaleForCalculations}; ${scaleSize}*${UTranslationDirectionNbTimes};" | bc -l)

minDropletRadius=$(echo "scale=${scaleForFolderNames}; 0.1 * ${unitCellSize};" | bc -l)
maxDropletRadius=$(echo "scale=${scaleForFolderNames}; 0.5 * ${unitCellSize};" | bc -l)
dropletRadiusStep=$(echo "scale=${scaleForFolderNames}; 0.05 * ${unitCellSize};" | bc -l)

UTranslationDirectionNbTimes=1
VTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}
WTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}

work=($(find . -type d -name "theta*"))
workNum=0

#Functions section
# function min(){
#     if (( $(echo "${1} <= ${2}" | bc -l) )); then
#         echo ${1}
#     else
#         echo ${2}
#     fi
# }
# 
# function arccos()
# {
#     scale=${scaleForCalculations}
#     if (( $(echo "${1} == 0" | bc -l) )); then
#         echo "scale=${scaleForCalculations}; 0.5*${pi}" | bc -l
#     elif (( $(echo "(-1 <= ${1}) && (${1} < 0)" | bc -l) )); then
#         echo "scale=${scaleForCalculations}; ${pi} - a(sqrt((1/(${1}^2)) - 1))" | bc -l
#     elif (( $(echo "(0 < ${1}) && (${1} <= 1)" | bc -l) )); then
#         echo "scale=${scaleForCalculations}; a(sqrt((1/(${1}^2)) - 1))" | bc -l
#     else
#         echo "input out of range"
#     fi
# }
# 
# function porosity()
# {
#     IP=${1}
#     theta=$(echo "scale=${scaleForCalculations}; ${2}/180*${pi};" | bc -l)
# 
# #     touchIP=$(echo "scale=${scaleForCalculations}; 1 - 0.5/s(${theta}/2);" | bc -l)
# 
#     IP2=$(echo "scale=${scaleForCalculations}; 1 - 2*(1 - ${IP})*s(${theta}/2);" | bc -l)
# 
#     if (( $(echo "scale=${scaleForCalculations}; ${IP2} > 0;" | bc -l) )); then
#         IP2=${IP2}
#     else
#         IP2=0
#     fi
# #     IP2=$(echo "scale=${scaleForCalculations}; 1 - 2*(1 - ${IP})*s(${theta}/2);" | bc -l)
# #     IP2=$(echo "scale=${scaleForCalculations}; if (${touchIP} > 0) ${IP2} else 0;" | bc -l)
#     porosity=$(echo "scale=${scaleForCalculations}; 1 - ${pi}/12*(2 - 3*${IP}*${IP}*(3 - ${IP}) - 3*${IP2}*${IP2}*(3 - ${IP2}))/((1 - ${IP})*(1 - ${IP})*(1 - ${IP})*(1 - c(${theta}))*sqrt(1 + 2*c(${theta})));" | bc -l)
#     echo ${porosity}
# }
# 
# function IP()
# {
#     #m - first argument, theta - second argument
#     porosity=${1}
#     theta=$(echo "scale=${scaleForCalculations}; ${2}/180*${pi};" | bc -l)
# 
# #     critIP=$(echo "scale=${scaleForCalculations}; 1 - c(${theta}/2);" | bc -l)
# #     critPorosity=$(porosity ${critIP} ${2})
# 
#     touchIP=$(echo "scale=${scaleForCalculations}; 1 - 0.5/s(${theta}/2);" | bc -l)
#     touchPorosity=$(porosity ${touchIP} ${2})
# #     IP2=$(echo "scale=${scaleForCalculations}; 1 - 2*(1 - ${IP})*s(${theta}/2);" | bc -l)
# #     IP2=$(echo "scale=${scaleForCalculations}; if (${touchIP} > 0) ${IP2} else 0;" | bc -l)
# 
#     F=$(echo "scale=${scaleForCalculations}; 12/${pi}*(1 - ${porosity})*(1 - c(${theta}))*sqrt(1 + 2*c(${theta}));" | bc -l)
# 
#     if (( $(echo "scale=${scaleForCalculations}; ${touchPorosity} > ${porosity};" | bc -l) )); then
#         a3=$(echo "scale=${scaleForCalculations}; ${F} + 24*s(${theta}/2)*s(${theta}/2)*s(${theta}/2) + 3;" | bc -l)
#         a2=$(echo "scale=${scaleForCalculations}; (-3*${F} - 72*s(${theta}/2)^3 - 9)/${a3};" | bc -l)
#         a1=$(echo "scale=${scaleForCalculations}; (3*${F} + 72*s(${theta}/2)*s(${theta}/2)*s(${theta}/2) - 18*s(${theta}/2))/${a3};" | bc -l)
#         a0=$(echo "scale=${scaleForCalculations}; (-${F} - 24*s(${theta}/2)*s(${theta}/2)*s(${theta}/2) + 18*s(${theta}/2) - 4)/${a3};" | bc -l)
#     else
#         a3=$(echo "scale=${scaleForCalculations}; ${F} + 3;" | bc -l)
#         a2=$(echo "scale=${scaleForCalculations}; (-3*${F} - 9)/${a3};" | bc -l)
#         a1=$(echo "scale=${scaleForCalculations}; 3*${F}/${a3};" | bc -l)
#         a0=$(echo "scale=${scaleForCalculations}; (2 - ${F})/${a3};" | bc -l)
#     fi
# 
#     Q=$(echo "scale=${scaleForCalculations}; (${a2}*${a2} - 3*${a1})/9;" | bc -l)
#     R=$(echo "scale=${scaleForCalculations}; (2*${a2}*${a2}*${a2} - 9*${a2}*${a1} + 27*${a0})/54;" | bc -l)
# 
#     S=$(echo "scale=${scaleForCalculations}; ${Q}*${Q}*${Q} - ${R}*${R};" | bc -l)
# 
#     if (( $(echo "scale=${scaleForCalculations}; ${S} > 0;" | bc -l) )); then
#         #Three real roots
#         arg=$(echo "scale=${scaleForCalculations}; ${R}/sqrt(${Q}*${Q}*${Q});" | bc -l)
#         phi=$(echo "scale=${scaleForCalculations}; 1/3*$(arccos ${arg});" | bc -l)
# 
#         r1=$(echo "scale=${scaleForCalculations}; -2*sqrt(${Q})*c(${phi}) - ${a2}/3;" | bc -l)
#         r2=$(echo "scale=${scaleForCalculations}; -2*sqrt(${Q})*c(${phi} + 2/3*${pi}) - ${a2}/3;" | bc -l)
#         r3=$(echo "scale=${scaleForCalculations}; -2*sqrt(${Q})*c(${phi} - 2/3*${pi}) - ${a2}/3;" | bc -l)
# 
#         echo $(min $(min ${r1} ${r2}) ${r3})
#     elif (( $(echo "scale=${scaleForCalculations}; ${S} < 0;" | bc -l) )); then
#         #One real root
#         if (( $(echo "scale=${scaleForCalculations}; ${Q} > 0;" | bc -l) )); then
#             arg=$(echo "scale=${scaleForCalculations}; sqrt(${R}*${R})/sqrt(${Q}*${Q}*${Q});" | bc -l)
#             phi=$(echo "scale=${scaleForCalculations}; 1/3*l(${arg} + sqrt(${arg}*${arg} - 1))" | bc -l)            
# 
#             r1=$(echo "scale=${scaleForCalculations}; -2*${R}/sqrt(${R}*${R})*sqrt(${Q})*0.5*(e(${phi}) + e(-${phi})) - ${a2}/3;" | bc -l)
#             echo ${r1}
#         elif (( $(echo "scale=${scaleForCalculations}; ${Q} < 0;" | bc -l) )); then
#             arg=$(echo "scale=${scaleForCalculations}; sqrt(${R}*${R})/e(l(sqrt(${Q}*${Q}))*3/2);" | bc -l)
#             phi=$(echo "scale=${scaleForCalculations}; 1/3*l(${arg} + sqrt(${arg}*${arg} + 1))" | bc -l)
# 
#             r1=$(echo "scale=${scaleForCalculations}; -2*${R}/sqrt(${R}*${R})*sqrt(${Q})*0.5*(e(${phi}) - e(-${phi})) - ${a2}/3;" | bc -l)
#             echo ${r1}
#         else
#             r1=$(echo "scale=${scaleForCalculations}; -e(l(${a0} - ${a2}*${a2}*${a2}/27)*1/3) - ${a2}/3;" | bc -l)
#             echo ${r1}
#         fi
#     else
#         #Two real roots
#         r1=$(echo "scale=${scaleForCalculations}; -2*${R}/sqrt(${R}*${R})*sqrt(${Q}) - ${a2}/3;" | bc -l)
#         r2=$(echo "scale=${scaleForCalculations}; ${R}/sqrt(${R}*${R})*sqrt(${Q}) - ${a2}/3;" | bc -l)
#         echo $(min ${r1} ${r2})
#     fi
# }

#Compile some functions
# cp -r src/functionObjects/ functionObjects/
# ./functionObjects/Allwmake

theta=${minTheta}

#Intersection parameter based loop
# while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
# do
#     intersectionParameter=${minIntersectionParameter}
# 
#     while [[ $(echo "scale=${scaleForFolderNames}; ${intersectionParameter} <= 0.9 - c(0.5*${theta}/180*4*a(1));" | bc -l) -ne 0 ]]
#     do

#Porosity based loop
while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
do
#     porosity=${minPorosity}

#     while [[ $(echo "scale=${scaleForFolderNames}; ${porosity} <= ${maxPorosity};" | bc -l) -ne 0 ]]
#     do
    for POROSITY_DIR in theta${theta}*/ ; do
        startTime=$(date +%s.%N)

        let "workNum += 1"

#         CASE_DIR=theta$(echo "scale=${scaleForFolderNames}; ${theta};" | bc -l)IP$(IP ${porosity} ${theta})

        CASE_DIR=${POROSITY_DIR}

        echo -e "\n###################################"
        echo -e "Case: "${CASE_DIR}
        echo -e "Calculation progress: " ${workNum} / ${#work[*]}
        echo -e "###################################\n"

        cd ${CASE_DIR}

        #Preparation case directories
#         mkdir -p firstDirection secondDirection thirdDirection

        onePoreCentreOfMassX=$(sed 's/.*(\(.*\)).*/\1/' onePoreCentreOfMass | sed 's/^[ \t]*//;s/[ \t]*$//g' | sed 's/[ ]\+/\n/g' | sed -n '1p')
        onePoreCentreOfMassY=$(sed 's/.*(\(.*\)).*/\1/' onePoreCentreOfMass | sed 's/^[ \t]*//;s/[ \t]*$//g' | sed 's/[ ]\+/\n/g' | sed -n '2p')
        onePoreCentreOfMassZ=$(sed 's/.*(\(.*\)).*/\1/' onePoreCentreOfMass | sed 's/^[ \t]*//;s/[ \t]*$//g' | sed 's/[ ]\+/\n/g' | sed -n '3p')

#         echo "centre ("\
#             $(echo "scale=${scaleForCalculations}; ${scaledScaleFactor}*${onePoreCentreOfMassX};" | bc -l) \
#             $(echo "scale=${scaleForCalculations}; ${scaledScaleFactor}*${onePoreCentreOfMassY};" | bc -l) \
#             $(echo "scale=${scaleForCalculations}; ${scaledScaleFactor}*${onePoreCentreOfMassZ};" | bc -l)");" > centreOfMassScaled

        dropletRadius=${minDropletRadius}

        builtMeshSwitcher=0

        while [[ $(echo "${dropletRadius} <= ${maxDropletRadius};" | bc -l) -ne 0 ]]
        do
            SUB_CASE_DIR=dropletRadius$(echo "scale=${scaleForFolderNames}; ${dropletRadius};" | bc -l)

            mkdir -p ${SUB_CASE_DIR}

            cd ${SUB_CASE_DIR}

            #Copy sources to directory for first direction
#             cp -r ../src/firstDirection/constant/ firstDirection/
#             cp -r ../src/firstDirection/system/ firstDirection/
            cp -r ../../src/firstDirection/constant/ .
            cp -r ../../src/firstDirection/system/ .
            
            cp system/decomposeParDict.interFoam system/decomposeParDict

            cp system/controlDict.interFoam system/controlDict
            cp system/fvSchemes.interFoam system/fvSchemes
            cp system/fvSolution.interFoam system/fvSolution

#             cp system/controlDict.interFlow system/controlDict
#             cp system/fvSchemes.interFlow system/fvSchemes
#             cp system/fvSolution.interFlow system/fvSolution

#             cp ../centreOfMassScaled system/
            for (( i=1; i <= ${UTranslationDirectionNbTimes}; i++ ))
            do
                for (( j=1; j <= ${VTranslationDirectionNbTimes}; j++ ))
                do
                    for (( k=1; k <= ${WTranslationDirectionNbTimes}; k++ ))
                    do
                        echo -e \
                            "sphereToCell\n"\
                            "{\n"\
                            "   centre ("\
                                    $(echo "scale=${scaleForCalculations}; ${scaledScaleFactor}*${onePoreCentreOfMassX} + ${unitCellSize}*(${i} - 1);" | bc -l) \
                                    $(echo "scale=${scaleForCalculations}; ${scaledScaleFactor}*${onePoreCentreOfMassY} + ${unitCellSize}*(${j} - 1);" | bc -l) \
                                    $(echo "scale=${scaleForCalculations}; ${scaledScaleFactor}*${onePoreCentreOfMassZ} + ${unitCellSize}*(${k} - 1);" | bc -l)");\n"\
                            "   radius $(echo "scale=${scaleForCalculations}; ${dropletRadius};" | bc -l);\n"\
                            "   fieldValues\n"\
                            "   (\n"\
                            "       volScalarFieldValue alpha.phase1 1\n"\
                            "   );\n"\
                            "}\n" >> system/droplets
                    done
                done
            done
    
#             #Change the case properties
#             sed -r -i "/sphereToCell/{N;N;N;s/centre.*/centre                  $(echo $(cat ../centreOfMass) | sed -r 's/.*centre //')/;/.*centre.*/a\
#             radius                  $(echo "scale=${scaleForCalculations}; ${dropletRadius};" | bc -l);\n
#             }" system/setFieldsDict

            #Change the case properties
#             sed -r -i "/sphereToCell/{N;N;N;s/radius.*/radius                  $(echo "scale=${scaleForCalculations}; ${dropletRadius};" | bc -l);/;}" system/setFieldsDict

            if [[ ${builtMeshSwitcher} -eq 0 ]]
            then
                let "builtMeshSwitcher += 1"

                PREVIOUS_SUB_CASE_DIR=${SUB_CASE_DIR}

                ###############################
                ##First direction calculation##
                ###############################

                #Preparation STL-files for snappyHexMesh
#                 cd firstDirection
                mv ../*.unv .

                #Preparation mesh for SMESH
                runApplication ideasUnvToFoam porousCell.unv
                runApplication checkMesh -allTopology -allGeometry
                mv log.checkMesh log.checkMesh.originMesh

                #Run the polyDualMesh utility for creation dual mesh
                runApplication polyDualMesh -overwrite -concaveMultiCells 180
                runApplication checkMesh -allTopology -allGeometry
                mv log.checkMesh log.checkMesh.polyDualMesh

                runApplication transformPoints -scale "(${scaledScaleFactor} ${scaledScaleFactor} ${scaledScaleFactor})"

                #Copy mesh for another directions calculation
#                 cp -r constant/ ../secondDirection/constant/
#                 cp -r constant/ ../thirdDirection/constant/

                #Copy 0 to directory for first direction
#                 cp -r ../../src/firstDirection/0/ .

                cp -r ../../src/firstDirection/0/ .

                #Change the case properties
                sed -r -i "s/.*thetaAngle.*/thetaAngle                      thetaAngle                      [0 0 0 0 0 0 0]     ${theta};/" constant/transportProperties
                sed -r -i "s/.*UTranslationDirectionNbTimes.*/UTranslationDirectionNbTimes    UTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${UTranslationDirectionNbTimes};/" constant/transportProperties
                sed -r -i "s/.*VTranslationDirectionNbTimes.*/VTranslationDirectionNbTimes    VTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${VTranslationDirectionNbTimes};/" constant/transportProperties
                sed -r -i "s/.*WTranslationDirectionNbTimes.*/WTranslationDirectionNbTimes    WTranslationDirectionNbTimes    [0 0 0 0 0 0 0]     ${WTranslationDirectionNbTimes};/" constant/transportProperties
                sed -r -i "s/.*unitCellSize.*/unitCellSize                    unitCellSize                    [0 0 0 0 0 0 0]     ${unitCellSize};/" constant/transportProperties
                sed -r -i "s/.*L.*/L                               L                               [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${unitCellSize}*${UTranslationDirectionNbTimes};" | bc -l);/" constant/transportProperties

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
            else
                cp -r ../../src/firstDirection/0/ .

                #Copy mesh for next calculation
                cp -r ../${PREVIOUS_SUB_CASE_DIR}/constant/polyMesh/ constant/
            fi

            #interFoam calculation
            runApplication setFields
            runApplication decomposePar
#             runParallel interFlow
            runParallel interFoam

            cd ..

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
    
            dropletRadius=$(echo "${dropletRadius} + ${dropletRadiusStep}" | bc -l)
        done

        cd ..

    #     PREVIOUS_CASE_DIR=${CASE_DIR}

        caseExecutionTime=$(echo "$(date +%s.%N) - ${startTime}" | bc)

        echo -e "\nCase execution time:" ${caseExecutionTime} "seconds\n"

#         intersectionParameter=$(echo "scale=${scaleForCalculations}; ${intersectionParameter} + ${intersectionParameterStep};" | bc -l)
#         porosity=$(echo "scale=${scaleForCalculations}; ${porosity} + ${porosityStep};" | bc -l)
    done

    theta=$(echo "scale=${scaleForCalculations}; ${theta} + ${thetaStep};" | bc -l)
done

exit 0
