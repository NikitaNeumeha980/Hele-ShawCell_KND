 #!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
scaleForFolderNames=6
scaleForCalculations=20
scaleSize=$(echo "scale=${scaleForCalculations}; 0.000001/1.0;" | bc -l)
caseName=theta135

cd ${caseName}

#Mesh construction
#Preparation STL-files for snappyHexMesh
mkdir -p constant/triSurface/
mv *.stl constant/triSurface/

runApplication  surfaceFeatures -dict system/snappyHexMeshSrc/
runApplication  surfaceFeatureConvert constant/triSurface/*.eMesh constant/triSurface/edges.vtk

#Preparation background mesh for snappyHexMesh
runApplication  blockMesh

#Run snappyHexMesh
foamDictionary  -entry numberOfSubdomains -set 4 system/decomposeParDict
foamDictionary  -entry method -set scotch system/decomposeParDict
runApplication  decomposePar -copyZero
runParallel     snappyHexMesh -overwrite -dict system/snappyHexMeshSrc/
runParallel     extrudeMesh
runParallel     checkMesh -allTopology -allGeometry
runParallel     foamDictionary constant/polyMesh/boundary -entry entry0/back/type -set empty
mv log.foamDictionary log.foamDictionary1
runParallel     foamDictionary constant/polyMesh/boundary -entry entry0/front/type -set empty
mv log.foamDictionary log.foamDictionary2

#Run the topoSet utility for creation viscous sub-layer
runParallel     topoSet

#Relaxation step
cp\
    constant/src/transportProperties.0\
    constant/transportProperties

runParallel     setFields
runParallel     transformPoints "scale=(${scaleSize} ${scaleSize} ${scaleSize})"
foamDictionary  system/controlDict -set\
    "
        application=interFoam,
        endTime=$(echo "$(foamListTimes -processor -latestTime -withZero | tail -1) + 0.5" | bc -l),
        deltaT=0.000001,
        writeControl=adjustableRunTime,
        writeInterval=0.1,
        purgeWrite=0
    "
foamDictionary  -entry ddtSchemes/default -set Euler system/fvSchemes
foamDictionary  -entry relaxationFactors/fields/p -set 1 system/fvSolution
foamDictionary  -entry relaxationFactors/equations/U -set 1 system/fvSolution
runParallel     interFoam
mv log.interFoam log.interFoam.relaxationStep

#Calculation step
cp\
    constant/src/transportProperties.1\
    constant/transportProperties

foamDictionary  system/controlDict -set\
    "
        endTime=$(echo "$(foamListTimes -processor -latestTime -withZero | tail -1) + 10.0" | bc -l),
        writeInterval=0.01
    "
# runParallel     foamDictionary -entry boundaryField/left/type -set "zeroGradient" $(foamListTimes -processor -latestTime -withZero | tail -1)/U
# mv log.foamDictionary log.foamDictionary3
runParallel     foamDictionary -entry boundaryField/right -set\
    "{  
        type                flowRateInletVelocity;
        volumetricFlowRate  2e-12;
        extrapolateProfile  yes;
        value               uniform (0 0 0);
    }
    "\
    $(foamListTimes -processor -latestTime | tail -1)/U
mv log.foamDictionary log.foamDictionary4
# runParallel     foamDictionary -entry boundaryField/left/type -set "zeroGradient" $(foamListTimes -processor -latestTime | tail -1)/p_rgh
# mv log.foamDictionary log.foamDictionary5
runParallel     foamDictionary -entry boundaryField/right/type -set "zeroGradient" $(foamListTimes -processor -latestTime | tail -1)/p_rgh
mv log.foamDictionary log.foamDictionary6
# runParallel     foamDictionary -entry boundaryField/left -set\
#     "{
#         type            pressureNormalInletOutletVelocity;
#         phi             phi;
#         rho             rho;
#         value           uniform (0 0 0);
#     }"\
#     $(foamListTimes -processor -latestTime -withZero | tail -1)/U
# mv log.foamDictionary log.foamDictionary3
# runParallel     foamDictionary -entry boundaryField/right -set\
#     "{
#         type            pressureNormalInletOutletVelocity;
#         phi             phi;
#         rho             rho;
#         value           uniform (0 0 0);
#     }"\
#     $(foamListTimes -processor -latestTime -withZero | tail -1)/U
# mv log.foamDictionary log.foamDictionary4
# runParallel     foamDictionary -entry boundaryField/left -set\
#     "{
#         type            totalPressure;
#         p0              uniform 0;
#     }"\
#     $(foamListTimes -processor -latestTime -withZero | tail -1)/p_rgh
# mv log.foamDictionary log.foamDictionary5
# runParallel     foamDictionary -entry boundaryField/right -set\
#     "{
#         type            totalPressure;
#         p0              uniform 300;
#     }"\
#     $(foamListTimes -processor -latestTime -withZero | tail -1)/p_rgh
# mv log.foamDictionary log.foamDictionary6
runParallel     foamDictionary  $(foamListTimes -processor -latestTime | tail -1)/uniform/time -set\
    "
        deltaT=0.00000001,
        deltaT0=0.00000001
    "
mv log.foamDictionary log.foamDictionary7
# runParallel     foamDictionary -entry "deltaT" -set "0.00000001" $(foamListTimes -processor -latestTime | tail -1)/uniform/time
# mv log.foamDictionary log.foamDictionary7
# runParallel     foamDictionary -entry "deltaT0" -set "0.00000001" $(foamListTimes -processor -latestTime | tail -1)/uniform/time
# mv log.foamDictionary log.foamDictionary8

runParallel     interFoam

# SOLUTION_DIR=$(pwd)
# 
# scaleForFolderNames=6
# scaleForCalculations=20
# 
# channelLength0=$(echo "scale=${scaleForFolderNames}; 20000.0/1.0;" | bc -l)

# Rmin0=$(echo "scale=${scaleForFolderNames}; 500.0/1.0;" | bc -l)
# Rmax0=$(echo "scale=${scaleForFolderNames}; 2000.0/1.0;" | bc -l)
# RMinStep=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)
# RMaxStep=$(echo "scale=${scaleForFolderNames}; 20000.0/1.0;" | bc -l)
# 
# minDropletRadius=$(echo "scale=${scaleForFolderNames}; 1000.0/1.0;" | bc -l)
# maxDropletRadius=$(echo "scale=${scaleForFolderNames}; 10000.0/1.0;" | bc -l)
# dropletRadiusStep=$(echo "scale=${scaleForFolderNames}; 1000.0/1.0;" | bc -l)
# 
# scaleSize=$(echo "scale=${scaleForCalculations}; 0.000001/1.0;" | bc -l)
# 
# work=($(find . -type d -name "Rmax*"))
# workNum=0
# 
# rhoPhase1=$(echo "scale=${scaleForCalculations}; 900.0/1.0;" | bc -l)
# rhoPhase2=$(echo "scale=${scaleForCalculations}; 1000.0/1.0;" | bc -l)
# sigma=$(echo "scale=${scaleForCalculations}; 0.046/1.0;" | bc -l)
# pRef=$(echo "scale=${scaleForCalculations}; 100000.0/1.0;" | bc -l)
# 
# #Compile some functions and BC
# cp -r src/functionObjects/ functionObjects/
# ./functionObjects/Allwmake
# 
# cp -r src/BCs/adaptiveTotalPressure/ adaptiveTotalPressure/
# cd adaptiveTotalPressure/
# wclean
# wmake
# cd ..
# 
# Rmax=${Rmax0}
# 
# while [[ $(echo "scale=${scaleForFolderNames}; ${Rmax} >= ${Rmin0};" | bc -l) -ne 0 ]]
# do
#     Rmin=${Rmin0}
# 
#     while [[ $(echo "scale=${scaleForFolderNames}; ${Rmin} <= ${Rmax};" | bc -l) -ne 0 ]]
#     do
#         startTime=$(date +%s.%N)
# 
#         let "workNum += 1"
# 
#         CASE_DIR=Rmax$(echo "scale=${scaleForFolderNames}; ${Rmax};" | bc -l)Rmin$(echo "scale=${scaleForFolderNames}; ${Rmin};" | bc -l)
# 
#         echo -e "\n###################################"
#         echo -e "Case: "${CASE_DIR}
#         echo -e "Calculation progress: " ${workNum} / ${#work[*]}
#         echo -e "###################################\n"
# 
#         cd ${CASE_DIR}
# 
#         dropletRadius=${minDropletRadius}
# 
#         builtMeshSwitcher=0
# 
#         while [[ $(echo "${dropletRadius} <= ${maxDropletRadius};" | bc -l) -ne 0 ]]
#         do
#             SUB_CASE_DIR=dropletRadius$(echo "scale=${scaleForFolderNames}; ${dropletRadius};" | bc -l)
# 
#             mkdir -p ${SUB_CASE_DIR}
# 
#             cd ${SUB_CASE_DIR}
# 
#             #Copy sources to directory
#             cp -r ../../src/0/ .
#             cp -r ../../src/constant/ .
#             cp -r ../../src/system/ .
# 
#             criticalVolumetricFlowRate=$(echo "scale=${scaleForCalculations}; 0.25*4.0*a(1)*${scaleSize}*${scaleSize}*${Rmin}*${Rmin}/360.0*sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l)
#             p0High=$(echo "scale=${scaleForCalculations}; 20.0*2.0*${sigma}/(${scaleSize}*${Rmin}) + ${pRef};" | bc -l)
#             p0Low=$(echo "scale=${scaleForCalculations}; 0.1*2.0*${sigma}/(${scaleSize}*${Rmin}) + ${pRef};" | bc -l)
# 
#             foamDictionary -entry "boundaryField.inlet.criticalVolumetricFlowRate" -set "${criticalVolumetricFlowRate}" 0/p_rgh
#             foamDictionary -entry "boundaryField.inlet.tc" -set "$(echo "scale=${scaleForCalculations}; 1.0*${scaleSize}*${channelLength0}/sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l)" 0/p_rgh
#             foamDictionary -entry "boundaryField.inlet.p0High" -set "uniform ${p0High}" 0/p_rgh
#             foamDictionary -entry "boundaryField.inlet.p0Low" -set "uniform ${p0Low}" 0/p_rgh
#             foamDictionary -entry "endTime" -set "$(echo "scale=${scaleForCalculations}; 500.0*${scaleSize}*${channelLength0}/sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l)" system/controlDict
#             foamDictionary -entry "writeInterval" -set "$(echo "scale=${scaleForCalculations}; 0.5*${scaleSize}*${channelLength0}/sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l)" system/controlDict
# #             foamDictionary -entry "channelLength" -set "channelLength [0 1 0 0 0 0 0] $(echo "scale=${scaleForCalculations}; ${scaleSize}*${channelLength0}/1.0;" | bc -l)" constant/transportProperties
# 
# #             sed -r -i "s/.*channelLength.*/channelLength   channelLength   [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${scaleSize}*${channelLength0}/1.0;" | bc -l);/" constant/transportProperties
#             sed -r -i "/sphereToCell/{N;N;N;s/centre.*/centre                  (0 0 $(echo "scale=${scaleForCalculations}; 0.5*${scaleSize}*${channelLength0};" | bc -l));/;/.*centre.*/a\
#             radius                  $(echo "scale=${scaleForCalculations}; ${scaleSize}*${dropletRadius};" | bc -l);\n
#             }" system/setFieldsDict
# 
# #             foamDictionary  constant/polyMesh/boundary -entry entry0.empty_1.type -set empty
#             
#             if [[ ${builtMeshSwitcher} -eq 0 ]]
#             then
#                 let "builtMeshSwitcher += 1"
# 
#                 PREVIOUS_SUB_CASE_DIR=${SUB_CASE_DIR}
# 
#                 #Build mesh by snappyHexMesh
#                 mv ../*.unv .
#                 mv ../*.stl constant/triSurface/
#                 mv ../locationInMesh system/
#                 cp system/snappyHexMeshSrc/* system/
#                 cp system/decomposeParDict.snappyHexMesh system/decomposeParDict
#                 #runApplication surfaceOrient constant/triSurface/planeChannel3D.stl constant/triSurface/planeChannel3D.stl "(1e10 1e10 1e10)"
#                 runApplication  surfaceCheck constant/triSurface/planeChannel3D.stl
#                 runApplication  ideasUnvToFoam backgroundMesh.unv
#                 runApplication  surfaceFeatureExtract
#                 runApplication  surfaceFeatureConvert constant/triSurface/*.eMesh constant/triSurface/edges.vtk
#                 runApplication  decomposePar -copyZero
#                 runParallel     snappyHexMesh -overwrite
#                 runApplication  reconstructParMesh -constant
#                 rm -r processor* log.decomposePar
#                 runApplication  extrudeMesh
# #                 runApplication  transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
#                 runApplication  checkMesh -allGeometry
# 
# #                 #Build mesh by SMESH
# #                 cd ${CASE_DIR}
# #                 runApplication ideasUnvToFoam channel.unv
# #                 runApplication checkMesh -allGeometry
# #                 runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
# 
# #                 #Assign the translation vectors for cyclicAMI
# #                 #Assign the translation vector for inlet/outlet boundary
# #                 Tx=0
# #                 Ty=0
# #                 Tz=$(echo "scale=${scaleForCalculations}; ${scaleSize}*${channelLength0};" | bc -l)
# # 
# #                 sed -r -i "/inlet/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
# #                 neighbourPatch cyclicOutlet;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
# #                 }" system/createPatchDict
# # 
# #                 Tx=0
# #                 Ty=0
# #                 Tz=$(echo "scale=${scaleForCalculations}; -${scaleSize}*${channelLength0};" | bc -l)
# # 
# #                 sed -r -i "/outlet/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
# #                 neighbourPatch cyclicInlet;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
# #                 }" system/createPatchDict
# 
#                 #Run the createPatch utility for creation cyclicAMI
#                 runApplication  createPatch -overwrite
# 
#                 #Run the topoSet utility for creation AMI-groups
#                 runApplication  topoSet
#                 
#                 runApplication  transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
# 
#             else
#                 #Copy mesh for next calculation
#                 cp -r ../${PREVIOUS_SUB_CASE_DIR}/constant/polyMesh/ constant/
#             fi
#             
#             cp system/decomposeParDict.interFoam system/decomposeParDict
# 
# #             foamDictionary -entry "internalField" -set "uniform (0 0 $(echo "scale=${scaleForCalculations}; 4.0*sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l))" 0/U
# 
#             runApplication  decomposePar
# 
# #             for directory in processor* ; do
# #                 foamDictionary -entry "boundaryField.inlet.type" -set "fixedValue" ${directory}/0/U
# #                 foamDictionary -entry "boundaryField.inlet" -add "value uniform (0 0 $(echo "scale=${scaleForCalculations}; 4.0*sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l))" ${directory}/0/U
# # #                 foamDictionary -entry "boundaryField.inlet.value" -set "uniform (0 0 $(echo "scale=${scaleForCalculations}; 4.0*sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l))" ${directory}/0/U
# # #                 foamDictionary -entry "boundaryField.outlet.value" -set "uniform (0 0 $(echo "scale=${scaleForCalculations}; 0.5*sqrt(4.0*${sigma}/(${rhoPhase2}*${scaleSize}*${Rmin}));" | bc -l))" ${directory}/0/U
# #             done
# # 
# #             runParallel     potentialFoam
# # 
# #             for directory in processor* ; do
# #                 foamDictionary -entry "boundaryField.inlet.type" -set "zeroGradient" ${directory}/0/U
# # #                 foamDictionary -entry "boundaryField.outlet.type" -set "zeroGradient" ${directory}/0/U
# #             done
# 
#             #interFoam calculation
#             runParallel     setFields
#             runParallel     interFoam
# 
#             cd ..
# 
#             dropletRadius=$(echo "${dropletRadius} + ${dropletRadiusStep}" | bc -l)
#         done
# 
#         cd ..
# 
#         caseExecutionTime=$(echo "$(date +%s.%N) - ${startTime}" | bc)
# 
#         echo -e "\nCase execution time:" ${caseExecutionTime} "seconds\n"
# 
#         Rmin=$(echo "${Rmin} + ${RMinStep}" | bc -l)
#     done
# 
#     Rmax=$(echo "${Rmax} - ${RMaxStep}" | bc -l)
# done

exit 0
