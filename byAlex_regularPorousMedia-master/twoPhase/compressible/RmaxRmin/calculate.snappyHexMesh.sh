 #!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
SOLUTION_DIR=$(pwd)

scaleForFolderNames=6
scaleForCalculations=20

channelLength0=$(echo "scale=${scaleForFolderNames}; 20000.0/1.0;" | bc -l)

Rmin0=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)
Rmax0=$(echo "scale=${scaleForFolderNames}; 2000.0/1.0;" | bc -l)
RMinStep=$(echo "scale=${scaleForFolderNames}; 100.0/1.0;" | bc -l)
RMaxStep=$(echo "scale=${scaleForFolderNames}; 20000.0/1.0;" | bc -l)

minDropletRadius=$(echo "scale=${scaleForFolderNames}; 500.0/1.0;" | bc -l)
maxDropletRadius=$(echo "scale=${scaleForFolderNames}; 10000.0/1.0;" | bc -l)
dropletRadiusStep=$(echo "scale=${scaleForFolderNames}; 500.0/1.0;" | bc -l)

scaleSize=$(echo "scale=${scaleForCalculations}; 0.000001/1.0;" | bc -l)

work=($(find . -type d -name "Rmax*"))
workNum=0

#Compile some functions
cp -r src/functionObjects/ functionObjects/
./functionObjects/Allwmake

Rmax=${Rmax0}

while [[ $(echo "scale=${scaleForFolderNames}; ${Rmax} >= ${Rmin0};" | bc -l) -ne 0 ]]
do
    Rmin=${Rmin0}

    while [[ $(echo "scale=${scaleForFolderNames}; ${Rmin} <= ${Rmax};" | bc -l) -ne 0 ]]
    do
        startTime=$(date +%s.%N)

        let "workNum += 1"

        CASE_DIR=Rmax$(echo "scale=${scaleForFolderNames}; ${Rmax};" | bc -l)Rmin$(echo "scale=${scaleForFolderNames}; ${Rmin};" | bc -l)

        echo -e "\n###################################"
        echo -e "Case: "${CASE_DIR}
        echo -e "Calculation progress: " ${workNum} / ${#work[*]}
        echo -e "###################################\n"

        cd ${CASE_DIR}

        dropletRadius=${minDropletRadius}

        builtMeshSwitcher=0

        while [[ $(echo "${dropletRadius} <= ${maxDropletRadius};" | bc -l) -ne 0 ]]
        do
            SUB_CASE_DIR=dropletRadius$(echo "scale=${scaleForFolderNames}; ${dropletRadius};" | bc -l)

            mkdir -p ${SUB_CASE_DIR}

            cd ${SUB_CASE_DIR}

            #Copy sources to directory
            cp -r ../../src/0/ .
            cp -r ../../src/constant/ .
            cp -r ../../src/system/ .

            sed -r -i "s/.*channelLength.*/channelLength   channelLength   [0 1 0 0 0 0 0]     $(echo "scale=${scaleForCalculations}; ${scaleSize}*${channelLength0}/1.0;" | bc -l);/" constant/thermophysicalProperties
            sed -r -i "/sphereToCell/{N;N;N;s/centre.*/centre                  (0 0 $(echo "scale=${scaleForCalculations}; 0.5*${scaleSize}*${channelLength0};" | bc -l));/;/.*centre.*/a\
            radius                  $(echo "scale=${scaleForCalculations}; ${scaleSize}*${dropletRadius};" | bc -l);\n
            }" system/setFieldsDict

            if [[ ${builtMeshSwitcher} -eq 0 ]]
            then
                let "builtMeshSwitcher += 1"

                PREVIOUS_SUB_CASE_DIR=${SUB_CASE_DIR}

                #Build mesh by snappyHexMesh
                mv ../*.unv .
                mv ../*.stl constant/triSurface/
                mv ../locationInMesh system/
                cp system/snappyHexMeshSrc/* system/
                cp system/decomposeParDict.snappyHexMesh system/decomposeParDict
                #runApplication surfaceOrient constant/triSurface/planeChannel3D.stl constant/triSurface/planeChannel3D.stl "(1e10 1e10 1e10)"
                runApplication surfaceCheck constant/triSurface/planeChannel3D.stl
                runApplication ideasUnvToFoam backgroundMesh.unv
                runApplication surfaceFeatureExtract
                runApplication surfaceFeatureConvert constant/triSurface/*.eMesh constant/triSurface/edges.vtk
                runApplication decomposePar -copyZero
                runParallel snappyHexMesh -overwrite
                runApplication reconstructParMesh -constant
                rm -r processor* log.decomposePar
                runApplication extrudeMesh
                runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"
                runApplication checkMesh -allGeometry

        #         #Build mesh by SMESH
        #         cd ${CASE_DIR}
        #         runApplication ideasUnvToFoam channel.unv
        #         runApplication checkMesh -allGeometry
        #         runApplication transformPoints -scale "(${scaleSize} ${scaleSize} ${scaleSize})"

                #Assign the translation vectors for cyclicAMI
                #Assign the translation vector for inlet/outlet boundary
                Tx=0
                Ty=0
                Tz=$(echo "scale=${scaleForCalculations}; ${scaleSize}*${channelLength0};" | bc -l)

                sed -r -i "/cyclicInlet/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
                neighbourPatch cyclicOutlet;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
                }" system/createPatchDict

                Tx=0
                Ty=0
                Tz=$(echo "scale=${scaleForCalculations}; -${scaleSize}*${channelLength0};" | bc -l)

                sed -r -i "/cyclicOutlet/{N;N;N;N;s/type.*patch/type cyclicAMI/;/.*type.*/a\
                neighbourPatch cyclicInlet;\ntransform translational;\nseparationVector (${Tx} ${Ty} ${Tz});\nmatchTolerance 1e-4;
                }" system/createPatchDict

                #Run the createPatch utility for creation cyclicAMI
                runApplication createPatch -overwrite

                #Run the topoSet utility for creation AMI-groups
                runApplication topoSet

            else
                #Copy mesh for next calculation
                cp -r ../${PREVIOUS_SUB_CASE_DIR}/constant/polyMesh/ constant/
            fi
            
            cp system/decomposeParDict.compressibleInterFoam system/decomposeParDict

            #compressibleInterFoam calculation
            runApplication setFields
            runApplication decomposePar
            runParallel compressibleInterFoam

            cd ..

            dropletRadius=$(echo "${dropletRadius} + ${dropletRadiusStep}" | bc -l)
        done

        cd ..

        caseExecutionTime=$(echo "$(date +%s.%N) - ${startTime}" | bc)

        echo -e "\nCase execution time:" ${caseExecutionTime} "seconds\n"

        Rmin=$(echo "${Rmin} + ${RMinStep}" | bc -l)
    done

    Rmax=$(echo "${Rmax} - ${RMaxStep}" | bc -l)
done

exit 0
