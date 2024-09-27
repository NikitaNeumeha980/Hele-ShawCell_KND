#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#For avoid mpi buffer size error
export MPI_BUFFER_SIZE=2000000000

#Variables section
PRECISION=20
SCALE_COEFF=$(echo "scale=${PRECISION}; 0.001;" | bc -l)
PROBLEM_DIR_PATH=$(pwd)
CASE_DIR_PATH=${PROBLEM_DIR_PATH}/case2

mkdir -p ${CASE_DIR_PATH}

cd ${CASE_DIR_PATH}

#Copy sources to case directory
cp -r ${PROBLEM_DIR_PATH}/src/0/ .
cp -r ${PROBLEM_DIR_PATH}/src/constant/ .
cp -r ${PROBLEM_DIR_PATH}/src/system/ .

#Preparation STL-files for snappyHexMesh
mkdir -p constant/triSurface/
cp ${PROBLEM_DIR_PATH}/src/model/*.stl constant/triSurface/

runApplication -o\
    surfaceClean -noClean constant/triSurface/model.stl constant/triSurface/model.stl 0 0
runApplication -o\
    surfaceFeatures -dict system/snappyHexMeshSrc/surfaceFeaturesDict
runApplication -o\
    surfaceFeatureConvert\
        constant/triSurface/model.eMesh\
        constant/triSurface/model.vtk

#Preparation background mesh for snappyHexMesh
runApplication -o\
    blockMesh

#Run snappyHexMesh
foamDictionary  -entry numberOfSubdomains -set 480 system/decomposeParDict
foamDictionary  -entry method -set scotch system/decomposeParDict
# foamDictionary  -entry method -set simple system/decomposeParDict
# foamDictionary  -entry simpleCoeffs -set "{n (4 4 4); delta 0.001;}" system/decomposeParDict

runApplication -o\
    decomposePar -force -copyZero

# runApplication -o\
#     decomposePar -force -noZero

# for directory in processor* ; do
#     echo "${directory}"
#     mkdir -p ${directory}/constant/triSurface
#     ln -sr constant/triSurface/smoothedRockSkeleton.obj -t ${directory}/constant/triSurface/
# done
# 
# runParallel -o\
#     surfaceRedistributePar -keepNonMapped smoothedRockSkeleton.obj follow

exit 0
