#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
PRECISION=20
SCALE_COEFF=$(echo "scale=${PRECISION}; 0.001;" | bc -l)
PROBLEM_DIR_PATH=$(pwd)
CASE_DIR_PATH=${PROBLEM_DIR_PATH}/case1

cd ${CASE_DIR_PATH}

#Calculation
# foamDictionary  -entry numberOfSubdomains -set 60 system/decomposeParDict
# foamDictionary  -entry method -set scotch system/decomposeParDict

runParallel -o\
    setFields
# runApplication -o\
#     transformPoints "scale=(${SCALE_COEFF} ${SCALE_COEFF} ${SCALE_COEFF})"
# rm log.transformPoints
runParallel -o\
    transformPoints "scale=(${SCALE_COEFF} ${SCALE_COEFF} ${SCALE_COEFF})"

# foamDictionary  -entry endTime -set $(echo "0.05" | bc -l) system/controlDict
# foamDictionary  -entry deltaT -set $(echo "0.000001" | bc -l) system/controlDict
# foamDictionary  -entry writeInterval -set $(echo "0.000001" | bc -l) system/controlDict
# foamDictionary  -entry adjustTimeStep -set "yes" system/controlDict
# foamDictionary  -entry maxCo -set $(echo "0.1" | bc -l) system/controlDict
# foamDictionary  -entry maxAlphaCo -set $(echo "0.1" | bc -l) system/controlDict

foamDictionary  -entry endTime -set $(echo "2.5" | bc -l) system/controlDict
foamDictionary  -entry deltaT -set $(echo "0.000001" | bc -l) system/controlDict
foamDictionary  -entry writeInterval -set $(echo "0.0005" | bc -l) system/controlDict
foamDictionary  -entry adjustTimeStep -set "yes" system/controlDict
foamDictionary  -entry maxCo -set $(echo "0.1" | bc -l) system/controlDict
foamDictionary  -entry maxAlphaCo -set $(echo "0.1" | bc -l) system/controlDict

# runParallel     interFoam
runParallel -o\
    interSSFFoam
# runParallel     interFSFFoam
# runParallel     compressibleInterSSFFoam

exit 0
