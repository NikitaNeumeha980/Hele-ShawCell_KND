#!/bin/bash

#SBATCH --partition=all
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=12
#SBATCH --output=log.calculate2
#SBATCH --error=log.calculateErr2
#SBATCH --input=/dev/null
#SBATCH --export=ALL

#Variables section
PRECISION=20
SCALE_COEFF=$(echo "scale=${PRECISION}; 0.001;" | bc -l)
PROBLEM_DIR_PATH=$(pwd)
CASE_DIR_PATH=${PROBLEM_DIR_PATH}/case2

cd ${CASE_DIR_PATH}

#Calculation
# foamDictionary  -entry numberOfSubdomains -set 60 system/decomposeParDict
# foamDictionary  -entry method -set scotch system/decomposeParDict

srun \
    setFields -parallel
# runApplication -o\
#     transformPoints "scale=(${SCALE_COEFF} ${SCALE_COEFF} ${SCALE_COEFF})"
# rm log.transformPoints
srun \
    transformPoints -parallel "scale=(${SCALE_COEFF} ${SCALE_COEFF} ${SCALE_COEFF})"

foamDictionary  -entry endTime -set $(echo "2.5" | bc -l) system/controlDict
foamDictionary  -entry deltaT -set $(echo "0.000001" | bc -l) system/controlDict
foamDictionary  -entry writeInterval -set $(echo "0.0005" | bc -l) system/controlDict
foamDictionary  -entry adjustTimeStep -set "yes" system/controlDict
foamDictionary  -entry maxCo -set $(echo "0.1" | bc -l) system/controlDict
foamDictionary  -entry maxAlphaCo -set $(echo "0.1" | bc -l) system/controlDict

# srun \
#     interFoam -parallel
srun \
    interSSFFoam -parallel
# srun \
#     interFSFFoam -parallel
# srun \
#     compressibleInterSSFFoam -parallel

exit 0
