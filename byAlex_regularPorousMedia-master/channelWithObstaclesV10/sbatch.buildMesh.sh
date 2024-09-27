#!/bin/bash

#SBATCH --partition=all
#SBATCH --nodes=40
#SBATCH --ntasks-per-node=12
#SBATCH --output=log.buildMesh2
#SBATCH --error=log.buildMeshErr2
#SBATCH --input=/dev/null
#SBATCH --export=ALL

#Variables section
PRECISION=20
SCALE_COEFF=$(echo "scale=${PRECISION}; 0.001;" | bc -l)
PROBLEM_DIR_PATH=$(pwd)
CASE_DIR_PATH=${PROBLEM_DIR_PATH}/case2

cd ${CASE_DIR_PATH}

# Run the parallel MPI executable
srun \
    snappyHexMesh -parallel -overwrite -dict system/snappyHexMeshSrc
srun \
    renumberMesh -parallel -noZero -overwrite
srun \
    patchSummary -parallel
srun \
    checkMesh -parallel -allTopology -allGeometry

exit 0
