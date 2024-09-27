rm  permeabilityTensor\
    I1\
    I2\
    I3\
    permeabilitiesFirstDirection\
    permeabilitiesSecondDirection\
    permeabilitiesThirdDirection\
    massFlowsFirstDirection\
    massFlowsSecondDirection\
    massFlowsThirdDirection

g++ -O3 src/permeabilityTensorCalculation.cpp -o permeabilityTensorCalculation

scaleForFolderNames=6
scaleForCalculations=20

minTheta=$(echo "scale=${scaleForFolderNames}; 60.0/1.0;" | bc -l)
maxTheta=$(echo "scale=${scaleForFolderNames}; 90.0/1.0;" | bc -l)
thetaStep=$(echo "scale=${scaleForFolderNames}; 0.25/1.0;" | bc -l)

minIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.01/1.0;" | bc -l)
#maxIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 1.0 - c(0.5*${maxTheta}/180*4*a(1));" | bc -l)
intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.0025/1.0;" | bc -l)

theta=${minTheta}

while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
do
    while [[ $(echo "scale=${scaleForFolderNames}; ${intersectionParameter} <= 0.99 - c(0.5*${theta}/180*4*a(1));" | bc -l) -ne 0 ]]
    do
        CASE_DIR=theta$(echo "scale=${scaleForFolderNames}; ${theta};" | bc -l)IP$(echo "scale=${scaleForFolderNames}; ${intersectionParameter};" | bc -l)

        #First direction post-processing
        cd ${CASE_DIR}/firstDirection

        echo "${theta} ${intersectionParameter} $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesFirstDirection
        echo "${theta} ${intersectionParameter} $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsFirstDirection
        echo "$(tail -n 4 mAGradPdevNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
        echo "$(tail -n 4 mAGradPdevNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

        cd ..

        #Second direction post-processing
        cd secondDirection

        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesSecondDirection
        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsSecondDirection
        echo "$(tail -n 4 mAGradPdevNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
        echo "$(tail -n 4 mAGradPdevNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

        cd ..

        #Third direction post-processing
        cd thirdDirection

        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(sed '$!d' permeabilityLUFO)" >> ../../permeabilitiesThirdDirection
        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/left.*right.*relativeDiff/p' | sed -r 's/left.*right.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/opposite.*front.*relativeDiff/p' | sed -r 's/opposite.*front.*relativeDiff//' | sed -r 's/^[ \t]*//')" "$(tail -n 3 inputOutputMassFlowLUFO | sed -r -n '/top.*bottom.*relativeDiff/p' | sed -r 's/top.*bottom.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../../massFlowsThirdDirection
        echo "$(tail -n 4 mAGradPdevNuAndAULUFO | sed -n '3p' | sed -r 's/ /\n/g')" >> ../mAGradPdevNu
        echo "$(tail -n 4 mAGradPdevNuAndAULUFO | sed -n '4p' | sed -r 's/ /\n/g')" >> ../AU

        cd ../../

        ./permeabilityTensorCalculation ${CASE_DIR}

        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l)" >> permeabilityTensor
        echo "$(sed -n '1,3p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> permeabilityTensor
        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(sed -n '4p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I1
        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(sed -n '5p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I2
        echo "$(echo "${theta}" | bc -l) $(echo "${intersectionParameter}" | bc -l) $(sed -n '6p' ${CASE_DIR}/permeabilityTensorI1I2I3)" >> I3

        intersectionParameter=$(echo "${intersectionParameter} + ${intersectionParameterStep}" | bc -l)
    done

    theta=$(echo "${theta} + ${thetaStep}" | bc -l)
done
