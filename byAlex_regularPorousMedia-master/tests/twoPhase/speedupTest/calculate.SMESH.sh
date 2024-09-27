#!/bin/bash

. $WM_PROJECT_DIR/bin/tools/RunFunctions

#Variables section
minNumberOfSubdomains=4
maxNumberOfSubdomains=18
numberOfSubdomainsStep=2

numberOfSubdomains=${minNumberOfSubdomains}

work=$(echo "((${maxNumberOfSubdomains} - ${minNumberOfSubdomains}) / ${numberOfSubdomainsStep} + 1);" | bc)

while [[ $(echo "${numberOfSubdomains} <= ${maxNumberOfSubdomains};" | bc) -ne 0 ]]
do
    startTime=$(date +%s.%N)

    let "workNum += 1"

    CASE_DIR=numberOfSubdomains$(echo "${numberOfSubdomains};" | bc)

    echo -e "\n###################################"
    echo -e "Case: "${CASE_DIR}
    echo -e "Calculation progress: " ${workNum} / ${work}
    echo -e "###################################\n"

    mkdir -p ${CASE_DIR}

    cd ${CASE_DIR}

    #Copy sources to directory for first direction
    cp -r ../src/0/ .
    cp -r ../src/constant/ .
    cp -r ../src/system/ .

    sed -r -i "s/numberOfSubdomains.*/numberOfSubdomains  ${numberOfSubdomains};/" system/decomposeParDict

    cp system/controlDict.interFoam system/controlDict
    cp system/fvSchemes.interFoam system/fvSchemes
    cp system/fvSolution.interFoam system/fvSolution

#     cp system/controlDict.interFlow system/controlDict
#     cp system/fvSchemes.interFlow system/fvSchemes
#     cp system/fvSolution.interFlow system/fvSolution

    #interFoam calculation
    runApplication setFields
    runApplication decomposePar
#     runParallel interFlow
    runParallel interFoam

    cd ..

    echo -e "\nCase execution time:" $(echo "$(date +%s.%N) - ${startTime}" | bc) "seconds\n"

    numberOfSubdomains=$(echo "${numberOfSubdomains} + ${numberOfSubdomainsStep};" | bc -l)
done

exit 0
