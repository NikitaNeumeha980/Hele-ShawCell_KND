#!/bin/bash

#Variables and constant section
SALOME_INSTALL_DIR=/home/tb-itam-alexshtil/SALOME-8.5.0-DB08
SOLUTION_DIR=$(pwd)
# SALOME_SCRIPT=${SOLUTION_DIR}/src/makePorousCell.snappyHexMesh.py
SALOME_SCRIPT=${SOLUTION_DIR}/src/makePorousCell.SMESH.py

numThreads=18
threadNum=0
initialSalomePort=2810
salomePort=${initialSalomePort}

scaleForFolderNames=6
scaleForCalculations=20

pi=$(echo "scale=${scaleForCalculations}; 4*a(1);" | bc -l)

minTheta=$(echo "scale=${scaleForFolderNames}; 60.0/1.0;" | bc -l)
maxTheta=$(echo "scale=${scaleForFolderNames}; 90.0/1.0;" | bc -l)
thetaStep=$(echo "scale=${scaleForFolderNames}; 1.0/1.0;" | bc -l)

# minIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)
# maxIntersectionParameter=$(echo "scale=${scaleForFolderNames}; 0.9 - c(0.5*${maxTheta}/180*${pi});" | bc -l)
# intersectionParameterStep=$(echo "scale=${scaleForFolderNames}; 0.1/1.0;" | bc -l)

minPorosity=$(echo "scale=${scaleForFolderNames}; 0.15/1.0;" | bc -l)
maxPorosity=$(echo "scale=${scaleForFolderNames}; 0.2/1.0;" | bc -l)
porosityStep=$(echo "scale=${scaleForFolderNames}; 0.05/1.0;" | bc -l)

# scaleFactor=$(echo "scale=${scaleForCalculations}; 0.001/1.0;" | bc -l)
increasingScaleFactor=$(echo "scale=${scaleForCalculations}; 100.0/1.0;" | bc -l)

UTranslationDirectionNbTimes=1
VTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}
WTranslationDirectionNbTimes=${UTranslationDirectionNbTimes}

# numStructures=$(echo "((${maxTheta} - ${minTheta}) / ${thetaStep} + 1) * ((${maxIntersectionParameter} - ${minIntersectionParameter}) / ${intersectionParameterStep} + 1);" | bc)
numStructures=$(echo "((${maxTheta} - ${minTheta}) / ${thetaStep} + 1) * ((${maxPorosity} - ${minPorosity}) / ${porosityStep} + 1);" | bc)
structureNum=0

#Functions section

function min(){
    if (( $(echo "${1} <= ${2}" | bc -l) )); then
        echo ${1}
    else
        echo ${2}
    fi
}

function arccos()
{
    scale=${scaleForCalculations}
    if (( $(echo "${1} == 0" | bc -l) )); then
        echo "scale=${scaleForCalculations}; 0.5*${pi}" | bc -l
    elif (( $(echo "(-1 <= ${1}) && (${1} < 0)" | bc -l) )); then
        echo "scale=${scaleForCalculations}; ${pi} - a(sqrt((1/(${1}^2)) - 1))" | bc -l
    elif (( $(echo "(0 < ${1}) && (${1} <= 1)" | bc -l) )); then
        echo "scale=${scaleForCalculations}; a(sqrt((1/(${1}^2)) - 1))" | bc -l
    else
        echo "input out of range"
    fi
}

function porosity()
{
    IP=${1}
    theta=$(echo "scale=${scaleForCalculations}; ${2}/180*${pi};" | bc -l)

#     touchIP=$(echo "scale=${scaleForCalculations}; 1 - 0.5/s(${theta}/2);" | bc -l)

    IP2=$(echo "scale=${scaleForCalculations}; 1 - 2*(1 - ${IP})*s(${theta}/2);" | bc -l)

    if (( $(echo "scale=${scaleForCalculations}; ${IP2} > 0;" | bc -l) )); then
        IP2=${IP2}
    else
        IP2=0
    fi
#     IP2=$(echo "scale=${scaleForCalculations}; 1 - 2*(1 - ${IP})*s(${theta}/2);" | bc -l)
#     IP2=$(echo "scale=${scaleForCalculations}; if (${touchIP} > 0) ${IP2} else 0;" | bc -l)
    porosity=$(echo "scale=${scaleForCalculations}; 1 - ${pi}/12*(2 - 3*${IP}*${IP}*(3 - ${IP}) - 3*${IP2}*${IP2}*(3 - ${IP2}))/((1 - ${IP})*(1 - ${IP})*(1 - ${IP})*(1 - c(${theta}))*sqrt(1 + 2*c(${theta})));" | bc -l)
    echo ${porosity}
}

function IP()
{
    #m - first argument, theta - second argument
    porosity=${1}
    theta=$(echo "scale=${scaleForCalculations}; ${2}/180*${pi};" | bc -l)

#     critIP=$(echo "scale=${scaleForCalculations}; 1 - c(${theta}/2);" | bc -l)
#     critPorosity=$(porosity ${critIP} ${2})

    touchIP=$(echo "scale=${scaleForCalculations}; 1 - 0.5/s(${theta}/2);" | bc -l)
    touchPorosity=$(porosity ${touchIP} ${2})
#     IP2=$(echo "scale=${scaleForCalculations}; 1 - 2*(1 - ${IP})*s(${theta}/2);" | bc -l)
#     IP2=$(echo "scale=${scaleForCalculations}; if (${touchIP} > 0) ${IP2} else 0;" | bc -l)

    F=$(echo "scale=${scaleForCalculations}; 12/${pi}*(1 - ${porosity})*(1 - c(${theta}))*sqrt(1 + 2*c(${theta}));" | bc -l)

    if (( $(echo "scale=${scaleForCalculations}; ${touchPorosity} > ${porosity};" | bc -l) )); then
        a3=$(echo "scale=${scaleForCalculations}; ${F} + 24*s(${theta}/2)*s(${theta}/2)*s(${theta}/2) + 3;" | bc -l)
        a2=$(echo "scale=${scaleForCalculations}; (-3*${F} - 72*s(${theta}/2)^3 - 9)/${a3};" | bc -l)
        a1=$(echo "scale=${scaleForCalculations}; (3*${F} + 72*s(${theta}/2)*s(${theta}/2)*s(${theta}/2) - 18*s(${theta}/2))/${a3};" | bc -l)
        a0=$(echo "scale=${scaleForCalculations}; (-${F} - 24*s(${theta}/2)*s(${theta}/2)*s(${theta}/2) + 18*s(${theta}/2) - 4)/${a3};" | bc -l)
    else
        a3=$(echo "scale=${scaleForCalculations}; ${F} + 3;" | bc -l)
        a2=$(echo "scale=${scaleForCalculations}; (-3*${F} - 9)/${a3};" | bc -l)
        a1=$(echo "scale=${scaleForCalculations}; 3*${F}/${a3};" | bc -l)
        a0=$(echo "scale=${scaleForCalculations}; (2 - ${F})/${a3};" | bc -l)
    fi

    Q=$(echo "scale=${scaleForCalculations}; (${a2}*${a2} - 3*${a1})/9;" | bc -l)
    R=$(echo "scale=${scaleForCalculations}; (2*${a2}*${a2}*${a2} - 9*${a2}*${a1} + 27*${a0})/54;" | bc -l)

    S=$(echo "scale=${scaleForCalculations}; ${Q}*${Q}*${Q} - ${R}*${R};" | bc -l)

    if (( $(echo "scale=${scaleForCalculations}; ${S} > 0;" | bc -l) )); then
        #Three real roots
        arg=$(echo "scale=${scaleForCalculations}; ${R}/sqrt(${Q}*${Q}*${Q});" | bc -l)
        phi=$(echo "scale=${scaleForCalculations}; 1/3*$(arccos ${arg});" | bc -l)

        r1=$(echo "scale=${scaleForCalculations}; -2*sqrt(${Q})*c(${phi}) - ${a2}/3;" | bc -l)
        r2=$(echo "scale=${scaleForCalculations}; -2*sqrt(${Q})*c(${phi} + 2/3*${pi}) - ${a2}/3;" | bc -l)
        r3=$(echo "scale=${scaleForCalculations}; -2*sqrt(${Q})*c(${phi} - 2/3*${pi}) - ${a2}/3;" | bc -l)

        echo $(min $(min ${r1} ${r2}) ${r3})
    elif (( $(echo "scale=${scaleForCalculations}; ${S} < 0;" | bc -l) )); then
        #One real root
        if (( $(echo "scale=${scaleForCalculations}; ${Q} > 0;" | bc -l) )); then
            arg=$(echo "scale=${scaleForCalculations}; sqrt(${R}*${R})/sqrt(${Q}*${Q}*${Q});" | bc -l)
            phi=$(echo "scale=${scaleForCalculations}; 1/3*l(${arg} + sqrt(${arg}*${arg} - 1))" | bc -l)            

            r1=$(echo "scale=${scaleForCalculations}; -2*${R}/sqrt(${R}*${R})*sqrt(${Q})*0.5*(e(${phi}) + e(-${phi})) - ${a2}/3;" | bc -l)
            echo ${r1}
        elif (( $(echo "scale=${scaleForCalculations}; ${Q} < 0;" | bc -l) )); then
            arg=$(echo "scale=${scaleForCalculations}; sqrt(${R}*${R})/e(l(sqrt(${Q}*${Q}))*3/2);" | bc -l)
            phi=$(echo "scale=${scaleForCalculations}; 1/3*l(${arg} + sqrt(${arg}*${arg} + 1))" | bc -l)

            r1=$(echo "scale=${scaleForCalculations}; -2*${R}/sqrt(${R}*${R})*sqrt(${Q})*0.5*(e(${phi}) - e(-${phi})) - ${a2}/3;" | bc -l)
            echo ${r1}
        else
            r1=$(echo "scale=${scaleForCalculations}; -e(l(${a0} - ${a2}*${a2}*${a2}/27)*1/3) - ${a2}/3;" | bc -l)
            echo ${r1}
        fi
    else
        #Two real roots
        r1=$(echo "scale=${scaleForCalculations}; -2*${R}/sqrt(${R}*${R})*sqrt(${Q}) - ${a2}/3;" | bc -l)
        r2=$(echo "scale=${scaleForCalculations}; ${R}/sqrt(${R}*${R})*sqrt(${Q}) - ${a2}/3;" | bc -l)
        echo $(min ${r1} ${r2})
    fi
}

function make_porous_cell(){
    local CASE_DIR=theta$(echo "scale=${scaleForFolderNames}; $2;" | bc -l)IP$(echo "scale=${scaleForFolderNames}; $3;" | bc -l)

    #Preparation case directories
    mkdir -p ${CASE_DIR}

    #Execution SALOME's script
    cd ${SALOME_INSTALL_DIR}
        ./salome kill $1
        ./salome start -t -b --port=$1 ${SALOME_SCRIPT} args:${SOLUTION_DIR}/${CASE_DIR},$2,$3,${increasingScaleFactor},${UTranslationDirectionNbTimes},${VTranslationDirectionNbTimes},${WTranslationDirectionNbTimes}
        ./salome kill $1
    cd ${SOLUTION_DIR}
}

theta=${minTheta}

#Intersection parameter based loop
# while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
# do
#     intersectionParameter=${minIntersectionParameter}
# 
#     while [[ $(echo "scale=${scaleForFolderNames}; ${intersectionParameter} <= 0.9 - c(0.5*${theta}/180*${pi});" | bc -l) -ne 0 ]]
#     do
#         make_porous_cell ${salomePort} ${theta} ${intersectionParameter} &
# 
#         let "salomePort += 1"
#         let "threadNum += 1"
#         let "structureNum += 1"
# 
#         if [[ $(echo "${threadNum} % ${numThreads};" | bc) -eq 0 ]]
#         then
#             salomePort=${initialSalomePort}
#             threadNum=0
#             echo -e "\nGeneration progress: " ${structureNum} / ${numStructures} "\n"
#             wait
#         fi
# 
#         intersectionParameter=$(echo "scale=${scaleForCalculations}; ${intersectionParameter} + ${intersectionParameterStep};" | bc -l)
#     done
# 
#     theta=$(echo "scale=${scaleForCalculations}; ${theta} + ${thetaStep};" | bc -l)
# done

#Porosity based loop
while [[ $(echo "scale=${scaleForFolderNames}; ${theta} <= ${maxTheta};" | bc -l) -ne 0 ]]
do
    porosity=${minPorosity}

    while [[ $(echo "scale=${scaleForFolderNames}; ${porosity} <= ${maxPorosity};" | bc -l) -ne 0 ]]
    do
        make_porous_cell ${salomePort} ${theta} $(IP ${porosity} ${theta}) &

        let "salomePort += 1"
        let "threadNum += 1"
        let "structureNum += 1"

        if [[ $(echo "${threadNum} % ${numThreads};" | bc) -eq 0 ]]
        then
            salomePort=${initialSalomePort}
            threadNum=0
            echo -e "\nGeneration progress: " ${structureNum} / ${numStructures} "\n"
            wait
        fi

        porosity=$(echo "scale=${scaleForCalculations}; ${porosity} + ${porosityStep};" | bc -l)
    done

    theta=$(echo "scale=${scaleForCalculations}; ${theta} + ${thetaStep};" | bc -l)
done

wait

echo -e "\nGenerated: " ${structureNum} / ${numStructures}
echo -e "\nGeneration done!\n"

exit 0

