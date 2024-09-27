#!/bin/bash

scaleForCalculations=10
pi=$(echo "scale=${scaleForCalculations}; 4*a(1);" | bc -l)

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
    m=$(echo "scale=${scaleForCalculations}; 1 - ${pi}/12*(2 - 3*${IP}*${IP}*(3 - ${IP}) - 3*${IP2}*${IP2}*(3 - ${IP2}))/((1 - ${IP})*(1 - ${IP})*(1 - ${IP})*(1 - c(${theta}))*sqrt(1 + 2*c(${theta})));" | bc -l)
    echo ${m}
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

echo $(IP 0.08 60.0)
echo $(IP 0.08 63.0)
echo $(IP 0.08 66.0)
echo $(IP 0.08 69.0)
echo $(IP 0.08 72.0)
echo $(IP 0.08 75.0)
echo $(IP 0.08 78.0)
echo $(IP 0.08 81.0)
echo $(IP 0.08 84.0)
echo $(IP 0.08 87.0)
echo $(IP 0.08 90.0)
