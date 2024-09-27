rm  permeabilities\
    volumetricFlows\
    massFlows

scaleForFolderNames=6

Rmin0=$(echo "scale=${scaleForFolderNames}; 1.0/1.0;" | bc -l)
Rmax0=$(echo "scale=${scaleForFolderNames}; 10.0/1.0;" | bc -l)
Rstep=$(echo "scale=${scaleForFolderNames}; 0.2/1.0;" | bc -l)

Rmax=${Rmin0}

while [[ $(echo "scale=${scaleForFolderNames}; ${Rmax} <= ${Rmax0};" | bc -l) -ne 0 ]]
do
    Rmin=${Rmin0}

    while [[ $(echo "scale=${scaleForFolderNames}; ${Rmin} <= ${Rmax};" | bc -l) -ne 0 ]]
    do
        echo -e "\n###################################"
        echo "Rmax="${Rmax}
        echo "Rmin="${Rmin}
        echo -e "###################################\n"

        CASE_DIR=Rmax$(echo "scale=${scaleForFolderNames}; ${Rmax};" | bc -l)Rmin$(echo "scale=${scaleForFolderNames}; ${Rmin};" | bc -l)

        cd ${CASE_DIR}

        #Post-processing
        echo "${Rmax} ${Rmin} $(sed '$!d' permeabilityLUFO)" >> ../permeabilities
        echo "${Rmax} ${Rmin} $(tail -n 1 postProcessing/flowRatePatch\(name\=inlet\)/0/surfaceFieldValue.dat | sed -r 's/\t//' | sed -r -e 's/[ ]+/ /g')" >> ../volumetricFlows
#         echo "${Rmax} ${Rmin} $(tail -n 1 inputOutputMassFlowLUFO | sed -r -n '/inlet.*outlet.*relativeDiff/p' | sed -r 's/inlet.*outlet.*relativeDiff//' | sed -r 's/^[ \t]*//')" >> ../massFlows

        cd ..

        Rmin=$(echo "${Rmin} + ${Rstep}" | bc -l)
    done

    Rmax=$(echo "${Rmax} + ${Rstep}" | bc -l)
done
