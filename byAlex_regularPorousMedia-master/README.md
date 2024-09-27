# regularPorousMedia
Single phase rhombohedron structure problem.

Work order:
1. Should to execute the cleanSolution.sh - script for cleaning solution directory;
2. Run the generateRhombohedronStructures.SMESH.sh - script for generation a set of structures with the specified parameters. Directories have name like "theta<value of angle>IP<value of intersection>";
3. Run the calculate.SMESH.sh - script for calculation a set of structures by next steps:
    1. Copy source directories;
    2. Convertion of UNV-mesh to openfoam format;
    3. Creation dual-mesh;
    4. Copy mesh for another cases (seconDirection, thirdDirection);
    5. Cheange some text files by sed for setting local parameters for during structure;
    6. Creation AMI-groups;
    7. Calculation.
4. Run the postProcessing.sh - script for postprocessing. This script produce several text files which contain:
    1. permeabilityTensor - permeability tensor;
    2. I1 - first invariant of permiability tensor;
    3. I2 - second invariant of permiability tensor;
    4. I3 - third` invariant of permiability tensor;
    5. permeabilitiesFirstDirection - table of permeabilities for calculation on first direction (through leftSide/rightSide) in angles and intersection parameters;
    6. permeabilitiesSecondDirection - table of permeabilities for calculation on second direction (through oppositeSide/frontSide) in angles and intersection parameters;
    7. permeabilitiesThirdDirection - table of permeabilities for calculation on third direction (though topSide/bottomSide)in angles and intersection parameters;
    8. massFlowsFirstDirection - table of mass flows for calculation on first direction (through leftSide/rightSide) in angles and intersection parameters;
    9. massFlowsSecondDirection - table of mass flows for calculation on second direction (through oppositeSide/frontSide) in angles and intersection parameters;
    10. massFlowsThirdDirection - table of mass flows for calculation on third direction (though topSide/bottomSide) in angles and intersection parameters.
