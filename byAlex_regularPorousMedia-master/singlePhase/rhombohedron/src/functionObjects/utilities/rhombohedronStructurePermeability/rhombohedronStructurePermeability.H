/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::rhombohedronStructurePermeability

Description
    Evaluates and outputs permeability for rhombohedron structure.

    Example of function object specification:
    \verbatim
    permeability1
    {
        type            rhombohedronStructurePermeability;
        libs            ("rhombohedronStructureLibutilityFunctionObjects.so");

        inflowPatches   (\<inflowPatchesSet\>);
        outflowPatches  (\<outflowPatchesSet\>);
    }
    \endverbatim

    Output data is written to the dir postProcessing/permeability1/\<timeDir\>/.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::timeControl

SourceFiles
    rhombohedronStructurePermeability.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_rhombohedronStructurePermeability_H
#define functionObjects_rhombohedronStructurePermeability_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
            Class rhombohedronStructurePermeability Declaration
\*---------------------------------------------------------------------------*/

class rhombohedronStructurePermeability
:
    public fvMeshFunctionObject,
    public logFiles
{
protected:

    // Protected data
    // Read from dictionary

        //- Inflow patches
        labelHashSet inflowPatchSet_;

        //- Outflow Patches
        labelHashSet outflowPatchSet_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


private:

        //- Disallow default bitwise copy construct
        rhombohedronStructurePermeability(const rhombohedronStructurePermeability&);

        //- Disallow default bitwise assignment
        void operator=(const rhombohedronStructurePermeability&);


public:

    //- Runtime type information
    TypeName("rhombohedronStructurePermeability");


    // Constructors

        //- Construct from Time and dictionary
        rhombohedronStructurePermeability
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rhombohedronStructurePermeability();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate and write the permeability
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
