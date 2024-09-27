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
    Foam::functionObjects::mAveragedGradPDivNu

Description
    Evaluates and outputs mAveragedGradPDivNu for rhombohedron structure.

    Example of function object specification:
    \verbatim
    mAveragedGradPDivNu1
    {
        type        mAveragedGradPDivNu;
        libs        ("rhombohedronStructureLibutilityFunctionObjects.so");
    }
    \endverbatim

Output data is written to the dir postProcessing/mAveragedGradPDivNu1/\<timeDir\>/.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::timeControl

SourceFiles
    mAveragedGradPDivNu.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_mAveragedGradPDivNu_H
#define functionObjects_mAveragedGradPDivNu_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class mAveragedGradPDivNu Declaration
\*---------------------------------------------------------------------------*/

class mAveragedGradPDivNu
:
    public fvMeshFunctionObject,
    public logFiles
{
protected:

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


private:

        //- Disallow default bitwise copy construct
        mAveragedGradPDivNu(const mAveragedGradPDivNu&);

        //- Disallow default bitwise assignment
        void operator=(const mAveragedGradPDivNu&);


public:

    //- Runtime type information
    TypeName("mAveragedGradPDivNu");


    // Constructors

        //- Construct from Time and dictionary
        mAveragedGradPDivNu
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~mAveragedGradPDivNu();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate and write the mAveragedGradPDivNu
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
