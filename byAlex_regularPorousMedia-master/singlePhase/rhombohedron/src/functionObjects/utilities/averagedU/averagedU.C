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

\*---------------------------------------------------------------------------*/

#include "averagedU.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(averagedU, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        averagedU,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::averagedU::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "averagedU");
        writeCommented(file(), "Time");
        writeTabbed(file(), "x");
        writeTabbed(file(), "y");
        writeTabbed(file(), "z");
        file() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::averagedU::averagedU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::averagedU::~averagedU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::averagedU::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::averagedU::execute()
{
    return true;
}


bool Foam::functionObjects::averagedU::write()
{
    const dictionary&
        caseProperties = mesh_.lookupObject<IOdictionary>("transportProperties");

    const volVectorField&
        U = mesh_.lookupObject<volVectorField>("U");

    dimensionedScalar
        thetaAngle(caseProperties.lookup("thetaAngle")),
        UTranslationDirectionNbTimes(caseProperties.lookup("UTranslationDirectionNbTimes")),
        VTranslationDirectionNbTimes(caseProperties.lookup("VTranslationDirectionNbTimes")),
        WTranslationDirectionNbTimes(caseProperties.lookup("WTranslationDirectionNbTimes")),
        unitCellSize(caseProperties.lookup("unitCellSize"));

    scalar
        alphaAngle =
            Foam::acos
            (
                Foam::cos(scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
              / Foam::cos(0.5 * scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
            );

    vector
        volumeAveragedU =
            gSum
            (
                (U.internalField() * mesh_.V())()
            )
          / (
                scalar(unitCellSize.value()) * scalar(UTranslationDirectionNbTimes.value())
              * scalar(unitCellSize.value()) * scalar(VTranslationDirectionNbTimes.value())
              * scalar(unitCellSize.value()) * scalar(WTranslationDirectionNbTimes.value())
              * Foam::sin(alphaAngle)
              * Foam::sin(scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
            );

    Info<< "averagedU = ("
        <<  volumeAveragedU.component(0) << " "
        <<  volumeAveragedU.component(1) << " "
        <<  volumeAveragedU.component(2) << " "
        << ") [m/s]" << nl;

    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());
        file()
            << tab << volumeAveragedU.component(0)
            << tab << volumeAveragedU.component(1)
            << tab << volumeAveragedU.component(2)
            << endl;
    }

    return true;
}


// ************************************************************************* //
