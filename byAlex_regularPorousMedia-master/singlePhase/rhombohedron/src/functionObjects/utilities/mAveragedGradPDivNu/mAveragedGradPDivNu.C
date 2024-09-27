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

#include "mAveragedGradPDivNu.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mAveragedGradPDivNu, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mAveragedGradPDivNu,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mAveragedGradPDivNu::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "mAveragedGradPDivNu");
        writeCommented(file(), "Time");
        writeTabbed(file(), "x");
        writeTabbed(file(), "y");
        writeTabbed(file(), "z");
        file() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mAveragedGradPDivNu::mAveragedGradPDivNu
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

Foam::functionObjects::mAveragedGradPDivNu::~mAveragedGradPDivNu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mAveragedGradPDivNu::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::mAveragedGradPDivNu::execute()
{
    return true;
}


bool Foam::functionObjects::mAveragedGradPDivNu::write()
{
    const dictionary&
        caseProperties = mesh_.lookupObject<IOdictionary>("transportProperties");

    const volScalarField&
        p = mesh_.lookupObject<volScalarField>("p");

    const volVectorField
        gradP = fvc::grad(p);

    dimensionedScalar
        nu(caseProperties.lookup("nu")),
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
        mVolumeAveragedGradPdivNu =
          - gSum
            (
                (gradP.internalField() * mesh_.V())()
            )
          / (
                scalar(unitCellSize.value()) * scalar(UTranslationDirectionNbTimes.value())
              * scalar(unitCellSize.value()) * scalar(VTranslationDirectionNbTimes.value())
              * scalar(unitCellSize.value()) * scalar(WTranslationDirectionNbTimes.value())
              * Foam::sin(alphaAngle)
              * Foam::sin(scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
            )
          / scalar(nu.value());

    Info<< "averagedGradP = ("
        <<  mVolumeAveragedGradPdivNu.component(0) << " "
        <<  mVolumeAveragedGradPdivNu.component(1) << " "
        <<  mVolumeAveragedGradPdivNu.component(2) << " "
        << ") [m/s]" << nl;

    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());
        file()
            << tab << mVolumeAveragedGradPdivNu.component(0)
            << tab << mVolumeAveragedGradPdivNu.component(1)
            << tab << mVolumeAveragedGradPdivNu.component(2)
            << endl;
    }

    return true;
}


// ************************************************************************* //
