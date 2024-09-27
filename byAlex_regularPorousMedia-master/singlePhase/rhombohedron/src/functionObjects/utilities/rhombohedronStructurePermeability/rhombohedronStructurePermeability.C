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

#include "rhombohedronStructurePermeability.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rhombohedronStructurePermeability, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        rhombohedronStructurePermeability,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::rhombohedronStructurePermeability::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "rhombohedronStructurePermeability");
        writeCommented(file(), "Time");
        writeTabbed(file(), "rhombohedronStructurePermeability");
        file() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rhombohedronStructurePermeability::rhombohedronStructurePermeability
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    inflowPatchSet_(),
    outflowPatchSet_()
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rhombohedronStructurePermeability::~rhombohedronStructurePermeability()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rhombohedronStructurePermeability::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    inflowPatchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("inflowPatches", wordReList()))
        );

    outflowPatchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("outflowPatches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (inflowPatchSet_.empty() || outflowPatchSet_.empty())
    {
        if (Pstream::master())
        {
            FatalErrorInFunction
                << "Unable to find inflowPatches/outflowPatches"
                << exit(FatalError);
        }
    }

    Info<< "    inflowPatches =";

    forAllConstIter(labelHashSet, inflowPatchSet_, iter)
    {
        label patchi = iter.key();

        Info<< " " << pbm[patchi].name();
    }

    Info<< nl;

    Info<< "    outflowPatches =";

    forAllConstIter(labelHashSet, outflowPatchSet_, iter)
    {
        label patchi = iter.key();

        Info<< " " << pbm[patchi].name();
    }

    Info<< nl << nl;

    return true;
}


bool Foam::functionObjects::rhombohedronStructurePermeability::execute()
{
    return true;
}


bool Foam::functionObjects::rhombohedronStructurePermeability::write()
{
    const dictionary&
        caseProperties = mesh_.lookupObject<IOdictionary>("transportProperties");

    const volScalarField&
        p = mesh_.lookupObject<volScalarField>("p");

    const surfaceScalarField&
        phid = mesh_.lookupObject<surfaceScalarField>("phi");

    dimensionedScalar
        rho(caseProperties.lookup("rho")),
        nu(caseProperties.lookup("nu")),
        thetaAngle(caseProperties.lookup("thetaAngle")),
        L(caseProperties.lookup("L"));

    scalar
        inflowPatchesVolumetricFlow = 0,
        outflowPatchesVolumetricFlow = 0,
        inflowPatchesAveragePressure = 0,
        outflowPatchesAveragePressure = 0;

    forAllConstIter(labelHashSet, inflowPatchSet_, iter)
    {
        label patchi = iter.key();

        inflowPatchesVolumetricFlow += gSum(phid.boundaryField()[patchi]);
        inflowPatchesAveragePressure +=
            gSum
            (
                p.boundaryField()[patchi]
              * mesh_.magSf().boundaryField()[patchi]
            )
          / gSum
            (
                mesh_.magSf().boundaryField()[patchi]
            );
    }

    forAllConstIter(labelHashSet, outflowPatchSet_, iter)
    {
        label patchi = iter.key();

        outflowPatchesVolumetricFlow += gSum(phid.boundaryField()[patchi]);
        outflowPatchesAveragePressure +=
            gSum
            (
                p.boundaryField()[patchi]
              * mesh_.magSf().boundaryField()[patchi]
            )
          / gSum
            (
                mesh_.magSf().boundaryField()[patchi]
            );
    }

    scalar
        alphaAngle =
            Foam::acos
            (
                Foam::cos(scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
              / Foam::cos(0.5 * scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
            ),
        permeability =
            fabs
            (
                nu.value()
                * rho.value()
                * inflowPatchesVolumetricFlow
                * Foam::sin(alphaAngle)
                / (
                    rho.value()
                    * (outflowPatchesAveragePressure - inflowPatchesAveragePressure)
                    * L.value()
                    * Foam::sin(scalar(thetaAngle.value()) * constant::mathematical::pi / 180.0)
                )
            );

    Info<< "rho * inflowPatchesAveragePressure = "
        << rho.value() * inflowPatchesAveragePressure << " [Pa]" << nl
        << "rho * outflowPatchesAveragePressure = "
        << rho.value() * outflowPatchesAveragePressure << " [Pa]" << nl
        << "permeability = " << permeability << " [m^2]" << nl;

    logFiles::write();

    if (Pstream::master())
    {
        writeTime(file());
        file()
            << tab << permeability
            << endl;
    }

    return true;
}


// ************************************************************************* //
