/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "actuatorDiskSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(actuatorDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        actuatorDiskSource,
        dictionary
    );

}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::actuatorDiskSource::checkData() const
{

    if (mag(extR_) < VSMALL)
    {
        FatalErrorInFunction
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }

}

bool Foam::fv::actuatorDiskSource::pointIsInDisk(const point& cellPoint, vector& axialVector, vector& circVector) const {

    // Check if a given point is located in the actuator disk region.
    vector posVector(cellPoint - diskLoc_);

    axialVector = (posVector & diskDir_)*diskDir_;
    circVector = posVector - axialVector;

    // Check if the point is inside the actuator disk in the axial direction
    if(mag(axialVector)>=thickness_/2.) {
        return false;
    }

    scalar radialDist = mag(circVector);
    // Check if the point is inside the actuator disk in the radial direction
    return (radialDist <= extR_ && radialDist >= intR_);

}


Foam::scalar Foam::fv::actuatorDiskSource::calcAxialForce(const scalar &radialDist) const {
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute the force component in the axial direction. The force is computed from a simple equation
    // resulting in a force that varies with the radial distance.
    // If you have a better model of a rotor, comment the four lines below and add your own calculation
    // of the axial force.
    // Do not forget to also change the calculation of the tangential force (calcCircForce())
    // below.
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    scalar axialForce = 0.0;
    scalar radiusScaled = (radialDist/extR_ - intR_/extR_)/(1.0 - intR_/extR_);
    scalar Ax = (105.0/8.0)*thrust_/(thickness_*PI_*(3.0*intR_+4.0*extR_)*(extR_-intR_));
    axialForce = Ax*radiusScaled*sqrt(1.0 - radiusScaled);

    return axialForce;
}

Foam::scalar Foam::fv::actuatorDiskSource::calcCircForce(const scalar &radialDist) const {
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute the force component in the tangential direction. The force is computed from a simple equation
    // resulting in a force that varies with the radial distance.
    // Change the four lines below if you have a better model.
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    scalar tangentialForce = 0.0;
    scalar radiusScaled = (radialDist/extR_ - intR_/extR_)/(1.0 - intR_/extR_);
    scalar At = (105.0/8.0)*torque_/(thickness_*PI_*extR_*(extR_-intR_)*(3.0*extR_+4.0*intR_));
    tangentialForce = (At*radiusScaled*sqrt(1.0-radiusScaled)/(radiusScaled*(1.0-intR_/extR_)+intR_/extR_));

    return tangentialForce;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::actuatorDiskSource::actuatorDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    diskLoc_(coeffs_.lookup("centerPoint")),
    diskDir_(coeffs_.lookup("diskDir")),
    intR_(readScalar(coeffs_.lookup("interiorRadius"))),
    extR_(readScalar(coeffs_.lookup("exteriorRadius"))),
    thickness_(readScalar(coeffs_.lookup("thickness"))),
    thrust_(readScalar(coeffs_.lookup("thrust"))),
    torque_(readScalar(coeffs_.lookup("torque"))),
    rho_(readScalar(coeffs_.lookup("density"))),
    forceField_
    (
        IOobject
        (
            "force." + name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "force",
            dimForce/dimVolume/dimDensity,
            vector::zero
        )
    )

{
    coeffs_.lookup("fields") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation disk zone: "
        << this->name() << endl;

    if (fabs(mag(diskDir_)-1)>1e-2)
        diskDir_ = diskDir_/mag(diskDir_);

    {
        Info << "Actuator disk values loaded from fvSolution:\n";
        Info << "centerPoint: " << diskLoc_ << "\n";
        Info << "diskDir: " << diskDir_ << "\n";
        Info << "interiorRadius: " << intR_ << "\n";
        Info << "exteriorRadius: " << extR_ << "\n";
        Info << "thrust: " << thrust_ << "\n";
        Info << "torque: " << torque_ << "\n";
        Info << "density: " << rho_ << "\n";
    }

    checkData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::actuatorDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const scalarField& cellsV = mesh_.V();

    const vectorField& cellsC = mesh_.C();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addActuatorDiskSource
        (
            eqn,
            cells_,
            cellsV,
            cellsC,
            geometricOneField(),
            U
        );
    }
}


void Foam::fv::actuatorDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const scalarField& cellsV = mesh_.V();
    const vectorField& cellsC = mesh_.C();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        addActuatorDiskSource
        (
            eqn,
            cells_,
            cellsV,
            cellsC,
            rho,
            U
        );
    }
}


bool Foam::fv::actuatorDiskSource::read(const dictionary& dict)
{

    if (cellSetOption::read(dict))
    {
        coeffs_.readIfPresent("centerPoint", diskLoc_);
        coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("interiorRadius", intR_);
        coeffs_.readIfPresent("exteriorRadius", extR_);
        coeffs_.readIfPresent("thickness", thickness_);
        coeffs_.readIfPresent("thrust", thrust_);
        coeffs_.readIfPresent("torque", torque_);
        coeffs_.readIfPresent("density", rho_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}




// ************************************************************************* //
