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
#include "volFields.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::actuatorDiskSource::addActuatorDiskSource
(
    fvMatrix<vector>& eqn,
    const labelList& cells,
    const scalarField& Vcells,
    const vectorField& Ccells,
    const RhoFieldType& rho,
    const vectorField& U
)
{

    scalar radialDist;
    vector axialVector;
    vector circVector;

    vector totalForce(0.0, 0.0, 0.0);
    scalar totalTorque = 0.0;

    vector meanVel(0.0, 0.0, 0.0);

    vector unitForce(-10, 0, 0);

    scalar diskVolume = 0;

    vectorField& Usource = eqn.source();
    // Zero out force field
    forceField_ *= 0;

    forAll(cells, i)
    {
        if (pointIsInDisk(Ccells[cells[i]], axialVector, circVector))
        {
            radialDist = mag(circVector);
            vector axialForce = diskDir_*calcAxialForce(radialDist)/rho_;

            forceField_[cells[i]] += axialForce;
            Usource[cells[i]] += axialForce*Vcells[cells[i]];

            // compute the total force added to the actuator disk, this is just for control
            totalForce += axialForce*Vcells[cells[i]];

            vector circForce = -circVector/radialDist*calcCircForce(radialDist)/rho_;
            forceField_[cells[i]] += circForce;
            Usource[cells[i]] += circForce*Vcells[cells[i]];

            totalTorque += (calcCircForce(radialDist)/rho_)*radialDist*Vcells[cells[i]];
            diskVolume += Vcells[cells[i]];

            // compute the mean velocity in the actuator disk
            meanVel += U[cells[i]]*Vcells[cells[i]];
        }  
    }

    // Add source to eqn
    //eqn += forceField_;
    
    reduce(meanVel, sumOp<vector>());
    reduce(diskVolume, sumOp<scalar>());
    reduce(totalForce, sumOp<vector>());
    reduce(totalTorque, sumOp<scalar>());

    meanVel = meanVel/diskVolume;

    Info << "Total axial force: " << totalForce << "\n";
    Info << "Total torque: " << totalTorque << "\n";
    Info << "Total disk volume: " << diskVolume << "\n";
    Info << "Mean disk speed: " << meanVel << "\n";

    if (mesh_.time().writeTime())
    {
        forceField_.write();
    }
}


// ************************************************************************* //
