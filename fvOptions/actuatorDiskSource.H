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

Class
    Foam::fv::actuatorDiskSource

Description
    Actuation disk source

    Constant values for momentum source for actuation disk
    \f[
        T = 2 \rho A U_{o}^2 a (1-a)
    \f]
    and
    \f[
        U_1 = (1 - a)U_{o}
    \f]

    where:
    \vartable
        A   | disk area
        U_o | upstream velocity
        a   | 1 - Cp/Ct
        U_1 | velocity at the disk
    \endvartable

Usage
    Example usage:
    \verbatim
    actuatorDiskSourceCoeffs
    {
        fields          (U);        // names of fields to apply source
        diskDir         (-1 0 0);   // disk direction
        Cp              0.1;        // power coefficient
        Ct              0.5;        // thrust coefficient
        diskArea        5.0;        // disk area
        upstreamPoint   (0 0 0);    // upstream point
    }
    \endverbatim


SourceFiles
    actuatorDiskSource.C
    actuatorDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorDiskSource_H
#define actuatorDiskSource_H

#include "cellSetOption.H"

#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvMatrices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class actuatorDiskSource Declaration
\*---------------------------------------------------------------------------*/

class actuatorDiskSource
:
    public cellSetOption
{

protected:

    // Disk center point
      point diskLoc_;

    // Disk area normal
      vector diskDir_;

    // Interior and exterior radius
      scalar intR_;
      scalar extR_;

    // Thickness
      scalar thickness_;

    // Thrust and torque
      scalar thrust_;
      scalar torque_;

    // Density
      scalar rho_;

    //- Force field from all elements
      volVectorField forceField_;

      static constexpr double PI_ = 3.141592654;


private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addActuatorDiskSource
        (
            fvMatrix<vector>& eqn,
            const labelList& cells,
            const scalarField& V,
            const vectorField& C,
            const RhoFieldType& rho,
            const vectorField& U
        );

        //- Disallow default bitwise copy construct
        actuatorDiskSource(const actuatorDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const actuatorDiskSource&);

        //- Return whether point in the disk, its axial and circumferential vector to the disk center
        bool pointIsInDisk(const point& cellPoint, vector& axialVector, vector& circVector) const;

        //- Calculate axial force
        scalar calcAxialForce(const scalar& radialDist) const;

        //- Calculate circumferential force
        scalar calcCircForce(const scalar& radialDist) const;

public:

    //- Runtime type information
    TypeName("actuatorDiskSource");


    // Constructors

        //- Construct from components
        actuatorDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~actuatorDiskSource()
    {}


    // Member Functions

        // Access


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "actuatorDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
