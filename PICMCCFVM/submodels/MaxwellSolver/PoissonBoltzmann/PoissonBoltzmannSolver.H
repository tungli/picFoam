/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2020 picFoam
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
    Foam::PoissonBoltzmannSolver

Description
    Poisson Boltzmann electrons field solver. Solves nonlinear Poisson's
    equation with equilibrium (Boltzmann) electron density.
\*---------------------------------------------------------------------------*/

#ifndef PoissonBoltzmann_H
#define PoissonBoltzmann_H

#include "MaxwellSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class PoissonBoltzmannSolver Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PoissonBoltzmannSolver
:
    public MaxwellSolver<CloudType>
{
private:
    dimensionedScalar ne_Inf;
    dimensionedScalar Te;
    dimensionedScalar phiENewtonTol;
    label maxNewtonIter;

public:

    //- Runtime type information
    TypeName("PoissonBoltzmann");


    // Constructors

        //- Construct from dictionary
        PoissonBoltzmannSolver
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~PoissonBoltzmannSolver();


    // Member Functions

        //- Solve Poisson equation
        virtual void solveFields();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PoissonBoltzmannSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
