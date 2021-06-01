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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PoissonBoltzmannSolver.H"
#include "electromagneticConstants.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PoissonBoltzmannSolver<CloudType>::PoissonBoltzmannSolver
(
    const dictionary& dict,
    CloudType& cloud
)
:
    MaxwellSolver<CloudType>(cloud),
    ne_Inf(dimless/dimVolume, readScalar(dict.subDict("PoissonBoltzmann").lookup("ne_Inf"))),
    Te(dimTemperature, readScalar(dict.subDict("PoissonBoltzmann").lookup("Te"))),
    phiENewtonTol(dimMass * dimArea / dimTime / dimCurrent, readScalar(dict.subDict("PoissonBoltzmann").lookup("phiENewtonTol"))),
    maxNewtonIter(readLabel(dict.subDict("PoissonBoltzmann").lookup("maxNewtonIter")))
{
    Info << "\nBolzmann electrons with:" << nl
         << "ne_Inf          " << ne_Inf << nl
         << "Te              " << Te << nl
         << "phiENewtonTol   " << phiENewtonTol << nl
         << "maxNewtonIter   " << maxNewtonIter << nl
         << endl;
    Info << "\nSolving Poisson's equation..." << endl;

    label elTypeId = cloud.electronTypeId();
    volScalarField& phiE(cloud.elpotentialField());
    volScalarField& rhoCharge(cloud.rhoCharge());
    volScalarField chargeWithoutElectrons = rhoCharge - cloud.rhoChargeSpecies()[elTypeId];
    dimensionedScalar e = constant::electromagnetic::e;
    dimensionedScalar eps0 = constant::electromagnetic::epsilon0;

    fvScalarMatrix poissonEqn
    (
        fvm::laplacian(eps0, phiE) == e * ne_Inf - chargeWithoutElectrons
    );
    poissonEqn.solve();


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PoissonBoltzmannSolver<CloudType>::~PoissonBoltzmannSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PoissonBoltzmannSolver<CloudType>::solveFields()
{
    /* This function is called after `CloudType::calculateFields()` and before data write.
     * Therefore this method can overwrite the density from particle data, no problemo.
     */ 

    CloudType& cloud(this->owner());

    volVectorField& E(cloud.electricField());
    volScalarField& phiE(cloud.elpotentialField());
    volScalarField& rhoCharge(cloud.rhoCharge());

    dimensionedScalar k = constant::physicoChemical::k;
    dimensionedScalar e = constant::electromagnetic::e;
    dimensionedScalar eps0 = constant::electromagnetic::epsilon0;
    dimensionedScalar lambda_D2 = eps0 * k * Te / e / e / ne_Inf;
    dimensionedScalar ekT = e / k / Te;

    label elTypeId = cloud.electronTypeId();
    volScalarField chargeWithoutElectrons = rhoCharge - cloud.rhoChargeSpecies()[elTypeId];

    Info<< "\nSolving Poisson-Boltzmann eq.\n" << endl;

    volScalarField phiEPrevIter(phiE);
    dimensionedScalar phiEMaxDiff(phiE.dimensions(), 1.0);
    label iterN = 0;
    while ((iterN < maxNewtonIter) && (phiEMaxDiff > phiENewtonTol)) {
        iterN += 1;
        Info << "max(ekt * phiE)   " << Foam::max(ekT * phiE) << endl;
        Info << "max(ekt * phiPrevIterE)   " << Foam::max(ekT * phiE) << endl;
        // while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pbEqn
            (
                fvm::laplacian(eps0, phiE) 
                - fvm::Sp(eps0 / lambda_D2 * Foam::exp( ekT * phiE ),  phiE)
                + chargeWithoutElectrons
                - (e * ne_Inf - eps0 * phiEPrevIter / lambda_D2) * Foam::exp( ekT * phiEPrevIter )
            );
            pbEqn.solve();
            phiE.correctBoundaryConditions();
        }

        phiEMaxDiff = Foam::sqrt(Foam::max(Foam::pow(phiE - phiEPrevIter, 2)));
        phiEPrevIter = phiE;
    }
    Info << "\nFinished in " << iterN << " (Newton) iterations." << endl;
    Info << "Last max. potential difference: " << phiEMaxDiff.value() << endl;

    //Update the electric field
    E = -fvc::grad(phiE);
    E.correctBoundaryConditions();
    cloud.eFieldWeighting().update();//FieldWeighting
}


// ************************************************************************* //
