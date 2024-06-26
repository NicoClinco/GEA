/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "basicNumericFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::basicNumericFlux> Foam::basicNumericFlux::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    basicThermo& thermo
)
{
    const dictionary& subDict =
        rho.mesh().schemesDict().subDict("divSchemes").subDict("dbns");

    word name = word(subDict.lookup("flux")) + "Flux"
        + word(subDict.lookup("limiter")) + "Limiter";

    Info<< "Selecting numericFlux " << name << endl;

    BasicStateConstructorTable::iterator cstrIter =
        BasicStateConstructorTablePtr_->find(name);

    if (cstrIter == BasicStateConstructorTablePtr_->end())
    {
        FatalErrorIn("basicNumericFlux::New(const fvMesh&)")
            << "Unknown basicNumericFlux type " << name << nl << nl
            << "Valid basicNumericFlux types are:" << nl
            << BasicStateConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<basicNumericFlux>(cstrIter()(rho, U, T, thermo));
}

// ************************************************************************* //
