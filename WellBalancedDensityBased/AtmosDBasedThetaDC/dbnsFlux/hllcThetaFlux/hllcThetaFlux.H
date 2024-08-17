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

Application
    hllcThetaFlux

Description
    HLLC flux class, computes numerical fluxes based on HLLC scheme. Uses
    p,rho,U as inputs variables and return the fluxes rhoU, rhoUU,rhoThetaU

    Toro, E.I. and Spruce, M. and Spears, W.: Restoration of the Contact
    Surface in the HLL-Riemann Solver, Shock Waves, 4:25-34, 1996

Author
    Nicola Clinco

SourceFiles
    hllcThetaFlux.H

\*---------------------------------------------------------------------------*/

#ifndef hllcThetaFlux_H
#define hllcThetaFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hllcFlux Declaration
\*---------------------------------------------------------------------------*/

class hllcThetaFlux
{
public:

    // Constructors

        //- Constructor null
        hllcThetaFlux()
        {}


    // Destructor - default


    // Member Functions

        //- Update fluxes based on current left and right states
	
	//- Potential Temperature constructor
        void evaluateFlux
        (
            scalar& rhoFlux,
            vector& rhoUFlux,
            scalar& rhoThetaFlux,
	    const scalar& rhoLeft,
	    const scalar& rhoRight,
            const scalar& pLeft,
            const scalar& pRight,
            const vector& ULeft,
            const vector& URight,
            const scalar& rhoThetaLeft,
            const scalar& rhoThetaRight,
            const scalar& R,
            const scalar& Cv,
            const vector& Sf,
            const scalar& magSf
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
