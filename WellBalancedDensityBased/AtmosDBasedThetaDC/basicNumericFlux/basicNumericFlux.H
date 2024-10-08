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

Class
    Foam::basicNumericFlux

Description
    Base class for runtime-selectable flux and limiters

SourceFiles
    basicNumericFlux.C
    newBasicNumericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef basicNumericFlux_H
#define basicNumericFlux_H

#include "fvMesh.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicNumericFlux Declaration
\*---------------------------------------------------------------------------*/

class basicNumericFlux
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;


protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        basicNumericFlux(const basicNumericFlux&);


public:

    //- Runtime type information
    TypeName("basicNumericFlux");


    //- Declare run-time constructor selection table
   

    // Construct for theta : Construct from rho,U,
    // and the conservative variable rhoTheta
    declareRunTimeSelectionTable
    (
     autoPtr,
     basicNumericFlux,
     PotentialTemperatureState,
     (
     	const volScalarField& rho,
	const volVectorField& U,
	const volScalarField& rhoTheta
     ),
     (rho, U, rhoTheta)
    );


    // Constructors

        //- Construct from components
        basicNumericFlux
        (
            const fvMesh& mesh
        );


    //- Selectors


    static autoPtr<basicNumericFlux> New
    (
     	const volScalarField& rho,
	const volVectorField& U,
	const volScalarField& rhoTheta
    );



    //- Destructor
    virtual ~basicNumericFlux();


    // Member functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        // Return fluxes

            //- Return density flux
            virtual const surfaceScalarField& rhoFlux() const = 0;

            //- Return velocity flux
            virtual const surfaceVectorField& rhoUFlux() const = 0;

            //- Return energy flux
            virtual const surfaceScalarField& rhoThetaFlux() const = 0;
            
	    //- Return Momentum Source
	    virtual const volVectorField& MomentumSource() const = 0;

	    //- Return Energy Source
	    virtual const volScalarField& EnergySource() const = 0;
	    
	    //- Return Pressure
	    virtual const volScalarField& p() const = 0;
	    //- Return the reference. 
	    virtual volScalarField& p() = 0;
       // Return residuals

            //- Return density equation residual
            virtual tmp<scalarField> rhoResidual() const = 0;

            //- Return momentum equation flux
            virtual tmp<vectorField> rhoUResidual() const = 0;

            //- Return energy equation flux
            virtual tmp<scalarField> rhoThetaResidual() const = 0;

	    virtual void 
	    InitializeHydroStaticField
	    (
	     const volScalarField& rho0,
	     const volScalarField& p0,
	     const volScalarField& T0
	    ) = 0;
		   

        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux() = 0;

	    virtual void EvalMomSource() = 0;

	    virtual void EvalEnergySource() = 0;
	   
	    virtual void EvaluateLimiter() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
