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

#include "numericFlux.H"
#include "MDLimiter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  /** @brief Construct from components:
  /* @param[in] rho: The density field
  /* @param[in] U :  The density-velocity field
  /* @param[in] rhoTheta: The density-pot.temperature field
  */
template<class Flux, class Limiter>
Foam::numericFlux<Flux, Limiter>::numericFlux
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& rhoTheta
)
:
    numericFluxBase<Flux>(rho.mesh()),
    rho_(rho),
    U_(U),
    rhoTheta_(rhoTheta),
    p_
    (
      IOobject
       (
	 "p",
	 this->mesh().time().timeName(),
         this->mesh(),
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
       ),
      (scalar(27.7372)*(Foam::pow((rhoTheta_),1.39861)))
    ),
    rhoFlux_
    (
        "phi",    
        (linearInterpolate(rho_*U_) & this->mesh().Sf())
    ),
    rhoUFlux_
    (
      "rhoUFlux",
      rhoFlux_*linearInterpolate(U_)
    ),
    rhoThetaFlux_
    (
        "rhoThetaFlux",  
        (linearInterpolate(rhoTheta_*U_) & this->mesh().Sf())
    ),
   MomentumSource_
   (
       IOobject
       (
         "MomSource",
         this->mesh().time().timeName(),
         this->mesh(),
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
       ),  
       rho_*dimensionedVector("g",dimAcceleration,vector(0,-9.81,0))*scalar(0.0),
       Foam::zeroGradientFvPatchField<Foam::vector>::typeName
   ),
   EnergySource_
   (
      "Esource",
      rho_*U_ & (dimensionedVector("g",dimAcceleration,vector(0,-9.81,0)))
   ),
   Dp_
   (
    "Dp",
    p_*dimensionedVector("z",dimensionSet(0,-1,0,0,0,0,0),vector(0.0,0.0,0.0))
   ),
   Drho_
   (
    "Drho",
    rho_*dimensionedVector("z",dimensionSet(0,-1,0,0,0,0,0),vector(0.0,0.0,0.0))
   )
{
 // Change the dimensions of p_:
 p_.dimensions().reset(dimensionSet(1,-1,-2,0,0,0,0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = this->mesh().owner();
    const unallocLabelList& neighbour = this->mesh().neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = this->mesh().Sf();
    const surfaceScalarField& magSf = this->mesh().magSf();

    const volVectorField& cellCentre = this->mesh().C();
    const surfaceVectorField& faceCentre = this->mesh().Cf();
    
    const volVectorField& Drho = this->Drho_;
    const volVectorField& Dp   = this->Dp_;
    // const volTensorField& DV   = this->DV_;

    const scalar Cv = 720.0;
    const scalar R  = 287.0;
    const scalar Cp = Cv + R;

    // Limiter for velocity
    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();
    // Standard Limiter for velocity:
    MDLimiter<vector, Limiter> vectorULimiter
    (
     this->U_,
     gradU
    );
    const volVectorField& ULimiter = vectorULimiter.phiLimiter();
    
    //Limiter for Potential temperature
    const volScalarField& rhoTheta(rhoTheta_);
    const tmp<volVectorField> tgradrhoTheta = fvc::grad(rhoTheta_);
    const volVectorField& gradrhoTheta = tgradrhoTheta();

    MDLimiter<scalar, Limiter> scalarrhoThetaLimiter
    (
     this->rhoTheta_,
     gradrhoTheta
    );
    const volScalarField& rhoThetaLimiter = scalarrhoThetaLimiter.phiLimiter();



    forAll (owner, faceI)
    {
        const label own = owner[faceI];          // own.
        const label nei = neighbour[faceI];      // nei.

        const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];
        
	scalar deltapRight   = Dp[nei] & (deltaRRight);
	scalar deltapLeft    = Dp[own] & (deltaRLeft);
	scalar deltarhoRight = Drho[nei] & (deltaRRight);
	scalar deltarhoLeft  = Drho[own] & (deltaRLeft);

	scalar rhoFNe = 0.0;
	scalar rhoFOwn = 0.0;
	scalar pFNe = 0.0;
	scalar pFOwn = 0.0;
	
	// Construct Local equilibrum for the face:
	
	GEA::ConstructLocalEquilibrum
	(
	  faceCentre[faceI].y(),
	  cellCentre[nei].y(),
	  cellCentre[own].y(),
	  rho_[own], 
	  rho_[nei],
	  p_[own],
	  p_[nei],
	  rhoFNe,
	  rhoFOwn,
	  pFNe,
	  pFOwn
	);
	
	 pFOwn += deltapLeft;
	 pFNe  += deltapRight;
	 rhoFOwn += deltarhoLeft;
         rhoFNe  += deltarhoRight;

	

	// Evaluate the standard flux.

        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoThetaFlux_[faceI],
            rhoFOwn,
            rhoFNe,
	    pFOwn,
	    pFNe,
            U_[own] + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
            U_[nei] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
            rhoTheta[own]+rhoThetaLimiter[own]*(deltaRLeft & gradrhoTheta[own]),
            rhoTheta[nei]+rhoThetaLimiter[nei]*(deltaRRight & gradrhoTheta[nei]),
            R,
            Cv,
            Sf[faceI],
            magSf[faceI]
        );
    } // end internal loop

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        //const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoThetaFlux = rhoThetaFlux_.boundaryField()[patchi];

        // Patch fields
        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const vectorField& pU = U_.boundaryField()[patchi];
        const scalarField& prhoTheta = rhoTheta_.boundaryField()[patchi];
	const labelUList& pFaceCells = p_.mesh().boundary()[patchi].faceCells();
        // const scalarField& pRho = rho_.boundaryField()[patchi];
       
        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
    	const volVectorField& cellCentre = this->mesh().C();
	const fvsPatchVectorField& pfC = faceCentre.boundaryField()[patchi];

        forAll (pp, facei)
        {
           scalar rhoFNe = 0.0;
           scalar rhoFOwn = 0.0;
           scalar pFNe = 0.0;
           scalar pFOwn = 0.0;

           // Construct Local equilibrum for the face:

           GEA::ConstructLocalEquilibrum
           (
            pfC[facei].y(),
            cellCentre[pFaceCells[facei]].y(),
            cellCentre[pFaceCells[facei]].y(),
            rho_[pFaceCells[facei]],
            rho_[pFaceCells[facei]],
            p_[pFaceCells[facei]],
            p_[pFaceCells[facei]],
            rhoFNe,
            rhoFOwn,
            pFNe,
            pFOwn
          );

         scalar TFown = pFOwn/(R*rhoFOwn);
         scalar TFnei = TFown;
         
	 scalar rhoThetaOwn = rhoFOwn*(TFown+9.81*pfC[facei].y()/Cp);
         scalar rhoThetaNe  = rhoThetaOwn;

          Flux::evaluateFlux
          (
             pRhoFlux[facei],
             pRhoUFlux[facei],
             pRhoThetaFlux[facei],
             rhoFOwn,
	     rhoFNe,	     
             pFOwn,
             pFOwn,
             pU[facei],
             pU[facei],
             rhoThetaOwn,
             rhoThetaNe,
             R,
             Cv,
             pSf[facei],
             pMagSf[facei]
          );
        }// end loop specific boundary
   }// end loop boundaries
}// end method

/**
 * @brief Evaluate the momentum source term.
 */
template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::EvalMomSource()
{
  volVectorField& Mom = this->MomentumSource_;
  Mom *= scalar(0.0);
  // At every timestep we freeze the background-state
  const volScalarField& p0 = this->p_;
  const volScalarField& rho0 = this->rho_;

  GEA::EvalBalancedSource
  (
    Mom,
    p0,
    rho0
  );

}

//Not used in this module
template<class Flux,class Limiter>
void Foam::numericFlux<Flux, Limiter>::EvalEnergySource()
{
  volScalarField& Esource = this->EnergySource_;
  dimensionedVector g("g",dimAcceleration,vector(0,-9.81,0));
  Esource.internalField() = (rho_*(g & U_));
}

/** @brief Evaluate a limiter for p,rho customizable (for fully orthogonal) mesh
 *  such as in the DC benchmark. 
 *  
 * @note
 *   The limiter is defined in the 'ReconstructionUtilities.H'
 * @endnote
 */
template<class Flux,class Limiter>
void Foam::numericFlux<Flux, Limiter>::EvaluateLimiter()
{
  const volScalarField& rho = this->rho_;
  const volScalarField&  p = this->p_;
  volVectorField& Drho = this->Drho_;
  volVectorField& Dp   = this->Dp_;
  
  Drho*=scalar(0.0);
  Dp*=scalar(0.0);   
  GEA::ConstructLimiter
  (
    Drho,
    Dp,
    rho,
    p
  );

}



// ************************************************************************* //
