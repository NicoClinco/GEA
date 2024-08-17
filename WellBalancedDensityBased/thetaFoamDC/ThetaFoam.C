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
    ThetaFoam

Description
    Density-based compressible explicit time-marching for
    atmospheric flows:
   
    The solver use 3 different set of equations:
    \rho
    \rho u_i
    \rho \theta


NOTE: 
 -This solver is used only for the density-current benchmark test-case.
 -The user should have in the run-directory for the test-case the base file
  that is the temperature and the density
 
Author
    Nicola Clinco

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "hllcThetaFlux.H"
#include "hllcAusmThetaFlux.H"
#include "MDLimiter.H"
#include "firstOrderLimiter.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"
#include "numericFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createTimeControls.H"

 
// File for post-processing about the steady-state
//ofstream FilePostVelocity;
//FilePostVelocity.open("FilePostVel.txt",std::ios::app);
//FilePostVelocity << "t" << "\t" << "Uymax" << "\t" << "Source-flux" << "\n";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Runge-Kutta coefficient
    scalarList beta(4);
    beta[0] = 0.1100;
    beta[1] = 0.2766;
    beta[2] = 0.5000;
    beta[3] = 1.0000;
    dimensionedVector gDim("g",dimAcceleration,vector(0,9.81,0));
    dimensionedScalar mu("mu",Foam::dimensionSet(1,-1,-1,0,0,0,0),scalar(75.0));
    // Switch off solver messages
    lduMatrix::debug = 0;

    while (runTime.run())
    {
     #include "readTimeControls.H"
   // #       include "readFieldBounds.H"
     #include "compressibleCourantNo.H"
     #include "setDeltaT.H"

        runTime++;

        Info<< "\n Time = " << runTime.value() << endl;

        // Low storage Runge-Kutta time integration
        forAll (beta, i)
        {
            // Solve the approximate Riemann problem for this time step
	    // Evaluate the source term

	    
	    dbnsFlux.EvalMomSource();
	    dbnsFlux.EvaluateLimiter();
            dbnsFlux.computeFlux();

	    
	    
	    //FilePostVelocity << runTime.value() << "\t" << Foam::max(mag(U.component(1))).value() << "\t" << 
	    //Foam::max(mag(fvc::div(dbnsFlux.rhoUFlux()) - dbnsFlux.MomentumSource())).value() << "\n";
	   

            // Time integration
	    //Mass
            solve
            (
              1.0/beta[i]*fvm::ddt(rho)
              + fvc::div(dbnsFlux.rhoFlux())
            );
	    //Momentum
            solve
            (
              1.0/beta[i]*fvm::ddt(rhoU)
              + fvc::div(dbnsFlux.rhoUFlux()) == dbnsFlux.MomentumSource() + fvc::laplacian(mu,U) 
            );
	    //Potential temperature
            solve
            (
              1.0/beta[i]*fvm::ddt(rhoTheta)
              + fvc::div(dbnsFlux.rhoThetaFlux())  == fvc::laplacian(mu,Theta)
            );

            #include "updateFields.H"
	   
        }
	
        runTime.write();
	
        Info<< "    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
