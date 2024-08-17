/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "AusmThetaFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::AusmThetaFlux::evaluateFlux
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
) const
{
// Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Ratio of specific heat capacities
    const scalar kappa = (R+Cv)/Cv;

    const scalar ThetaLeft = rhoThetaLeft/rhoLeft;
    const scalar ThetaRight = rhoThetaRight/rhoRight;


    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (ULeft & normalVector);
    const scalar qRight = (URight & normalVector);

    // Speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft =
        Foam::sqrt(max(0.0,kappa * pLeft/rhoLeft));

    const scalar aRight =
        Foam::sqrt(max(0.0,kappa * pRight/rhoRight));

    const scalar aTilde = 0.5*(aLeft+aRight); //Average

    scalar m_dot = 0.0;

    
    const scalar Ku    = 0.30;
    const scalar Kp    = 1e-2; 
    const scalar sigma = 1.0;
    const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));

    const scalar MaInf = 0.1;
    const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf)));
    const scalar MaZero    = Foam::sqrt(sqrMaZero);


    //const scalar fa = MaZero*(2.0-MaZero);

    const scalar fa = sqrt(sqr(1.0-sqrMaZero)*sqrMaDash+4.0*sqrMaZero)/(1.0+sqrMaZero);
    const scalar alpha = 3.0/16.0*(-4.0+5.0*sqr(fa));
    const scalar beta  = 1.0/8.0;

    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);


    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);
    
    // Version Ma4:
    const scalar rhoTilde = 0.5*(rhoLeft+rhoRight);
    const scalar MaP = -Kp*max(1.0-sigma*sqrMaDash,0.0)*((pRight-pLeft)/(fa*rhoTilde*sqr(aTilde)));
    
    const scalar Ma4betaPlusLeft = ((magMaRelLeft >= 1.0) ? Ma1PlusLeft :
		    (Ma2PlusLeft*(1-16.0*beta*Ma2MinusLeft )));
    const scalar Ma4betaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight :
		    (Ma2MinusRight*(1+16.0*beta*Ma2PlusRight )));
    
    scalar MaTilde = Ma4betaPlusLeft + Ma4betaMinusRight + MaP;

    m_dot = MaTilde*aTilde*((MaTilde> 0.0) ? rhoLeft : rhoRight);

    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
        (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) 
         -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
        (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)
         +16.0*alpha*MaRelRight*Ma2PlusRight)));

    const scalar pU = -Ku*P5alphaPlusLeft*P5alphaMinusRight*(rhoLeft+rhoRight)
                                 *(fa*aTilde)*(qRight-qLeft);  
 
    scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight + pU;


    if(m_dot>0)
    {
        rhoFlux = m_dot * magSf;
        rhoUFlux = (m_dot * ULeft + pTilde * normalVector) *magSf;
        rhoThetaFlux = (m_dot * ThetaLeft ) *magSf;
    }
    else
    {
        rhoFlux = m_dot * magSf;
        rhoUFlux = (m_dot * URight + pTilde * normalVector) *magSf;
        rhoThetaFlux = (m_dot * ThetaRight ) *magSf;
    }
  
}







