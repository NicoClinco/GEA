/*
  Script used for initializing the density current benchmark.
*/

scalar theta0(300.0);
scalar Pg(100000.0);
scalar g_(9.81);

scalar CvOR = 2.5087;
scalar CpTheta = Cp.value()*theta0;
scalar RTheta = R.value()*theta0;

vectorField VV = mesh.C();


// Initialize the hydrostatic field:
forAll(T,cI)
{
  // Density
  scalar Par = 1.0 - (g_*VV[cI].y())/CpTheta;
  scalar PgrRTheta = Pg/RTheta;
  rho[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  // Temperature 
 T[cI] = theta0 - (g_*VV[cI].y())/Cp.value(); 
  // The pressure is set into the class numericFlux,
  // calculated from the density and temperature
}

scalar rx = 4000.0;
scalar ry = 2000.0;
scalar pi = 3.14159;
forAll(T,cI)
{
  vector pos_ = mesh.C()[cI];
  scalar r = Foam::sqrt(Foam::sqr(pos_.x()/rx) + Foam::sqr((pos_.y()-3000)/ry));
  if(r<=1.001)
  {
  // Density 
  scalar thetacI = theta0-(15.0/2)*(1.0+Foam::cos(pi*r));
  scalar Par = 1.0 - (g_*pos_.y())/(Cp.value()*thetacI);
  scalar PgrRTheta = Pg/(R.value()*thetacI);

  rho[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  // Temperature 
  T[cI] = thetacI-g_*pos_.y()/Cp.value();
  thetaP[cI] = thetacI-theta0;
  }
  // The pressure is set into the class numericFlux  
}

T.correctBoundaryConditions();
thetaP.correctBoundaryConditions();

Info << "Done" << endl;
