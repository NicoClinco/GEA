{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

    // Update the boundary condition of rhoTheta (zeroGradient)
    rhoTheta.correctBoundaryConditions();
      
    // Update the potential temperature:

    Theta = rhoTheta/rho;
    Theta.correctBoundaryConditions();

    forAll(T,cI)
    {
      T[cI]  = Theta[cI]-g_*VV[cI].y()/Cp.value();
      thetaP[cI] = Theta[cI]-theta0;
    }
    thetaP.correctBoundaryConditions();
    
    p = rho*R*T;
    p.correctBoundaryConditions();
   
}
