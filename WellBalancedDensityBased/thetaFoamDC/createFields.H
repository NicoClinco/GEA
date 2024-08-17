Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> thermo
    (
        basicPsiThermo::New(mesh)
    );

    dimensionedScalar Cp("Cp",dimGasConstant,1007.0);
    dimensionedScalar Cv("Cv",dimGasConstant,720.0);
    dimensionedScalar R("R",dimGasConstant,287);
    
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
	IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
   
    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

// Potential temperature perturbation 
   volScalarField thetaP
   (
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T*scalar(0.0),
    Foam::zeroGradientFvPatchField<Foam::scalar>::typeName
   );
   
   Info << "Initializing the density and temperature\n" << endl;
// Initialize the density current bubble:
   #include "InitiateBubble.H"

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Conservative variables

    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    volScalarField Phi
    (
     IOobject
     (
       "Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        dimensionedScalar("g",dimAcceleration,9.81)*mesh.C().component(1)
    );
   // Potential temperature
    volScalarField Theta
    (
        IOobject
        (
            "Theta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thetaP+dimensionedScalar("Theta0",dimTemperature,300.0),
	Foam::zeroGradientFvPatchField<Foam::scalar>::typeName
    );
// Conservative variable for the third equation:
   volScalarField rhoTheta
    (
        IOobject
        (
            "rhoTheta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho*(Theta),
        Foam::zeroGradientFvPatchField<Foam::scalar>::typeName
    );


    // Create numeric flux
    autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New
    (
        rho,
        U,
        rhoTheta
    );
    basicNumericFlux& dbnsFlux = dbnsFluxPtr();

    // Create mass flux alias for easier coupling with other code components
    const surfaceScalarField& phi = dbnsFlux.rhoFlux();
    
    //- Reference to the pressure : 
    volScalarField& p = dbnsFlux.p(); // 
   

    
    
