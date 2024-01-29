Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu1
(
    "nu1",
    dimViscosity,
    transportProperties
);

dimensionedScalar nu2
(
    "nu2",
    dimViscosity,
    transportProperties
);

/*dimensionedScalar DT1
(
    "DT1",
    //dimViscosity,
    transportProperties
);

dimensionedScalar DT2
(
    "DT2",
    //dimViscosity,
    transportProperties
);

dimensionedScalar DTS
(
    "DTS",
    //dimViscosity,
    transportProperties
);*/

dimensionedScalar rhoc1
(
    "rhoc1",
    dimless,
    transportProperties
);

dimensionedScalar rhoc2
(
    "rhoc2",
    dimless,
    transportProperties
);

//dimensionedScalar DT
//(
//    "DT",
//    //dimViscosity,
//    transportProperties
//);


//grey control
dimensionedScalar q
(
    "q",
    dimless,
    transportProperties
);

//grey control
dimensionedScalar one
(
    "one",
    dimless,
    transportProperties
);

//minimum value of interpolation function
dimensionedScalar alphaMin
(
    "alphaMin",
    dimless/dimTime,
    transportProperties
);

//maximum value of interpolation function
dimensionedScalar alphaMax
(
    "alphaMax",
    dimless/dimTime,
    transportProperties
);

//Minimum value of design variable
dimensionedScalar zeroVf
(
    "zeroVf",
    dimless,
    transportProperties
);

//Maximum value of the design variable
dimensionedScalar vfMax
(
    "vfMax",
    dimless,
    transportProperties
);

volScalarField vf1
(
     IOobject
     (
         "vf1",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     mesh
);

volScalarField vf2
(
     IOobject
     (
         "vf2",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     mesh
);

Info<< "Recalculating vf1 and vf2\n" << endl;
const volVectorField& Cell = mesh.C();
forAll(Cell, cellI)
{
    #include "densityProjection.H"
}
vf1.write();
vf2.write();
vf1.correctBoundaryConditions();
vf2.correctBoundaryConditions();
//const volScalarField& vf3 = db().lookupObject<volScalarField>("vf2");
//Info<< vf3 << endl;
//Info<< "Mesh: " << Cell << endl;
Info<< "Projection Completed!" << endl;
//Info<< vf3 << endl;

volScalarField vf1Pen
(
     IOobject
     (
         "vf1Pen",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     vf1*(1 + q)/(vf1 + q)
);

volScalarField vf2Pen
(
     IOobject
     (
         "vf2Pen",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     vf2*(1 + q)/(vf2 + q)
); 

//vf1Pen=vf1*(1 + q)/(vf1 + q);
//vf2Pen=vf2*(1 + q)/(vf2 + q);

/*volScalarField Dvf1PenDvf1
(
     IOobject
     (
         "Dvf1PenDvf1",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     (1 + q)/(vf1 + q) - vf1*(1 + q)/((vf1 + q)*(vf1 + q))
);

volScalarField Dvf2PenDvf2
(
     IOobject
     (
         "Dvf2PenDvf2",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     (1 + q)/(vf2 + q) - vf2*(1 + q)/((vf2 + q)*(vf2 + q))
); */

/*volScalarField gsensVol
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1-vf1,
    zeroGradientFvPatchScalarField::typeName
);*/

volScalarField alpha1
(
     IOobject
     (
         "alpha1",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     //alphaMax + (alphaMin - alphaMax)*(vf1*(1 + q)/(vf1 + q))
     vf1Pen*(vf2Pen*alphaMax + (1 - vf2Pen)*alphaMin) + (1 - vf1Pen)*alphaMax
);

volScalarField alpha2
(
     IOobject
     (
         "alpha2",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
     //alphaMax + (alphaMin - alphaMax)*(vf2*(1 + q)/(vf2 + q))
     vf1Pen*(vf2Pen*alphaMin + (1 - vf2Pen)*alphaMax) + (1 - vf1Pen)*alphaMax
);


Info<< "Reading field p1\n" << endl;
volScalarField p1
(
    IOobject
    (
        "p1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p2\n" << endl;
volScalarField p2
(
    IOobject
    (
        "p2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U1\n" << endl;
volVectorField U1
(
    IOobject
    (
        "U1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U2\n" << endl;
volVectorField U2
(
    IOobject
    (
        "U2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

//#include "createPhi.H"
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading/calculating face flux field phi1\n" << endl;
surfaceScalarField phi1
(
    IOobject
    (
        "phi1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U1)
);

Info<< "Reading/calculating face flux field phi2\n" << endl;
surfaceScalarField phi2
(
    IOobject
    (
        "phi2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U2)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p1, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p1.name());


singlePhaseTransportModel laminarTransport(U1, phi1);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U1, phi1, laminarTransport)
);

setRefCell(p2, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p2.name());


//singlePhaseTransportModel laminarTransport(U2, phi2);

//autoPtr<incompressible::turbulenceModel> turbulence
//(
//    incompressible::turbulenceModel::New(U2, phi2, laminarTransport)
//);

#include "createMRF.H"
#include "createFvOptions.H"



