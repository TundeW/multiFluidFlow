{
    volScalarField rAU(1.0/U1Eqn.A());
    volVectorField HbyA(constrainHbyA(rAU*U1Eqn.H(), U1, p1));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U1, p1);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - U1Eqn.H1()); //Do not change H1
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p1)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p1);
    }

    tU1Eqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p1, U1, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix p1Eqn
        (
            fvm::laplacian(rAtU(), p1) == fvc::div(phiHbyA) //**Check**
        );

        p1Eqn.setReference(pRefCell, pRefValue); //**Check**

        p1Eqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p1Eqn.flux();
        }
    }

    #include "continuityErrs.H"
    phi1 = phi;

    // Explicitly relax pressure for momentum corrector
    p1.relax();

    // Momentum corrector
    U1 = HbyA - rAtU()*fvc::grad(p1); //**Check**
    U1.correctBoundaryConditions();
    //fvOptions.correct(U1);
}
