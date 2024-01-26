{
    volScalarField rAU(1.0/U2Eqn.A());
    volVectorField HbyA(constrainHbyA(rAU*U2Eqn.H(), U2, p2));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U2, p2);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - U2Eqn.H1()); //Do not change H1
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p2)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p2);
    }

    tU2Eqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p2, U2, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix p2Eqn
        (
            fvm::laplacian(rAtU(), p2) == fvc::div(phiHbyA)
        );

        p2Eqn.setReference(pRefCell, pRefValue);

        p2Eqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p2Eqn.flux();
        }
    }

    #include "continuityErrs.H"
    phi2 = phi;

    // Explicitly relax pressure for momentum corrector
    p2.relax();

    // Momentum corrector
    U2 = HbyA - rAtU()*fvc::grad(p2);
    U2.correctBoundaryConditions();
    //fvOptions.correct(U2);
}
