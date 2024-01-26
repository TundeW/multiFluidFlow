    // Momentum predictor

    alpha2 = vf1Pen*(vf2Pen*alphaMin + (1 - vf2Pen)*alphaMax) + (1 - vf1Pen)*alphaMax;
    
    //MRF.correctBoundaryVelocity(U2);

    tmp<fvVectorMatrix> tU2Eqn
    (
        fvm::div(phi2, U2)
      - fvm::laplacian(nu2, U2)
      + fvm::Sp(alpha2, U2)
    );
    fvVectorMatrix& U2Eqn = tU2Eqn.ref();

    U2Eqn.relax();

    //fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(U2Eqn == -fvc::grad(p2));

        //fvOptions.correct(U2);
    }
