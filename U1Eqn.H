    // Momentum predictor
    //#include "multiply.H"
    

    //alpha1 = alphaMax + (alphaMin - alphaMax)*vf1Pen;
    alpha1 = vf1Pen*(vf2Pen*alphaMax + (1 - vf2Pen)*alphaMin) + (1 - vf1Pen)*alphaMax; 

    //MRF.correctBoundaryVelocity(U1);

    tmp<fvVectorMatrix> tU1Eqn
    (
        fvm::div(phi1, U1)
      - fvm::laplacian(nu1, U1)
      + fvm::Sp(alpha1, U1)
      //fvm::div(phi1, U1)
    );
    fvVectorMatrix& U1Eqn = tU1Eqn.ref();

    //write("System", runTime.timeIndex(), U1Eqn, U1);

    //alpha1=1;

    U1Eqn.relax();

    //fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(U1Eqn == -fvc::grad(p1));

        //fvOptions.correct(U1);
    }
