Info<< "Reading thermal properties\n" << endl;

IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

bool thermalStress(thermalProperties.get<bool>("thermalStress"));

/*volScalarField threeKalpha
(
    IOobject
    (
        "threeKalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimensionSet(0, 2, -2 , -1, 0), Zero)
);*/

dimensionedScalar alphatherm_max("alphatherm_max",  inv(dimTemperature),readScalar(thermalProperties.lookup("alphatherm_max")));
dimensionedScalar alphatherm_min("alphatherm_min",  inv(dimTemperature),readScalar(thermalProperties.lookup("alphatherm_min")));
volScalarField alphatherm(qd*(1-vf1)/(qd+vf1)*(alphatherm_max - alphatherm_min) + alphatherm_min);
volScalarField threeKalpha(threeK*alphatherm);
