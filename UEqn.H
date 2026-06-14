fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
//- fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  + turbulence->divDevRhoReff(U)
  ==
    fvOptions(rho, U)
);

UEqn.relax();

fluid.addInterfacePorosity(UEqn);


// --- LPBF recoil and Marangoni momentum source ----------------------------

const volScalarField& alphaFeLPBF =
    mesh.lookupObject<volScalarField>("alpha.Fe");

const volScalarField& alphaSiLPBF =
    mesh.lookupObject<volScalarField>("alpha.Si");

const volScalarField& alphaSteelLPBF =
    mesh.lookupObject<volScalarField>("alpha.steel");

volScalarField alphaMetalMomentum
(
    IOobject
    (
        "alphaMetalMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaFeLPBF + alphaSiLPBF + alphaSteelLPBF
);

forAll(alphaMetalMomentum, celli)
{
    alphaMetalMomentum[celli] =
        min(max(alphaMetalMomentum[celli], scalar(0)), scalar(1));
}
alphaMetalMomentum.correctBoundaryConditions();

volVectorField gradAlphaMomentum
(
    IOobject
    (
        "gradAlphaMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(alphaMetalMomentum)
);

volScalarField magGradAlphaMomentum
(
    IOobject
    (
        "magGradAlphaMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(gradAlphaMomentum)
);

volVectorField nHatMomentum
(
    IOobject
    (
        "nHatMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    gradAlphaMomentum
   /(magGradAlphaMomentum + dimensionedScalar("smallGrad", dimless/dimLength, SMALL))
);

volVectorField gradTMomentum
(
    IOobject
    (
        "gradTMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(T)
);

volScalarField liquidFractionMomentum
(
    IOobject
    (
        "liquidFractionMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("liquidFractionMomentum", dimless, 0.0)
);

volScalarField recoilPressure
(
    IOobject
    (
        "recoilPressure",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "recoilPressure",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        0.0
    )
);

volVectorField recoilForce
(
    IOobject
    (
        "recoilForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "recoilForce",
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector::zero
    )
);

volVectorField marangoniForce
(
    IOobject
    (
        "marangoniForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "marangoniForce",
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector::zero
    )
);

const scalar RgasM = 8.314462618;

const scalar TsolFeM = 1700.0;
const scalar TliqFeM = 1811.0;
const scalar TboilFeM = 3134.0;
const scalar LvFeM = 349800.0;
const scalar MFeM = 0.05545;
const scalar dSigmaFe = -0.00049;

const scalar TsolSiM = 1600.0;
const scalar TliqSiM = 1687.0;
const scalar TboilSiM = 3538.0;
const scalar LvSiM = 1280000.0;
const scalar MSiM = 0.028;
const scalar dSigmaSi = -1.3e-05;

const scalar TsolSteelM = 1650.0;
const scalar TliqSteelM = 1800.0;
const scalar TboilSteelM = 3134.0;
const scalar LvSteelM = 349800.0;
const scalar MSteelM = 0.05545;
const scalar dSigmaSteel = -0.00049;

forAll(T, celli)
{
    const scalar metalRaw = alphaMetalMomentum[celli];

    if (metalRaw < 1e-8)
    {
        liquidFractionMomentum[celli] = 0.0;
        recoilPressure[celli] = 0.0;
        recoilForce[celli] = vector::zero;
        marangoniForce[celli] = vector::zero;
        continue;
    }

    const scalar metal = max(metalRaw, SMALL);

    const scalar fFe = alphaFeLPBF[celli]/metal;
    const scalar fSi = alphaSiLPBF[celli]/metal;
    const scalar fSteel = alphaSteelLPBF[celli]/metal;

    const scalar Tsol =
        fFe*TsolFeM + fSi*TsolSiM + fSteel*TsolSteelM;

    const scalar Tliq =
        fFe*TliqFeM + fSi*TliqSiM + fSteel*TliqSteelM;

    const scalar Tboil =
        fFe*TboilFeM + fSi*TboilSiM + fSteel*TboilSteelM;

    const scalar Lv =
        fFe*LvFeM + fSi*LvSiM + fSteel*LvSteelM;

    const scalar M =
        fFe*MFeM + fSi*MSiM + fSteel*MSteelM;

    const scalar dSigma =
        fFe*dSigmaFe + fSi*dSigmaSi + fSteel*dSigmaSteel;

    const scalar Ti = T[celli];

    if (Ti <= Tsol)
    {
        liquidFractionMomentum[celli] = 0.0;
    }
    else if (Ti >= Tliq)
    {
        liquidFractionMomentum[celli] = 1.0;
    }
    else
    {
        liquidFractionMomentum[celli] =
            (Ti - Tsol)/max(Tliq - Tsol, SMALL);
    }

    if (Ti > Tboil && liquidFractionMomentum[celli] > SMALL)
    {
        recoilPressure[celli] =
            0.54*101325.0
           *Foam::exp(scalar(Lv*M*(Ti - Tboil)/(RgasM*Ti*Tboil)));
    }

    recoilForce[celli] =
        recoilPressure[celli]
       *magGradAlphaMomentum[celli]
       *nHatMomentum[celli]
       *liquidFractionMomentum[celli];

    const vector gradTn =
        (gradTMomentum[celli] & nHatMomentum[celli])*nHatMomentum[celli];

    const vector gradTt = gradTMomentum[celli] - gradTn;

    marangoniForce[celli] =
        dSigma
       *gradTt
       *magGradAlphaMomentum[celli]
       *liquidFractionMomentum[celli];
}

volVectorField lpbfMomentumSource
(
    IOobject
    (
        "lpbfMomentumSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    recoilForce + marangoniForce
);

volScalarField marangoniForceMag
(
    IOobject
    (
        "marangoniForceMag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(marangoniForce)
);

Info<< "Diagnostic LPBF momentum: max recoilPressure = "
    << max(recoilPressure).value()
    << " , max recoilForce = " << max(mag(recoilForce)).value()
    << " , max marangoniForceMag = " << max(marangoniForceMag).value()
    << " , max lpbfMomentumSource = " << max(mag(lpbfMomentumSource)).value()
    << endl;

// --------------------------------------------------------------------------


if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
                fluid.surfaceTensionForce()
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
      + lpbfMomentumSource
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
