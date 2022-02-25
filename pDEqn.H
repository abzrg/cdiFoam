#pragma once

// Dummy coefficient for diminishing the effect of two dummy laplacian terms
dimensionedScalar dummyCoeff(1e-9);

// Coefficient field for first time derivative (implicit) discretization
volScalarField pDTemporalCoeff
(
    2 * src_coeff_c * (F / CapSt) * (c / VT) * cosh(pD/VT)
);

// Coefficient for (implicit/explicit and explicit) source terms
volScalarField pDSrcCoeff
(
    -2 * src_coeff_c * (F / CapSt) * exp(Mu) * fvc::ddt(c)
);

fvScalarMatrix pDEqn
(
     fvm::ddt((1 + pDTemporalCoeff), pD)
     + src_coeff_c * dummyCoeff * fvm::laplacian(De, pD) // dummy (implicit)
     ==
     src_coeff_c  *  dummyCoeff * fvc::laplacian(De, pD) // dummy (explicit)
     - fvc::ddt(p)
     + (1/VT) * fvm::SuSp(pDSrcCoeff, pD)
     + pDSrcCoeff * (pow(pD, 3) / pow(VT, 3) / 6)
);

pDEqn.relax();
pDEqn.solve();
