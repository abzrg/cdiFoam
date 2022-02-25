#pragma once

// Coefficient field for first time derivative (implicit) discretization
volScalarField cTemporalCoeff
(
    src_coeff_c * exp(Mu) * porosityRatio * cosh(pD/VT)
);

// Coefficient for (implicit/explicit) source term
volScalarField cSrcCoeff
(
    -1 * src_coeff_c * exp(Mu) * porosityRatio
       * sinh(pD/VT) * fvc::ddt(pD) / VT
);

fvScalarMatrix cEqn
(
     fvm::ddt((1 + cTemporalCoeff), c)
     + coeff_conv * fvm::div(phi, c)
     - coeff_diff_1 * fvm::laplacian(De, c) // spacer *2 : *0
     - coeff_diff_2 * fvm::laplacian(De, c) // electrode *0 : *1
     ==
     fvm::SuSp(cSrcCoeff, c)
);

cEqn.relax();
cEqn.solve();
