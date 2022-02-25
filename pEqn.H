#pragma once

fvScalarMatrix pEqn
(
     - fvm::laplacian(c, p)
     ==
      (-1) * src_coeff_p    // This must be a positive term
           * porosityRatio
           * exp(Mu)
           * (VT/De)
           * (
                 sinh(pD/VT) * fvc::ddt(c)
               + c * cosh(pD/VT) * fvc::ddt(pD) / VT
             )
);

pEqn.relax();
pEqn.solve();
