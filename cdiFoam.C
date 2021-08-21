/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    cdiFoam

Description
    Solves the steady or transient transport equation for a passive scalar.
    Here the passive scalars are concentration and electrical potential.

\*---------------------------------------------------------------------------*/

// Including necessary class definition and libraries ( Header files )
#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Creating objects from the above defined classes
    // These .H files are basically code snippets ( "inclusion" files )
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;


    // Trigger the objects to store the old value of theirs at each iteration
    c.oldTime();
    p.oldTime();
    pSt.oldTime();
    pD.oldTime();
    w.oldTime();
    q.oldTime();

    // Define constant coefficients before the loop to only be calculated once
    dimensionedScalar coeff1 = (Pmi/Pma);
    dimensionedScalar coeff2 = 1e-9 * coeff1 * (VT/De);
    dimensionedScalar coeff3 = -2 * F / CapSt;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Potential Equation
        fvScalarMatrix pEqn
        (
          - fvm::laplacian(c, p)
          ==
            coeff2 * fvc::ddt(src_coeff_p, q)
        );

        pEqn.relax();
        pEqn.solve();

        while (simple.correctNonOrthogonal())
        {
            // Concentration Equation
            fvScalarMatrix cEqn
            (
                fvm::ddt(c)
              + coeff_conv * fvm::div(phi, c)
              // - fvm::laplacian(De, c)
              - coeff_diff_1 * fvm::laplacian(De, c) // spacer *2 : *0
              - coeff_diff_2 * fvm::laplacian(De, c) // electrode *0 : *1
              ==
                coeff1 * fvc::ddt(src_coeff_c, w)
            );

            cEqn.relax();
            cEqn.solve();

        }

        // q = -c * exp(Mu) * sinh(pD/VT);
        // // Stern layer potential (Volt)
        // pSt = -2 * F * q / CapSt;
        // // Donan layer potential (Volt)
        // pD = -1 * VT *asinh(q/ (5 * c));
        // // Volumeteric ions concentration in micropores (mol/m^3)
        // w = sqrt(q * q + 4 * c * c);
        // // Electrical potential
        // // p = EV - pD - pSt;

        q = -c * exp(Mu) * sinh(pD/VT);
        // Stern layer potential (Volt)
        pSt = coeff3 * q;
        // Donan layer potential (Volt)
        pD = -1 * max(p) + EV - max(pSt);
        // Charge density in the micropores (mol/m^3)
        // Volumetric ions concentration in micropores (mol/m^3)
        //w = c * exp(Mu) * cosh(pD/VT);
        w = sqrt(q * q + 4 * c * c);


        runTime.write();
    }

    Info<< endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
