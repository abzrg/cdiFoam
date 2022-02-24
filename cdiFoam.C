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

// handler for printing stuff
#define print(X) Info << X << endl

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
    pD.oldTime();

    // // Define constant coefficients before the loop to only be calculated once
    dimensionedScalar porosityRatio = (Pmi/Pma);

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Electrical Potential (p)
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

        dimensionedScalar dummy_coeff(1e-9);
        volScalarField pD_coeff(2 * src_coeff_c * (F / CapSt) * (c / VT) * cosh(pD/VT));

        // Donnan Potential (pD)
        fvScalarMatrix pDEqn
        (
             fvm::ddt(1+pD_coeff,pD)
             // + src_coeff_c * dummy_coeff * fvm::laplacian(De,pD) // volt * m^2 / s / m^2
             ==
             // src_coeff_c  *  dummy_coeff * fvc::laplacian(De,pD)
             - fvc::ddt(p)
             -2 * src_coeff_c     // This must be a negative term
                  * F / CapSt * exp(Mu) * sinh(pD/VT) * fvc::ddt(c)
        );

        pDEqn.relax();
        pDEqn.solve();

        while (simple.correctNonOrthogonal())
        {
            volScalarField c_coeff(src_coeff_c * exp(Mu) * porosityRatio * cosh(pD/VT));

            // Ion concentration (c)
            fvScalarMatrix cEqn
            (
                 fvm::ddt(1+c_coeff,c)
                 + coeff_conv * fvm::div(phi, c)
                 - coeff_diff_1 * fvm::laplacian(De, c)         // spacer *2 : *0
                 - coeff_diff_2 * fvm::laplacian(De, c)         // electrode *0 : *1
                 ==
                 (-1) * src_coeff_c * exp(Mu) * porosityRatio   // This must be a negative term
                      * (
                          c * sinh(pD/VT) * fvc::ddt(pD) / VT
                        )
            );

            cEqn.relax();
            cEqn.solve();
        }

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
