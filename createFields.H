Info<< "Reading the field c\n" << endl;

// Concentration of ions in the electrode
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading the field p\n" << endl;

// Electrical potential of the macropore
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading the field pD\n" << endl;

// Donan layer potential
volScalarField pD
(
    IOobject
    (
        "pD",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading velocity vector U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField coeff_conv
(
    IOobject
    (
        "coeff_conv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField coeff_diff_1
(
    IOobject
    (
        "coeff_diff_1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField coeff_diff_2
(
    IOobject
    (
        "coeff_diff_2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField src_coeff_p
(
    IOobject
    (
        "src_coeff_p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField src_coeff_c
(
    IOobject
    (
        "src_coeff_c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading transportProperties: De, Pmi, and Pma\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Diffusivity in electrode
dimensionedScalar De
(
    transportProperties.lookup("De")
);

// dimensionedScalar D
// (
//     transportProperties.lookup("D")
// );

// Micro- and macro-pore porosity of the electrode
dimensionedScalar Pmi
(
    transportProperties.lookup("Pmi")
);

dimensionedScalar Pma
(
    transportProperties.lookup("Pma")
);

Info<< "Reading Electro-Chemical Properties: CapSt, VT, F\n" << endl;

IOdictionary electrochemicalProperties
(
    IOobject
    (
        "electrochemicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Electrode voltage
dimensionedScalar EV
(
    electrochemicalProperties.lookup("EV")
);

// Non-electrostatic adsorption parameter
dimensionedScalar Mu
(
    electrochemicalProperties.lookup("Mu")
);

// Stern layer (or microscopic) Capacitance
dimensionedScalar CapSt
(
    electrochemicalProperties.lookup("CapSt")
);

// Thermal voltage
dimensionedScalar VT
(
    electrochemicalProperties.lookup("VT")
);

// Faraday Constant
dimensionedScalar F
(
    electrochemicalProperties.lookup("F")
);

