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

Info<< "Reading the field q\n" << endl;

// Micropore charge
volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading the field w\n" << endl;

// Micropore Consentration
volScalarField w
(
    IOobject
    (
        "w",
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

Info<< "Reading the field pSt\n" << endl;

// Stern layer potential
volScalarField pSt
(
    IOobject
    (
        "pSt",
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

// vim: ft=foam
