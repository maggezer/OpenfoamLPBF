# OpenfoamLPBF

OpenfoamLPBF is an OpenFOAM-based solver setup for laser powder bed fusion (LPBF) and in-situ alloying simulations.

The current development target is a multiphase LPBF solver that can model laser heating, melting, evaporation-related losses, recoil pressure and Marangoni-driven flow during localized laser-material interaction.

> Development status: experimental / work in progress.

## Solver

The main solver executable is:

```text
lpbfIcoReactingMultiphaseInterFoam
```

The solver is based on a multiphase incompressible, non-isothermal formulation and is being extended for LPBF-specific physics.

Main source file:

```text
icoReactingMultiphaseInterFoam.C
```

Important equation files:

```text
TEqn.H
UEqn.H
pEqn.H
YEqns.H
createFields.H
createFieldRefs.H
```

## Current physics included

The current version includes development-stage implementations of:

* Multiphase VOF-style phase tracking
* Fe, Si, steel and argon phases
* Volumetric Gaussian laser heat source
* Temperature-dependent melting / liquid fraction treatment
* Radiation heat loss
* Evaporation-related heat loss
* Recoil pressure momentum source
* Marangoni momentum source
* Pressure-correction coupling for LPBF momentum contribution
* Example LPBF case under `case/`

## Repository structure

```text
.
├── Make/                         Solver build configuration
├── laserDTRM/                    Laser / radiation model library
├── case/                         Example LPBF case
│   ├── 0/                        Initial and boundary fields
│   ├── constant/                 Material, phase and physical properties
│   └── system/                   Numerical and mesh dictionaries
├── Allwmake                      Build script
├── Allwclean                     Build clean script
├── icoReactingMultiphaseInterFoam.C
├── TEqn.H
├── UEqn.H
├── pEqn.H
└── README.md
```

## Requirements

This repository is being developed with:

```text
OpenFOAM v2412
```

Before building, load your OpenFOAM environment. For example:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

Adjust the path according to your local OpenFOAM installation.

## Build

From the repository root:

```bash
./Allwmake
```

To clean the compiled solver and library:

```bash
./Allwclean
```

## Run the example case

From the repository root:

```bash
cd case
./Allrun
```

The current `case/Allrun` script performs:

```bash
blockMesh
setFields
lpbfIcoReactingMultiphaseInterFoam
```

To clean generated case files:

```bash
cd case
./Allclean
```

## Mesh files

Generated mesh files are not tracked in Git.

The following directory is intentionally ignored:

```text
case/constant/polyMesh/
```

The mesh should be regenerated locally using:

```bash
blockMesh
```

This keeps the repository lightweight and avoids committing large generated mesh files.

## Case phases

The current example case uses the following phases:

```text
argon
Fe
Si
steel
```

The corresponding thermophysical property files are stored under:

```text
case/constant/
```

Expected files:

```text
thermophysicalProperties.argon
thermophysicalProperties.Fe
thermophysicalProperties.Si
thermophysicalProperties.steel
```

## Development notes

This repository is currently being cleaned and refactored.

Current priorities:

1. Clean repository structure
2. Keep generated OpenFOAM files out of Git
3. Make the example case reproducible from dictionaries
4. Move hard-coded LPBF parameters into OpenFOAM dictionaries
5. Improve initialization of powder, substrate and gas regions
6. Validate the solver with simple benchmark cases
7. Improve documentation for physical models and numerical settings

## Known limitations

The current solver is experimental.

Some LPBF parameters are still hard-coded in the source files. These should be moved into dictionary files such as:

```text
case/constant/LPBFProperties
```

The example case should be checked carefully before physical interpretation, especially:

* phase initialization,
* laser parameters,
* material properties,
* boundary conditions,
* mesh resolution,
* time-step control,
* energy-source normalization.

## License

No license has been specified yet.

Before sharing or publishing the repository more widely, an explicit license should be added.
