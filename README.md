##Build Status

[![Build Status](https://travis-ci.org/nividic/SolvationToolkit.svg?branch=master)](https://travis-ci.org/nividic/SolvationToolkit)

# SolvationToolkit
Tools for setting up arbitrary (currently non-water) solute-solvent mixtures for simulation in GROMACS or AMBER formats, or for use in OpenMM. 

This provides the MixtureSystem class (in `solvationtoolkit.solvated_mixtures`) for definition and construction of mixtures as discussed below. Usage examples are below, as well as in `scripts/start.py`.

# Installing

## Prerequisites
* Auto-installable:
  * openmoltools v0.6.5 or later
  * mdtraj (via omnia channel on conda)
  * OpenEye tools
  * packmol (via omnia channel on conda)
  * ParmEd v2.5.1.10 (currently, parmed-dev) or later, such as via omnia channel on conda

## Installing
* Download the archive from GitHub, extract, and type `python setup.py install`
* Hopefully installation via conda is coming very soon

# Usage

## Example usage
```python
# Set up mixture of phenol, toluene, and cyclohexane
# using specified numbers of each component
mixture = MixtureSystem('mydata')
mixture.addComponent(label='phenol', smiles='c1ccc(cc1)O', number=1)
mixture.addComponent(label='toluene', smiles='Cc1ccccc1', number=10)
mixture.addComponent(label='cyclohexane', smiles='C1CCCCC1', number=500)
#Generate output files for AMBER
mixture.build(amber = True)

# Set up simple liquid of TIP3P water
liquid = MixtureSystem()
liquid.addComponent('water')
liquid.build()

# Set up binary mixture of water and methanol
binary_mixture = MixtureSystem()
binary_mixture.addComponent(name='water', mole_fraction=0.2)
# add a filling compound - assumed to be rest of mixture if no
# mole_fraction specified
binary_mixture.addComponent(name='methanol')
#Build for GROMACS
binary_mixture.build(gromacs = True)

# Set up ternary mixture of ethanol, methanol, and water
ternary_mixture = MixtureSystem()
ternary_mixture.addComponent(name='ethanol', mole_fraction=0.2)
ternary_mixture.addComponent(name='methanol', mole_fraction=0.2)
ternary_mixture.addComponent(name='water')

# Set up a system of phenol at infinite dilution (a single molecule)
# in otherwise pure water
infinite_dilution = MixtureSystem()
infinite_dilution.addComponent(name='phenol', mole_fraction=0.0)
infinite_dilution.addComponent(name='water')
```

## Other usage information

One or multiple components must be specified, and each component must have a
name or label (or both), as well as optionally SMILES strings and number of
molecules or mole fraction. Each component must end up with a string which can
be used to construct filenames (which will be taken from "label" or, if not
provided, "name"), and each component must also be unambiguously identified
via a name or SMILES string. SMILES strings are the preferred identifier, but
if SMILES are not provided, it attempts to process the name to a SMILES, and if
there is no name, the label is interpreted as a name.

Except in certain special cases, each component must also have a quantity,
which can either be a number of molecules, or a mole fraction. Numbers of
molecules and mole fractions cannot be mixed; system specification must either
occur via numbers or via mole fractions. Special cases are the case of a single
component, in which no quantity is necessary as the solution is pure and can be
generated based on the target system size, and the case of a filling component,
which is provided without quantity specification along with several other
components with specified mole fractions.

Systems are finalized by executing the build functionality, which has optional
arguments `amber = True` and `gromacs = True` for generating output for the
specified packages. For aid with solvation free energy calculations, build of
GROMACS systems takes an additional optional argument, `solute_index`,
specifying the component number (i.e. in the order ardded, counting from 0) to
treat as the solute.

## Storage of results/intermediate files

Data files are stored in subdirectories of the directory specified by the
constructor, i.e. MixtureSystem('path'); the default location is 'data'.
Within that, subdirectories contain .sdf and .mol2 files for monomers,
the generated .pdb format boxes from packmol, and AMBER and GROMACS files if
generated.

# Other Stuff

## Issues/problems/bug reports
* See issue tracker

## Disclaimers
* This code is currently in approximately beta release status. Use at your own risk. We will almost certainly be making changes to the API in the near future. (A very substantial API change was made in v0.4.0) 
