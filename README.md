# SolvationToolkit
In-house tools for setting up arbitrary (currently non-water) solute-solvent mixtures for simulation in GROMACS

## Current status

At this point what's here consists primarily of `solvated_mixtures.py` which sets up a class structure to use openmoltools (and related utilities) to automatically set up simulations of mixtures of solvents/solutes, and `test_drive.py` which gives an example driver script using this to set up simulations of several mixtures. 

At this point its usage is limited to non-water solutes/solvents, but other than that it should basically work for any solute/solvent covered by GAFF.  


## Prerequisites
* Auto-installable:
  * openmoltools v0.6.5 or later
  * mdtraj (via omnia channel on conda)
  * OpenEye tools
  * packmol (via omnia channel on conda)
* Must be installed manually
  * ParmEd (from source)

## Issues
* Currently difficult to single out one molecule as a solute (will require manual topology editing)
* Topologies are inefficient (defining each molecule in the full system regardless of whether it has already been defined)
