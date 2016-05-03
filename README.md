##Build Status

[![Build Status](https://travis-ci.org/nividic/SolvationToolkit.svg?branch=master)](https://travis-ci.org/nividic/SolvationToolkit)

# SolvationToolkit
Tools for setting up arbitrary (currently non-water) solute-solvent mixtures for simulation in GROMACS or AMBER formats. 

At this point what's here consists primarily of `solvationtoolkit.solvated_mixtures` which sets up a class structure to use openmoltools (and related utilities) to automatically set up simulations of mixtures of solvents/solutes, and `scripts/start.py` which gives an example driver script using this to set up simulations of several mixtures. Automated testing is also implemented. 


## Prerequisites
* Auto-installable:
  * openmoltools v0.6.5 or later
  * mdtraj (via omnia channel on conda)
  * OpenEye tools
  * packmol (via omnia channel on conda)
  * ParmEd

## Issues
* See issue tracker

## Disclaimers
* This code is currently in approximately beta release status. Use at your own risk. We will almost certainly be making changes to the API in the near future.  
