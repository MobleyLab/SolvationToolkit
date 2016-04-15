#!/usr/bin/env python
######################################################################
# SolvationToolkit: A toolkit for setting up molecular simulations of mixtures
# Copyright 2011-2015 UC Irvine and the Authors
#
# Authors: David Mobley and Gaetano Calabro
# With thanks to Kyle Beauchamp, whose liquid_tools.py provided an initial basis for this code
# (https://github.com/choderalab/LiquidBenchmark/blob/master/src/simulation/liquid_tools.py) in April 2015
#
#This library is free software; you can redistribute it and/or
#modify it under the terms of the GNU Lesser General Public
#License as published by the Free Software Foundation; either
#version 2.1 of the License, or (at your option) any later version.
#
#This library is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with this library; if not, see <http://www.gnu.org/licenses/>.
######################################################################


from solvated_mixtures import *
from openeye.oechem import *
from openeye.oeiupac import *

#In this particular instance I'll just look at six solutes/solvent mixtures (not an all-by-all combination) which are pre-specified
#solute names
solutes = ['phenol', 'toluene', 'benzene']
#Solvent names
solvents = ['cyclohexane', 'cyclohexane', 'cyclohexane']

#Generate SMILES strings for these, as they're used by the mixture class
def get_smiles(name):
    mol = OEMol()
    if not OEParseIUPACName( mol, name):
        raise ValueError("Error: The supplied name '%s' could not be parsed." % name )
    return OECreateIsoSmiString( mol )

solute_smiles = [ get_smiles(name) for name in solutes ]
solvent_smiles = [ get_smiles(name) for name in solvents ]

#Number of solute/solvent molecules
#NOTE: These are just placeholders and you should update based on your intended system size and solute/solvent size. Probably a better general strategy is to estimate the number of molecules you want based on a target box size, given an estimate of the density, etc. 
Nsolu = 1
Nsolv = 800

#Construct systems
for idx in range( len( solutes) ):
    builder = MixtureSystem( [ solutes[idx], solvents[idx]], [solute_smiles[idx], solvent_smiles[idx] ], [ Nsolu, Nsolv], 'data/', solute_index=0)  
    builder.run( just_build = True )
    builder.convert_to_gromacs()
