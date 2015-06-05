#!/bin/env python

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
Nsolu = 3
Nsolv = 100

#Construct systems
for idx in range( len( solutes) ):
    builder = MixtureSystem( [ solutes[idx], solvents[idx]], [solute_smiles[idx], solvent_smiles[idx] ], [ Nsolu, Nsolv], 'data/', solute_index=0)  
    builder.run( just_build = True )
    builder.convert_to_gromacs()
