#!/bin/env python

from solvated_mixtures import *
from openeye.oechem import *
from openeye.oeiupac import *

#In this particular instance I'll just look at six solutes/solvent mixtures (not an all-by-all combination) which are pre-specified
#solute names
solutes = ['phenol', 'toluene', 'benzene']
#Solvent names
solvents = ['cyclohexane', 'cyclohexane', 'octanol']

#Generate SMILES strings for these, as they're used by the mixture class
def get_smiles(name):
    mol = OEMol()
    if not OEParseIUPACName( mol, name):
        raise ValueError("Error: The supplied name '%s' could not be parsed." % name )
    return OECreateIsoSmiString( mol )

solute_smiles = [ get_smiles(name) for name in solutes ]
solvent_smiles = [ get_smiles(name) for name in solvents ]

test1 = ['octanol', 'octanol', 'octanol']
test2 = ['methane','methane','methane']

#Number of solute/solvent molecules
Nsolu = 3
Nsolv = 100
Ntest1 = 4
Ntest2 = 5

#Construct systems
for idx in range( len( solutes) ):
    builder = MixtureSystem( [ solutes[idx], solvents[idx], test1[idx],test2[idx]], [ Nsolu, Nsolv,Ntest1,Ntest2], 'data/', solute_index=3)  
    builder.run( just_build = True )
    builder.convert_to_gromacs()
