#!/bin/env python

from solvated_mixtures import *

#In this particular instance I'll just look at six solutes/solvent mixtures (not an all-by-all combination) which are pre-specified
#solute names
solutes = ['phenol', 'toluene', 'benzene', 'methane', 'ethanol', 'naphthalene']
#Solvent names
solvents = ['cyclohexane', 'cyclohexane', 'cyclohexane', 'octanol', 'octanol', 'octanol']

#The mixture class taxes CAS strings so I need to pull these

solute_CAS = [ openmoltools.cirpy.resolve( name, 'CAS') for name in solutes ]
solvent_CAS = [ openmoltools.cirpy.resolve( name, 'CAS') for name in solvents ]


#Number of solute/solvent molecules
Nsolu = 3
Nsolv = 100


#Construct systems
for idx in range( len( solutes) ):
    builder = MixtureSystem( [ solutes[idx], solvents[idx]], [ Nsolu, Nsolv], 'data/')  
    builder.run( just_build = True )
    builder.convert_to_gromacs()

#Maybe down here I should try one which is a single solute in a binary mixture