#Initially written as liquid_tools.py by Kyle Beauchamp.
#Downloaded 4/30/15 from https://github.com/choderalab/LiquidBenchmark/blob/master/src/simulation/liquid_tools.py
#Revised 4/30/15 by David Mobley to test automation of setting up our own solution-phase simulations.

import numpy as np
import os
import itertools
import mdtraj as md

import openmoltools

def make_path(filename):
    try:
        path = os.path.split(filename)[0]
        os.makedirs(path)
    except OSError:
        pass


class AmberMixtureSystem(object):
    """A pipeline for simulating liquid mixtures using amber parameter files.
    
    Parameters
    ----------
    cas_strings : list(str)
        CAS strings for each component of the mixture
    n_monomers: list(int)
        Number of each type of molecule
    data_path
        Directory tree in which to store the data
    """

    def __init__(self, cas_strings, n_monomers, DATA_PATH ):

        self.cas_strings = cas_strings
        self.n_monomers = n_monomers

        identifier = list(itertools.chain(cas_strings, [str(n) for n in n_monomers]) )
        self.identifier = '_'.join(identifier)        
        
        self.monomer_pdb_filenames = [DATA_PATH + "monomers/" + string + ".pdb" for string in self.cas_strings]
        self.box_pdb_filename = DATA_PATH + "packmol_boxes/" + self.identifier + ".pdb"
        
        self.inpcrd_filename = DATA_PATH + "tleap/" + self.identifier + ".inpcrd"
        self.prmtop_filename = DATA_PATH + "tleap/" + self.identifier + ".prmtop"
        self.top_filename = DATA_PATH + "gromacs/" + self.identifier + ".top"
        self.gro_filename = DATA_PATH + "gromacs/" + self.identifier + ".gro"
        
        
        self.gaff_mol2_filenames = [DATA_PATH + "monomers/" + string + ".mol2" for string in self.cas_strings]
        self.frcmod_filenames = [DATA_PATH + "monomers/" + string + ".frcmod" for string in self.cas_strings]
        
        make_path(DATA_PATH + 'monomers/')
        make_path(DATA_PATH + 'packmol_boxes/')
        make_path(DATA_PATH + 'tleap/')
        make_path(DATA_PATH + 'gromacs/') 

    @property
    def smiles_strings(self):
        self._smiles_strings = []
        for mlc in self.cas_strings:
            self._smiles_strings.append(openmoltools.cirpy.resolve(mlc, 'smiles'))
        
        return self._smiles_strings

    def run(self, just_build=False):
        """Build mol2 monomers, packmol boxes, inpcrd files, equilibrate, and run production."""
        self.build_monomers()
        self.build_boxes()
        if not just_build:
            self.equilibrate()
            self.production()

    def build_monomers(self):
        """Generate GAFF mol2 and frcmod files for each chemical."""
        for k, smiles_string in enumerate(self.smiles_strings):
            mol2_filename = self.gaff_mol2_filenames[k]
            frcmod_filename = self.frcmod_filenames[k]
            if not (os.path.exists(mol2_filename) and os.path.exists(frcmod_filename)):
                #Convert SMILES strings to mol2 and frcmod files for antechamber
                openmoltools.openeye.smiles_to_antechamber(smiles_string, mol2_filename, frcmod_filename)

        #Generate unique residue names for molecules in mol2 files
        openmoltools.utils.randomize_mol2_residue_names( self.gaff_mol2_filenames )


    def build_boxes(self):
        """Build an initial box with packmol and use it to generate AMBER files."""
        if not os.path.exists(self.box_pdb_filename):
            packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], self.n_monomers)
            packed_trj.save(self.box_pdb_filename)

        if not (os.path.exists(self.inpcrd_filename) and os.path.exists(self.prmtop_filename)):
            tleap_cmd = openmoltools.amber.build_mixture_prmtop(self.gaff_mol2_filenames, self.frcmod_filenames, self.box_pdb_filename, self.prmtop_filename, self.inpcrd_filename)

    def convert_to_gromacs(self):
        """Take AMBER format prmtop and crd files and convert to GROMACS format."""

        openmoltools.utils.convert_via_acpype( self.identifier, self.prmtop_filename, self.inpcrd_filename, self.top_filename, self.gro_filename ) 
