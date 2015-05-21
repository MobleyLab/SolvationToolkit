#Initially written as liquid_tools.py by Kyle Beauchamp.
#Downloaded 4/30/15 from https://github.com/choderalab/LiquidBenchmark/blob/master/src/simulation/liquid_tools.py
#Revised 4/30/15 by David Mobley to test automation of setting up our own solution-phase simulations.

import numpy as np
import os
import itertools
import mdtraj as md
import copy

import openmoltools

def make_path(filename):
    try:
        path = os.path.split(filename)[0]
        os.makedirs(path)
    except OSError:
        pass


class MixtureSystem(object):
    """A pipeline for simulating liquid mixtures using amber and gromacs parameter files.
    
    Parameters
    ----------
    labels: list(str)
        Labels or names for components of the system; full name will be constructed from the labels of the components plus the number of each component. 
    smiles_strings : list(str)
        SMILES strings (isomeric if chiral) for the components of the system
    n_monomers: list(int)
        Number of each type of molecule
    data_path : str
        Directory tree in which to store the data
    solute_index : int/str, optional. Default: "auto"
        Optional parameter to specify which of the components (in the list of specified components)
        will be treated as a solute in constructing GROMACS topology files 
        (which means that a single molecule of this component will be singled out as the 'solute'
        in the resulting GROMACS topology file). Valid options are 'auto' (pick the first component present with n_monomers = 1,
        otherwise the first component), None (don't pick any), or an integer (pick the component smiles_strings[solute_index].

    Notes
    -----
    Monomer components of the system will be stored with filenames 'labelN' where N is the component number, i.e. 0, 1, ... len(smiles_strings).

    Limitations
    -----------
    Existing files with the same name present in the data directory tree may be overwritten. This results in a limitation/failure in a small (and probably random) fraction of cases if multiple systems involving the same monomers are written into the same data directory. Specifically, openmoltools.amber.build_mixture_prmtop requires that each mol2 file for a component have a unique residue name, which is handled automatically by openmoltools when constructing monomers (each is assigned a unique random residue name). However, if these are overwritten with other monomers (i.e. if we set up, say, 'octanol' in the same directory twice) which by chance end up with non-unique residue names then amber.build_mixture_prmtop will fail with a ValueError. This can be avoided by ensuring that if you are constructing multiple MixtureSystems involving the same monomers, your data directories are different. This issue also will likely be fixed when openmoltools switches to topology merging via ParmEd rather than tleap, as unique residue names are built into ParmEd in a better way. 
    """

    def __init__(self, labels, smiles_strings, n_monomers, DATA_PATH, solute_index = 'auto' ):

        self.smiles_strings = smiles_strings
        self.n_monomers = n_monomers
        self.solute_index = solute_index
        self.n_components = len( smiles_strings )
        self.labels = labels

        assert len(smiles_strings) == len(n_monomers), "Number of provided smiles strings must equal the number of each type of molecule provided." 
        assert len(smiles_strings) == len(labels), "Number of provided smiles strings must equal the number of labels provided." 

        identifier = list(itertools.chain(self.labels, [str(n) for n in n_monomers]) )
        self.identifier = '_'.join(identifier)  
        
        self.monomer_pdb_filenames = [DATA_PATH + "monomers/" + string + ".pdb" for string in self.labels]
        self.box_pdb_filename = DATA_PATH + "packmol_boxes/" + self.identifier + ".pdb"
        
        #input_crd_filename and prmtop_filename stores the AMBER filenames of the solvated molecules
        self.inpcrd_filename = DATA_PATH + "tleap/" + self.identifier + ".inpcrd"
        self.prmtop_filename = DATA_PATH + "tleap/" + self.identifier + ".prmtop"
        
        #top_filename and gro_filename stores the GROMACS filename of the solvated molecules
        self.top_filename = DATA_PATH + "gromacs/" + self.identifier + ".top"
        self.gro_filename = DATA_PATH + "gromacs/" + self.identifier + ".gro"
        
        
        self.gaff_mol2_filenames = [DATA_PATH + "monomers/" + string + ".mol2" for string in self.labels]
        self.frcmod_filenames = [DATA_PATH + "monomers/" + string + ".frcmod" for string in self.labels]
        
        #input_crd_filenames and prmtop_filenames stores the AMBER filenames of the molecules without solvation
        self.inpcrd_filenames = [DATA_PATH + "tleap/" + string + ".inpcrd" for string in self.labels]
        self.prmtop_filenames = [DATA_PATH + "tleap/" + string + ".prmtop" for string in self.labels]
        
        #top_filenames and gro_filenames stores the GROMACS filenames of the molecules without solvation
        self.gro_filenames = [DATA_PATH + "gromacs/" + string + ".gro" for string in self.labels]
        self.top_filenames = [DATA_PATH + "gromacs/" + string + ".top" for string in self.labels]

        make_path(DATA_PATH + 'monomers/')
        make_path(DATA_PATH + 'packmol_boxes/')
        make_path(DATA_PATH + 'tleap/')
        make_path(DATA_PATH + 'gromacs/') 

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
            inpcrd_filename = self.inpcrd_filenames[k]
            prmtop_filename = self.prmtop_filenames[k]
            gro_filename = self.gro_filenames[k]
            top_filename = self.top_filenames[k]
            if not (os.path.exists(mol2_filename) and os.path.exists(frcmod_filename)):
                #Convert SMILES strings to mol2 and frcmod files for antechamber
                openmoltools.openeye.smiles_to_antechamber(smiles_string, mol2_filename, frcmod_filename)
                #Generate amber coordinate and topology files for the unsolvated molecules
                mol_name = os.path.basename(gro_filename).split('.')[0]
                openmoltools.utils.run_tleap(mol_name, mol2_filename,frcmod_filename, prmtop_filename, inpcrd_filename)
                #Generate gromacs coordinate and topology coordinate files for the unsovated molecules
                openmoltools.utils.convert_via_acpype(mol_name, prmtop_filename, inpcrd_filename, top_filename, gro_filename)
                
        #Generate unique residue names for molecules in mol2 files
        openmoltools.utils.randomize_mol2_residue_names( self.gaff_mol2_filenames )
        
    def build_boxes(self):
        """Build an initial box with packmol and use it to generate AMBER files."""
        if not os.path.exists(self.box_pdb_filename):
            size = openmoltools.packmol.approximate_volume_by_density( self.smiles_strings, self.n_monomers )
            packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], self.n_monomers, box_size = size)
            packed_trj.save(self.box_pdb_filename)

        if not (os.path.exists(self.inpcrd_filename) and os.path.exists(self.prmtop_filename)):
            tleap_cmd = openmoltools.amber.build_mixture_prmtop(self.gaff_mol2_filenames, self.frcmod_filenames, self.box_pdb_filename, self.prmtop_filename, self.inpcrd_filename)

    def convert_to_gromacs(self):
        """From AMBER-format prmtop and crd files, generate final solvated GROMACS topology and coordinate files.

        Notes
        -----
        
        Algorithmic notes:
        ------------------
        Separate pathways are used for the coordinate and topology files. 
        Particularly, SOLVATED AMBER prmtop and crd files (created by way of packmol) are 
        converted to generate the solvated GROMACS .gro file. The topology 
        file is created by taking the un-solvated components of the system in GROMACS format 
        and generating a combined GROMACS topology file consisting of the appropriate number of monomers."""

        #Generate solvated topology and coordinate file for the full system via acpype
        #The topology file created here will be overwritten below since we don't need it     
        openmoltools.utils.convert_via_acpype( self.identifier, self.prmtop_filename, self.inpcrd_filename, self.top_filename, self.gro_filename ) 

        #Figure out what we're treating as the solute (if anything)
        monomer_present = False
        if self.solute_index=='auto':
            #Check which of the molecules is present in qty 1
            try:
                self.solute_index = self.n_monomers.index(1)
                monomer_present = True
            except ValueError:
                #If none is present in qty 1, then use the first 
                self.solute_index = 0
            
        #If we aren't treating anything as the solute, construct a combined topology file
        if self.solute_index == None:
            openmoltools.gromacs.merge_topologies( self.top_filenames, self.top_filename, 'mixture', molecule_numbers = self.n_monomers )
        
        else:
            #Handle case where a particular molecule is specified as the solute 
            
            #Check that the passed solute index is correct
            check_solute_indices = range(0,len(self.n_monomers))
            assert self.solute_index in check_solute_indices and isinstance(self.solute_index, int), "Solute index must be an element of the list: %s. The value passed is: %s" % (check_solute_indices,self.solute_index) 
            
            #If monomer_present is True (if one was already a monomer) then we preserve the same number of components; 
            #otherwise we are increasing the number of components in the topology by one by splitting off a monomer
            if not monomer_present:
                #Increase the number of components and construct new input topologies list (we are making one topology be included twice under two different names)
                self.top_filenames = self.top_filenames[0:self.solute_index] + [self.top_filenames[self.solute_index]] + self.top_filenames[self.solute_index:]
                #Change number of components accordingly
                self.n_monomers[self.solute_index] = self.n_monomers[self.solute_index]-1
                self.n_monomers = self.n_monomers[0:self.solute_index] + [1] + self.n_monomers[self.solute_index:] 
                #Construct names - solute will be specified as such
                names = self.labels[0:self.solute_index] + ['solute'] + self.labels[self.solute_index:]
            #Otherwise we're just changing the name of one of the components and leaving everything else as is   
            else:
                #Only change names
                names = copy.copy( self.labels )
                names[ self.solute_index ] = 'solute'
                 
            #Now merge
            openmoltools.gromacs.merge_topologies(self.top_filenames, self.top_filename, 'mixture', molecule_names = names, molecule_numbers = self.n_monomers ) 
                     
