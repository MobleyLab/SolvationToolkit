######################################################################
# SolvationToolkit: A toolkit for setting up molecular simulations of mixtures
# Copyright 2011-2016 UC Irvine and the Authors
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


import numpy as np
import os,sys
import inspect
import itertools
import mdtraj as md
import parmed
import openmoltools
import mol2tosdf


# We require at least ParmEd 2.5.1 because of issues with the .mol2 writer (issue #691 on ParmEd) prior to that.
# Previously 2.0.4 or later was required due to issues with FudgeLJ/FudgeQQ in resulting GROMACS topologies in
#  earlier versions
try: #Try to get version tag
    ver = parmed.version
except: #If too old for version tag, it is too old
    oldParmEd = Exception('ERROR: ParmEd is too old, please upgrade to 2.5.1 or later')
    raise oldParmEd
if ver < (2,5,1):
    raise RuntimeError("ParmEd is too old, please upgrade to 2.5.1 or later")


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
        #check numbers of passed arguments 
        num_args = len(inspect.getargspec(MixtureSystem.__init__).args)
        assert num_args >= 5 , "The number of passed values must be at least 5, given %d" % num_args
        
        #check the types of passed arguments
        check_type_error = False
        if not all(isinstance(i,str) for i in labels):
            check_type_error = True
        if not all(isinstance(i,str) for i in smiles_strings):
            check_type_error = True
        if not all(isinstance(i,int) for i in n_monomers):
            check_type_error = True
        if not isinstance(DATA_PATH, str):
            check_type_error =True
        if not (isinstance(solute_index,str) or isinstance(solute_index,int) or isinstance(solute_index,type(None))):
            check_type_error = True
        if check_type_error :
            raise TypeError('One or more passed types are wrong')
            
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
        
        self.sdf_filenames = [DATA_PATH + "monomers/" + string + ".sdf" for string in self.labels]
        
        #input_crd_filenames and prmtop_filenames stores the AMBER filenames of the molecules without solvation
        self.inpcrd_filenames = [DATA_PATH + "tleap/" + string + ".inpcrd" for string in self.labels]
        self.prmtop_filenames = [DATA_PATH + "tleap/" + string + ".prmtop" for string in self.labels]
        
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
            sdf_filename = self.sdf_filenames[k]
            if not (os.path.exists(mol2_filename) and os.path.exists(frcmod_filename)):
                #Convert SMILES strings to mol2 and frcmod files for antechamber
                openmoltools.openeye.smiles_to_antechamber(smiles_string, mol2_filename, frcmod_filename)
                #Correct the mol2 file partial atom charges to have a total net integer molecule charge  
                mol2f = parmed.formats.Mol2File
                mol2f.write(parmed.load_file(mol2_filename).fix_charges(),mol2_filename)
            #Generate amber coordinate and topology files for the unsolvated molecules
            mol_name = os.path.basename(mol2_filename).split('.')[0]
            openmoltools.amber.run_tleap(mol_name, mol2_filename,frcmod_filename, prmtop_filename, inpcrd_filename)
    
            #Read Mol2 File and write SDF file
            mol2tosdf.writeSDF(mol2_filename, sdf_filename, mol_name)

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
        """From AMBER-format prmtop and crd files, generate final solvated GROMACS topology and coordinate files. Ensure that the desired "solute" (as per solute_index) has a single monomer treated via a unique residue name to allow treatment as a solute separate from other residues of the same name (if desired). The solute will be given residue name "solute" Also, check to see if there are "WAT" residues present, in which case tleap will have re-ordered them to the end of the data file. If so, update data structures accordingly and handle conversion appropriately. 

        Notes
        -----
        Currently, this function ensures that - after AMBER conversion reorders water molecules with residue names 'WAT' to occur last in the resulting parameter/coordinate files - the internal data structures are updated to have the correct order in the relevant lists (labels, smiles_strings, n_monomers). If for some reason GROMACS conversion were removed, these would need to be updated elsewhere. (Probably this should be done anyway, as this is not really a GROMACS issue.)

        """
        #Read in AMBER format parameter/coordinate file and convert in gromacs
        gromacs_topology = parmed.load_file( self.prmtop_filename, self.inpcrd_filename )
     
        #Split the topology into components and check that we have the right number of components
        components = gromacs_topology.split()
        assert len(components)==len(self.n_monomers), "Number of monomers and number of components in the combined topology do not match." 

        ####    HANDLE ORDERING OF WATER   ####
        #Check if any of the residues is named "WAT". If it is, antechamber will potentially have re-ordered it from where it was (it places residues named "WAT" at the end) so it may no longer appear in the order in which we expect.
        resnames = [ components[i][0].residues[0].name for i in range(len(components)) ]
        wat_present = False
        if 'WAT' in resnames:
            wat_present = True

        #Manage presence of WAT residues and possible re-ordering
        if wat_present:
            #If there is a water present, then we MIGHT have re-ordering. Check smiles to find out where it was originally.
            wat_orig_index = self.smiles_strings.index('O')
            #Where is it now?
            wat_new_index = resnames.index('WAT')
            #Reordered? If so, we have to adjust the ordering of n_monomers, smiles_strings, labels,
            # and potentially solute_index. Filenames will be preserved since these were already created
            if wat_orig_index <> wat_new_index:
                #tleap moves water to the end so if they aren't equal, we know where water will be...
                self.n_monomers = self.n_monomers[0:wat_orig_index] + self.n_monomers[wat_orig_index+1:] + [self.n_monomers[wat_orig_index]] 
                self.smiles_strings = self.smiles_strings[0:wat_orig_index] + self.smiles_strings[wat_orig_index+1:] + [self.smiles_strings[wat_orig_index]] 
                self.labels = self.labels[0:wat_orig_index] + self.labels[wat_orig_index+1:] + [self.labels[wat_orig_index] ]
                #Check solute_index and alter if needed
                if not self.solute_index=='auto' and not self.solute_index==None:
                    #Index unchanged if it's before the water
                    if self.solute_index < wat_orig_index:
                        pass
                    #If it is the water, now it is at the end
                    elif self.solute_index == wat_orig_index:
                        self.solute_index = len(self.n_monomers)-1
                    #If it was after the water, then it moved up one position
                    else:
                        self.solute_index -= 1
        ####    END HANDLING OF ORDERING OF WATER   ####
                
             
        #Figure out what we're treating as the solute (if anything)
        if self.solute_index=='auto':
            #Check which of the molecules is present in qty 1
            try:
                self.solute_index = self.n_monomers.index(1)
            except ValueError:
                #If none is present in qty 1, then use the first 
                self.solute_index = 0
            
        #Check that the passed solute index is correct
        check_solute_indices = range(0,len(self.n_monomers))
        assert self.solute_index in check_solute_indices and isinstance(self.solute_index, int) or self.solute_index == None, "Solute index must be an element of the list: %s or None. The value passed is: %s" % (check_solute_indices,self.solute_index)
            
        #Now all we have to do is to change the name of the solute molecule (residue, in ParmEd) and ParmEd will automatically make it a new molecule on write.
        #To do this, first build a list of the residue names we want, by molecule
        resnames = [ ]
        for i in range(len(self.n_monomers)):
            #If this is not the solute, just keep what we had
            if i!=self.solute_index:
                resnames += [ self.labels[i] ] * self.n_monomers[i] 
            #If it is the solute, make the first residue be named solute and the rest what they were already
            else: 
                resnames += [ 'solute' ] + [ self.labels[i]] * (self.n_monomers[i]-1)  

        #Make sure we didn't botch this
        assert len(resnames) == len( gromacs_topology.residues ), "Must have the same number of residues named as defined in the topology file."


        #Now we just go through and rename all the residues and we're done
        for i in range(len(resnames)):
            gromacs_topology.residues[i].name = resnames[i] 


        #Write GROMACS topology/coordinate files
        gromacs_topology.save(self.top_filename, format='gromacs')
        gromacs_topology.save(self.gro_filename)
