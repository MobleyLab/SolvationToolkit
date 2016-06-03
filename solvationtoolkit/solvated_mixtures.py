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
import solvationtoolkit.mol2tosdf as mol2tosdf
from openeye.oechem import *
from openeye.oeiupac import *
from simtk.unit import * 


# We require at least ParmEd 2.5.1 because of issues with the .mol2 writer (issue #691 on ParmEd) prior to that, and 2.5.1.10 because of OpenEye reader formatting bugs requireing compressed spacing in .mol2 files (added in ParmEd 2.5.1.10)
# Previously 2.0.4 or later was required due to issues with FudgeLJ/FudgeQQ in resulting GROMACS topologies in
#  earlier versions
try: #Try to get version tag
    ver = parmed.version
except: #If too old for version tag, it is too old
    oldParmEd = Exception('ERROR: ParmEd is too old, please upgrade to 2.5.1 or later')
    raise oldParmEd
if ver < (2,5,1,10):
    raise RuntimeError("ParmEd is too old, please upgrade to 2.5.1 or later")


def make_path(pathname):
    try:
        os.makedirs(pathname)
    except:
        raise IOError('It was not possible to create the directory %s'%pathname)


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

    class Component(object):

        def __init__(self, name=None, label=None, smile=None, numbers=None, mole_fraction=None):
            
            # Checking passed types
            if name:
                if not isinstance(name, str):
                    raise ValueError("Error: The component parameter name %s is not a string" % name)
            if label:
                if not isinstance(label, str):
                    raise ValueError("Error: The component parameter label %s is not a string" % label)
            if smile:
                if not isinstance(smile, str):
                    raise ValueError("Error: The component parameter smile %s is not a string" % smile)
                #TO DO: Check if a string is a valid smile string
            if numbers:
                if not isinstance(numbers, int):
                    raise ValueError("Error: The component parameter numbers %s is not an integer" % numbers)
                if numbers < 1:
                     raise ValueError("Error: The selected number of molecule %s must be a positive integer" % numbers)
            if mole_fraction:
                if not isinstance(mole_fraction, float):
                    raise ValueError("Error: The component parameter mole_fraction %s is not a float" % mole_fraction)
                if mole_fraction < 0.0:
                    raise ValueError("Error: The selected mole_fraction %s must be a positive integer" % mole_fraction)

            # Checking  name and label 
            if not name and not label:
                raise ValueError("Error: No component parameters name or label has been provided")

            if numbers and mole_fraction:
                if name:
                    raise ValueError("Error: numbers and mole fraction of compound name %s cannot be both specified" % name)
                elif label:
                    raise ValueError("Error: numbers and mole fraction of compound label %s cannot be both specified" % label)
                    
            if not smile:
                mol = OEMol()
                if name:
                    try:
                        OEParseIUPACName(mol, name)
                        smile = OECreateIsoSmiString(mol)
                    except:
                        raise ValueError("Error: The supplied name '%s' could not be parsed" % name)
                elif label:
                    try:
                        OEParseIUPACName(mol, label)
                        smile = OECreateIsoSmiString(mol)
                    except:
                        raise ValueError("Error: The supplied label '%s' could not be parsed" % label)

            self.name = name
            self.label = label
            self.smile = smile
            self.numbers = numbers
            self.mole_fraction = mole_fraction
            self.properties = {}
            self.openeye_mol = None
            
            return


        def __str__(self):
            return "name = %s\nlabel =  %s\nsmile = %s\nnumbers = %s\nmole_frac = %s\nproperties = %s\n" \
            %(self.name, self.label, self.smile, self.numbers, self.mole_fraction, self.properties)
                

    def __init__(self, directory='data'):
        
        self.data_path = directory
        self.data_path_monomers = os.path.join(self.data_path,'monomers')
        self.data_path_packmol = os.path.join(self.data_path,'packmol_boxes')
        self.data_path_amber = os.path.join(self.data_path,'amber')
        self.data_path_gromacs = os.path.join(self.data_path,'gromacs')
        
        self.component_list = []

        self.smile_strings = []
        self.n_monomers = []
        self.mole_fractions = []
        self.labels = []
        self.filling_compound = None

        self.gaff_mol2_filenames = []
        self.frcmod_filenames = []
        self.inpcrd_filenames = [] 
        self.prmtop_filenames = []
        self.sdf_filenames = []
        
        
        self.mix_fname = '' 
        self.pdb_filename = ''
        self.prmtop_filename = ''
        self.inpcrd_filename = ''
        self.top_filename = ''
        self.gro_filename = ''

        return

    def addComponent(self, name=None, **args):
    
        component=self.Component(name, **args)
        
        self.component_list.append(component)
        

    def build(self, amber=False, gromacs=False, solute_index='auto'):
        
        def build_monomers(self):
            """Generate GAFF mol2 and frcmod files for each chemical."""
            
            for comp in self.component_list:
                if comp.label:
                    mol2_filename = os.path.join(self.data_path_monomers, comp.label+'.mol2')
                    frcmod_filename = os.path.join(self.data_path_monomers, comp.label+'.frcmod')
                    inpcrd_filename  = os.path.join(self.data_path_monomers, comp.label+'.inpcrd')
                    prmtop_filename  = os.path.join(self.data_path_monomers, comp.label+'.prmtop')
                    sdf_filename = os.path.join(self.data_path_monomers, comp.label+'.sdf')
                    self.mix_fname = self.mix_fname + '_' + comp.label
                else:
                    mol2_filename = os.path.join(self.data_path_monomers, comp.name+'.mol2')
                    frcmod_filename = os.path.join(self.data_path_monomers, comp.name+'.frcmod')
                    inpcrd_filename  = os.path.join(self.data_path_monomers, comp.name+'.inpcrd')
                    prmtop_filename  = os.path.join(self.data_path_monomers, comp.name+'.prmtop')
                    sdf_filename = os.path.join(self.data_path_monomers, comp.name+'.sdf')
                    self.mix_fname = self.mix_fname + '_' + comp.name

                    
                if comp.numbers == None and comp.mole_fraction == None:
                    if self.filling_compound == None:
                        self.filling_compound = comp
                        self.mole_fractions.append(comp.mole_fraction)
                    else:
                        raise ValueError('Error: Two or more fillig compounds have been specified')
                    
                if comp.numbers:
                    self.n_monomers.append(comp.numbers)
                if comp.mole_fraction is not None:
                    self.mole_fractions.append(comp.mole_fraction)


                self.smile_strings.append(comp.smile)
                self.gaff_mol2_filenames.append(mol2_filename)
                self.frcmod_filenames.append(frcmod_filename)
                self.inpcrd_filenames.append(inpcrd_filename)
                self.prmtop_filenames.append(prmtop_filename)
                self.sdf_filenames.append(sdf_filename)


                if not (os.path.exists(mol2_filename) and os.path.exists(frcmod_filename)):
                     #Convert SMILES strings to mol2 and frcmod files for antechamber
                     openmoltools.openeye.smiles_to_antechamber(comp.smile, mol2_filename, frcmod_filename)
                     #Correct the mol2 file partial atom charges to have a total net integer molecule charge  
                     mol2f = parmed.formats.Mol2File
                     mol2f.write(parmed.load_file(mol2_filename).fix_charges(),mol2_filename, compress_whitespace=True)
             
                #Generate amber coordinate and topology files for the unsolvated molecules
                mol_name = os.path.basename(mol2_filename).split('.')[0]
                openmoltools.amber.run_tleap(mol_name, mol2_filename, frcmod_filename, prmtop_filename, inpcrd_filename)
    
                #Read Mol2 File and write SDF file
                mol2tosdf.writeSDF(mol2_filename, sdf_filename, mol_name)


            #Generate unique residue names for molecules in mol2 files
            openmoltools.utils.randomize_mol2_residue_names(self.gaff_mol2_filenames) 


        def build_boxes(self):
            """Build an initial box with packmol and use it to generate AMBER files."""

            def mole_fractions_to_n_monomers(self, density=grams/milliliter, cutoff=12*angstrom):
                
                def max_dist_mol(mol):
                    max_dist = 0.0
                    coords = mol.GetCoords() # Are the coords always in A in mol2 file?
                    for i in range(0, mol.NumAtoms()):
                        crdi = np.array([coords[i][0], coords[i][1], coords[i][2]])
                        for j in range(i+1, mol.NumAtoms()):
                            crdj = np.array([coords[j][0], coords[j][1], coords[j][2]])
                            dist = np.linalg.norm(crdi-crdj)
                            if dist > max_dist:
                                max_dist = dist
                    
                    return max_dist # In angstrom
                
                sum_fractions = sum([i for i in self.mole_fractions if i != None])
                
                if sum_fractions > 1.0:
                    raise ValueError('Error: The total molar fraction is greater than 1.0')               
                
                if sum_fractions == 1.0 and self.filling_compound:
                     raise ValueError('Error: The total molar fraction is 1.0 and it is not possible to add any filling compound to the solution') 
                
                if sum_fractions <1.0 and not self.filling_compound:
                    raise ValueError('Error: The total molar fraction is less than 1.0 and the filling compoind is missing')

                if self.filling_compound:
                    self.filling_compound.mole_fraction = 1.0 - sum_fractions
                    self.mole_fractions = [i if i != None else (1.0 - sum_fractions) for i in self.mole_fractions]


                max_dist_mols = 0.0
                delta_volume = 0.0 * angstrom**3
                sum_wgt_frac = 0.0 * grams/mole

         
                for i in range(0, len(self.sdf_filenames)):
                    istream = oemolistream(self.sdf_filenames[i])#gaff_mol2_files give wrong wgt because not sybyl format!
                    mol = oechem.OEMol()
                    
                    if not OEReadMolecule(istream, mol):
                        raise IOError('Error: It was not possible to create the OpenEye molecule object reading the file: %s' % self.gaff_mol2_filenames[i])
                    
                    self.openeye_mol = mol

                    self.component_list[i].properties['wgt'] = oechem.OECalculateMolecularWeight(mol) * grams/mole
                    
                    if self.component_list[i].mole_fraction == 0.0:
                        delta_volume = oechem.OECalculateMolecularWeight(mol) * angstrom**3

                    sum_wgt_frac = sum_wgt_frac + self.component_list[i].properties['wgt'] * self.component_list[i].mole_fraction
                    
                    self.component_list[i].properties['max_dist'] = max_dist_mol(mol) 

                    if self.component_list[i].properties['max_dist']  > max_dist_mols:
                         max_dist_mols = self.component_list[i].properties['max_dist'] 

                                                  
                cube_length = ((max_dist_mols * angstrom + 2*cutoff)**3 + delta_volume)**(1.0/3.0)   

                n_monomers = []

                self.n_monomers = [int(round(AVOGADRO_CONSTANT_NA * comp.mole_fraction * density * cube_length**3 / sum_wgt_frac)) if comp.mole_fraction !=0 else 1 for comp in self.component_list]


                return self.n_monomers, cube_length


            if not self.gaff_mol2_filenames:
                raise ValueError('The list of gaff mol2 molecules is empty')

            if self.n_monomers and self.mole_fractions:
                    raise ValueError('Error: For different compounds it is not possible to mix mole_fractions and number of molecules')

            if self.n_monomers:
                
                if self.filling_compound:
                    raise ValueError('Error: The filling compound cannot be mixed with components specified by defining the number of molecules')
 
                size = openmoltools.packmol.approximate_volume_by_density(self.smile_strings, self.n_monomers)
                packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], self.n_monomers, box_size = size)
                
                self.labels = self.mix_fname[1:].split('_')
                self.mix_fname = self.mix_fname[1:]  + ''.join(['_'+str(i) for i in self.n_monomers])
                self.pdb_filename = os.path.join(self.data_path_packmol, self.mix_fname+'.pdb')
                packed_trj.save(self.pdb_filename)
                

            elif self.mole_fractions:

                n_monomers, size = mole_fractions_to_n_monomers(self)
                
                
                # The size estimated with the mole_to_n_monomers function is underestimating the volume calculated by using openmoltools and for now we are using this estimate
                # Packmol is struggling to find convergence and introduces extra molecules into the best solution found (bug?)
                size = openmoltools.packmol.approximate_volume_by_density(self.smile_strings, self.n_monomers)
                packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], n_monomers, box_size = size)
                #packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], n_monomers, box_size = size/anstrom)
                
                self.labels = self.mix_fname[1:].split('_')
                self.mix_fname = self.mix_fname[1:] +''.join(['_'+str(i) for i in self.mole_fractions if i is not None])
                self.pdb_filename = os.path.join(self.data_path_packmol, self.mix_fname+'.pdb')
                packed_trj.save(self.pdb_filename)

            return



        def convert_to_gromacs(self, solute_index):
            """From AMBER-format prmtop and crd files, generate final solvated GROMACS topology and coordinate files. Ensure that the desired "solute" (as per solute_index) has a single monomer treated via a unique residue name to allow treatment as a solute separate from other residues of the same name (if desired). The solute will be given residue name "solute" Also, check to see if there are "WAT" residues present, in which case tleap will have re-ordered them to the end of the data file. If so, update data structures accordingly and handle conversion appropriately. 

            Notes
            -----
            Currently, this function ensures that - after AMBER conversion reorders water molecules with residue names 'WAT' to occur last in the resulting parameter/coordinate files - the internal data structures are updated to have the correct order in the relevant lists (labels, smiles_strings, n_monomers). If for some reason GROMACS conversion were removed, these would need to be updated elsewhere. (Probably this should be done anyway, as this is not really a GROMACS issue.)

            """
            #Read in AMBER format parameter/coordinate file and convert in gromacs
            gromacs_topology = parmed.load_file(self.prmtop_filename, self.inpcrd_filename )
     
            #Split the topology into components and check that we have the right number of components
            components = gromacs_topology.split()
            assert len(components)==len(self.n_monomers), "Number of monomers and number of components in the combined topology do not match." 


            ####    HANDLE ORDERING OF WATER   ####
            #Check if any of the residues is named "WAT". If it is, antechamber will potentially have re-ordered it from where it was (it places residues named "WAT" at the end) so it may no longer appear in the order in which we expect.
            resnames = [ components[i][0].residues[0].name for i in range(len(components)) ]
            wat_present = False

            
            #Manage presence of WAT residues and possible re-ordering
            if 'WAT' in resnames:
                #If there is a water present, then we MIGHT have re-ordering. Check smiles to find out where it was originally.
                wat_orig_index = self.smile_strings.index('O')
                #Where is it now?
                wat_new_index = resnames.index('WAT')
                #Reordered? If so, we have to adjust the ordering of n_monomers, smiles_strings, labels,
                # and potentially solute_index. Filenames will be preserved since these were already created
                if wat_orig_index != wat_new_index:
                    #tleap moves water to the end so if they aren't equal, we know where water will be...
                    self.n_monomers = self.n_monomers[0:wat_orig_index] + self.n_monomers[wat_orig_index+1:] + [self.n_monomers[wat_orig_index]] 
                    self.smile_strings = self.smile_strings[0:wat_orig_index] + self.smile_strings[wat_orig_index+1:] + [self.smile_strings[wat_orig_index]] 
                    self.labels = self.labels[0:wat_orig_index] + self.labels[wat_orig_index+1:] + [self.labels[wat_orig_index] ]
                    #Check solute_index and alter if needed
                if not solute_index=='auto' and not solute_index==None:
                    #Index unchanged if it's before the water
                    if solute_index < wat_orig_index:
                        pass
                    #If it is the water, now it is at the end
                    elif solute_index == wat_orig_index:
                        solute_index = len(self.n_monomers)-1
                    #If it was after the water, then it moved up one position
                    else:
                        solute_index -= 1
            ####    END HANDLING OF ORDERING OF WATER   ####
                
             
            #Figure out what we're treating as the solute (if anything)
            if solute_index=='auto':
                #Check which of the molecules is present in qty 1
                try:
                    solute_index = self.n_monomers.index(1)
                except ValueError:
                    #If none is present in qty 1, then use the first 
                    solute_index = 0
            
            #Check that the passed solute index is correct
            check_solute_indices = range(0,len(self.n_monomers))
            assert solute_index in check_solute_indices and isinstance(solute_index, int) or solute_index == None, "Solute index must be an element of the list: %s or None. The value passed is: %s" % (check_solute_indices,self.solute_index)
            
            #Now all we have to do is to change the name of the solute molecule (residue, in ParmEd) and ParmEd will automatically make it a new molecule on write.
            #To do this, first build a list of the residue names we want, by molecule
            resnames = [ ]
            for i in range(len(self.n_monomers)):
                #If this is not the solute, just keep what we had
                if i!=solute_index:
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

            return



        # if os.path.exists(self.data_path):
        #     raise IOError('The directory %s already exists' % self.data_path)
            

        make_path(os.path.join(self.data_path_monomers))
        make_path(os.path.join(self.data_path_packmol))

        build_monomers(self)
        build_boxes(self)
        
        
        if amber:
            make_path(os.path.join(self.data_path_amber))
            self.prmtop_filename = os.path.join(self.data_path_amber, self.mix_fname+'.prmtop')
            self.inpcrd_filename = os.path.join(self.data_path_amber, self.mix_fname+'.inpcrd')
            tleap_cmd = openmoltools.amber.build_mixture_prmtop(self.gaff_mol2_filenames, self.frcmod_filenames, self.pdb_filename, self.prmtop_filename, self.inpcrd_filename)
            if gromacs:
                make_path(os.path.join(self.data_path_gromacs))
                self.top_filename = os.path.join(self.data_path_gromacs, self.mix_fname+'.top')
                self.gro_filename = os.path.join(self.data_path_gromacs, self.mix_fname+'.gro')
                convert_to_gromacs(self,solute_index)

                
