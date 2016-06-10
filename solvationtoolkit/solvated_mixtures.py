
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


# We require at least ParmEd 2.5.1 because of issues with the .mol2 writer 
# (issue #691 on ParmEd) prior to that, and 2.5.1.10 because of OpenEye reader
# formatting bugs requireing compressed spacing in .mol2 files (added in 
# ParmEd 2.5.1.10)
# Previously 2.0.4 or later was required due to issues with FudgeLJ/FudgeQQ in
# resulting GROMACS topologies in
#  earlier versions

try: #Try to get version tag
    ver = parmed.version
except: #If too old for version tag, it is too old
    oldParmEd = Exception('ERROR: ParmEd is too old, please upgrade to 2.5.1.10 or later (parmed-dev if necessary)' )
    raise oldParmEd
if ver < (2,5,1,10):
    raise RuntimeError("ParmEd is too old, please upgrade to 2.5.1.10 or later (parmed-dev if necessary)")


def make_path(pathname):
    try:
        os.makedirs(pathname)
    except:
        pass

class MixtureSystem(object):
    """A MixtureSystem object for defining and preparing input files for a 
mixture of arbitrary organic solutes. Outputs to formats for several simulation
packages are available.

    Usage
    -----------


    Parameters
    -----------
    directory : str (optional)
        Directory tree in which to store the data. Default: "data" 
    

    Limitations
    -----------
    Existing files with the same name present in the data directory tree may be
    overwritten. This results in a limitation/failure in a small (and probably
    random) fraction of cases if multiple systems involving the same monomers
    are written into the same data directory. Specifically, 
    openmoltools.amber.build_mixture_prmtop requires that each mol2 file for a
    component have a unique residue name, which is handled automatically by 
    openmoltools when constructing monomers (each is assigned a unique random
    residue name). However, if these are overwritten with other monomers (i.e.
    if we set up, say, 'octanol' in the same directory twice) 
    which by chance end up with non-unique residue names then 
    amber.build_mixture_prmtop will fail with a ValueError. This can be avoided
    by ensuring that if you are constructing multiple MixtureSystems involving
    the same monomers, your data directories are different. This issue also
    will likely be fixed when openmoltools switches to topology merging 
    via ParmEd rather than tleap, as unique residue names are built into 
    ParmEd in a better way. 
    """

                
    def __init__(self, directory='data'):
        
        """
        Initialization of the Molecule Database Class
    
        Parameters
        ----------
        directory : str 
           the directory name used to save the data

        """
        
        # Set directory names 
        self.data_path = directory
        self.data_path_monomers = os.path.join(self.data_path,'monomers')
        self.data_path_packmol = os.path.join(self.data_path,'packmol_boxes')

        # List container of all the added components to the solution
        self.component_list = []

        # Index used to perform index selection by using __iter__ function
        self.__ci = 0

        return


    def __str__(self):
        """
        Printing object function
        """

        
        string = ''

        for i in self.component_list:
            string = string + str(i)
        
        return string

    
    def __iter__(self):
        """
        Index generator
        """
        return self

    
    def next(self): # Python 3: def __next__(self)
        """
        Select the molecule during an iteration
        """
        
        if self.__ci > len(self.component_list) - 1:
            self.__ci = 0
            raise StopIteration
        else:
            self.__ci = self.__ci + 1
            return self.component_list[self.__ci - 1] 
           
        
    def __getitem__(self, index):
        """
        Index selection function
        """
        
        return self.component_list[index]
    

    def __setitem__(self, index, component):
        """
        Index setting function
        
        Parameters
        ----------
        index : int 
           the component index
        component : Component obj
           the component to assign to the component in the mixture
           MixtureSystem[index] = component
        
        """
        
        if not isinstance(component, Component):
            raise ValueError('The passed component is not a Component class object')
        
        self.component_list[index] = component


    
    def addComponent(self, name=None, **args):
        """
        Add a component to the solution 
        
        Parameters
        ----------
        name : string 
           the name of the compound to add the solution
        **args : see class Component for a full description
        
        """

        # Component object creation
        component=Component(name, **args)

        # Add object to the component list
        self.component_list.append(component)
        

        
    def build(self, amber=False, gromacs=False, solute_index='auto'):
        """
        Build all the monomers and the amber or gromacs mixture files 
        
        Parameters
        ----------
        amber : bool 
           this flag is used to control if output or not the amber files
        gromacs : bool
           this flag is used to control if output or not the gromacs files
        solute_index : int/str, optional. Default: "auto"
           Optional parameter to specify which of the components (in the list of specified components)
           will be treated as a solute in constructing GROMACS topology files 
           (which means that a single molecule of this component will be singled out as the 'solute'
           in the resulting GROMACS topology file). Valid options are 'auto' 
           (pick the first component present with n_monomers = 1,
           otherwise the first component), None (don't pick any), or an integer 
           (pick the component smiles_strings[solute_index].
     
        """

        #If no components were specified, then we can't proceed
        if len(self.component_list)==0: 
            raise TypeError("One or more components must be specified via addComponent")

        #If we want GROMACS, we have to also build AMBER as we get GROMACS from AMBER
        if gromacs and not amber:
            print("Turning on build of AMBER parameter/coordinate files in order to obtain files for GROMACS via conversion.")
            amber = True

        #NOW GENERATE FINAL STORAGE FOR FILE NAMES/COMPONENT LISTS
        #Now that we're building, we can generate the full set of components
        #we will be using since we're done adding components
        # List of all the smiles strings
        self.smiles_strings = []

        # List of all the number of monomers
        self.n_monomers = []

        # List of all the mole fractions
        self.mole_fractions = []

        # List of all the effective compound names. If the compond name is None
        # than the compound label will be used in this list as compound name 
        self.labels = []

        # The filling compound is a compound with None molecule number and None
        # mole fraction. It is used to fill out the solution
        self.filling_compound = None

        # Lists of filenames related to gaff mol2 files, amber files and sdf 
        # file format
        self.gaff_mol2_filenames = []
        self.frcmod_filenames = []
        self.inpcrd_filenames = [] 
        self.prmtop_filenames = []
        self.sdf_filenames = []

        # Final strings for output filenames 
        self.mix_fname = '' 
        self.pdb_filename = ''
        self.prmtop_filename = ''
        self.inpcrd_filename = ''
        self.top_filename = ''
        self.gro_filename = ''
        

        #BUILD
    

        # Now begin building by building monomers
        def build_monomers(self):
            """
            Generate GAFF mol2 and frcmod files for each chemical
            """

            # Filenames generation
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

                # Filling compound selection
                if comp.number == None and comp.mole_fraction == None:
                    if self.filling_compound == None:
                        self.filling_compound = comp
                        self.mole_fractions.append(comp.mole_fraction)
                    else:
                        raise ValueError('Error: Two or more fillig compounds have been specified')

                # Number and mol fractions lists generation
                if comp.number:
                    self.n_monomers.append(comp.number)
                if comp.mole_fraction is not None:
                    self.mole_fractions.append(comp.mole_fraction)

                # Lists of filenames generation    
                self.smiles_strings.append(comp.smiles)
                self.gaff_mol2_filenames.append(mol2_filename)
                self.frcmod_filenames.append(frcmod_filename)
                self.inpcrd_filenames.append(inpcrd_filename)
                self.prmtop_filenames.append(prmtop_filename)
                self.sdf_filenames.append(sdf_filename)

                if not (os.path.exists(mol2_filename) and os.path.exists(frcmod_filename)):
                     #Convert SMILES strings to mol2 and frcmod files for antechamber
                     openmoltools.openeye.smiles_to_antechamber(comp.smiles, mol2_filename, frcmod_filename)
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
            """
            Build an initial box with packmol and use it to generate AMBER files
            """

            def mole_fractions_to_n_monomers(self, density= 1 * grams/milliliter, cutoff=12*angstrom):
                """
                This function is used to generate the number of molecules for 
                each compound in the solution from the mole fractions of each molecule. 
                
                Parameters
                ----------
                density : openmm units 
                   the solution density 
                cutoff : openmm units
                   the cutoff distance of the largest compound in the solution 
                
                Returns
                -------
                self.n_monomers : integer list
                   the list of molecule number for each compound in the solution

                size : float
                   the edge of the box volume 

                """
                # Calculate the maximum atomic distance in a molecule 
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
                    
                    return max_dist * angstrom 

                # The sum of all the mole fractions 
                sum_fractions = sum([i for i in self.mole_fractions if i != None])
                
                if sum_fractions > 1.0:
                    raise ValueError('Error: The total molar fraction is greater than 1.0')               
                
                if sum_fractions == 1.0 and self.filling_compound:
                     raise ValueError('Error: The total molar fraction is 1.0 and it is not possible to add any filling compound to the solution, but a filling compound was specified') 
                
                if sum_fractions < 1.0 and not self.filling_compound:
                    raise ValueError('Error: The total molar fraction is less than 1.0 and no filling compound (i.e. compound with unspecified mole fraction) is provided')

                if self.filling_compound:
                    self.filling_compound.mole_fraction = 1.0 - sum_fractions
                    self.mole_fractions = [i if i != None else (1.0 - sum_fractions) for i in self.mole_fractions]


                max_dist_mols = 0.0 * angstrom
                delta_volume = 0.0 * angstrom**3
                sum_wgt_frac = 0.0 * grams/mole

         
                for i in range(0, len(self.sdf_filenames)):
                    istream = oemolistream(self.sdf_filenames[i])#gaff_mol2_files give wrong wgt because not sybyl format!
                    mol = oechem.OEMol()
                    
                    if not OEReadMolecule(istream, mol):
                        raise IOError('Error: It was not possible to create the OpenEye molecule object by reading the file: %s' % self.gaff_mol2_filenames[i])
                    # Molecular weight
                    wgt = oechem.OECalculateMolecularWeight(mol) * grams/mole
                    
                    if self.component_list[i].mole_fraction == 0.0:
                        delta_volume = oechem.OECalculateMolecularWeight(mol) * angstrom**3

                    sum_wgt_frac = sum_wgt_frac + wgt * self.component_list[i].mole_fraction
                    
                    max_dist= max_dist_mol(mol) 

                    if max_dist  > max_dist_mols:
                         max_dist_mols = max_dist 

                                                  
                cube_length = ((max_dist_mols + 2*cutoff)**3 + delta_volume)**(1.0/3.0)   

                n_monomers = []


                # n_i = Volume * Density * mole_fraction_i/sum_j(wgt_j * mole_fraction_j)
                self.n_monomers = [int(round(AVOGADRO_CONSTANT_NA * comp.mole_fraction * density * cube_length**3 / sum_wgt_frac)) \
                                   if comp.mole_fraction !=0 else 1 for comp in self.component_list]


                return self.n_monomers, cube_length


            if not self.gaff_mol2_filenames:
                raise ValueError('The list of gaff mol2 molecules is empty')

            if self.n_monomers and self.mole_fractions:
                    print (self.n_monomers, self.mole_fractions)
                    raise ValueError('Error: For different compounds it is not possible to mix mole_fractions and number of molecules')


            # The solution has been specified by using number of molecules    
            if self.n_monomers:
                
                if self.filling_compound:
                    raise ValueError('Error: The filling compound cannot be mixed with components specified by defining the number of molecules')
 
                size = openmoltools.packmol.approximate_volume_by_density(self.smiles_strings, self.n_monomers)
                packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], self.n_monomers, box_size = size)
                
                self.labels = self.mix_fname[1:].split('_')
                self.mix_fname = self.mix_fname[1:]  + ''.join(['_'+str(i) for i in self.n_monomers])
                self.pdb_filename = os.path.join(self.data_path_packmol, self.mix_fname+'.pdb')
                packed_trj.save(self.pdb_filename)
                

            # The solutions has been specified by using mole fractions    
            elif self.mole_fractions:

                n_monomers, size = mole_fractions_to_n_monomers(self)
                
                # WARNING: The size estimated with the mole_to_n_monomers
                # function is underestimating the volume calculated by using
                # openmoltools and for now we are using this estimate. 
                # If the volume is underestimated, apparently Packmol struggles
                # to find convergence and introduces extra molecules
                # into the found best solutionx (bug?)

                size = openmoltools.packmol.approximate_volume_by_density(self.smiles_strings, self.n_monomers)
                packed_trj = openmoltools.packmol.pack_box([md.load(mol2) for mol2 in self.gaff_mol2_filenames], n_monomers, box_size = size)
                
                self.labels = self.mix_fname[1:].split('_')
                self.mix_fname = self.mix_fname[1:] +''.join(['_'+str(i) for i in self.mole_fractions if i is not None])
                self.pdb_filename = os.path.join(self.data_path_packmol, self.mix_fname+'.pdb')
                packed_trj.save(self.pdb_filename)

            return



        def convert_to_gromacs(self, solute_index):
            """From AMBER-format prmtop and crd files, generate final solvated
GROMACS topology and coordinate files. Ensure that the desired "solute" (as per
solute_index) has a single monomer treated via a unique residue name to allow
treatment as a solute separate from other residues of the same name (if 
desired). The solute will be given residue name "solute" Also, check to see if
there are "WAT" residues present, in which case tleap will have re-ordered
them to the end of the data file. If so, update data structures accordingly
and handle conversion appropriately. 

    Notes
    -----
    Currently, this function ensures that - after AMBER conversion reorders 
    water molecules with residue names 'WAT' to occur last in the resulting
    parameter/coordinate files - the internal data structures are updated to 
    have the correct order in the relevant lists (labels, smiles_strings,
    n_monomers). If for some reason GROMACS conversion were removed, these 
    would need to be updated elsewhere. (Probably this should be done anyway,
    as this is not really a GROMACS issue.)

            """
            # Read in AMBER format parameter/coordinate file and convert in gromacs
            gromacs_topology = parmed.load_file(self.prmtop_filename, self.inpcrd_filename )
     
            # Split the topology into components and check that we have the right number of components
            components = gromacs_topology.split()
            assert len(components)==len(self.n_monomers), "Number of monomers and number of components in the combined topology do not match." 


            ####    HANDLE ORDERING OF WATER   ####
            # Check if any of the residues is named "WAT". If it is, antechamber will potentially have re-ordered it from where it was (it places residues named "WAT" at the end) so it may no longer appear in the order in which we expect.
            resnames = [ components[i][0].residues[0].name for i in range(len(components)) ]
            wat_present = False

            
            # Manage presence of WAT residues and possible re-ordering
            if 'WAT' in resnames:

                # If there is a water present, then we MIGHT have re-ordering. 
                # Check smiles to find out where it was originally.
                wat_orig_index = self.smiles_strings.index('O')

                # Where is it now?
                wat_new_index = resnames.index('WAT')

                # Reordered? If so, we have to adjust the ordering of 
                # n_monomers, smiles_strings, labels, and potentially 
                # solute_index. Filenames will be preserved since these were 
                #already created
                if wat_orig_index != wat_new_index:
                    # tleap moves water to the end so if they aren't equal, we
                    # know where water will be...
                    self.n_monomers = self.n_monomers[0:wat_orig_index] + self.n_monomers[wat_orig_index+1:] + [self.n_monomers[wat_orig_index]] 
                    self.smiles_strings = self.smiles_strings[0:wat_orig_index] + self.smiles_strings[wat_orig_index+1:] + [self.smiles_strings[wat_orig_index]] 
                    self.labels = self.labels[0:wat_orig_index] + self.labels[wat_orig_index+1:] + [self.labels[wat_orig_index] ]
                    # Check solute_index and alter if needed
                if not solute_index=='auto' and not solute_index==None:
                    # Index unchanged if it's before the water
                    if solute_index < wat_orig_index:
                        pass
                    # If it is the water, now it is at the end
                    elif solute_index == wat_orig_index:
                        solute_index = len(self.n_monomers)-1
                    # If it was after the water, then it moved up one position
                    else:
                        solute_index -= 1
            ####    END HANDLING OF ORDERING OF WATER   ####
                
             
            # Figure out what we're treating as the solute (if anything)
            if solute_index=='auto':
                # Check which of the molecules is present in qty 1
                try:
                    solute_index = self.n_monomers.index(1)
                except ValueError:
                    # If none is present in qty 1, then use the first 
                    solute_index = 0
            
            # Check that the passed solute index is correct
            check_solute_indices = range(0,len(self.n_monomers))
            assert solute_index in check_solute_indices and isinstance(solute_index, int) or solute_index == None, "Solute index must be an element of the list: %s or None. The value passed is: %s" % (check_solute_indices, solute_index)
            
            # Now all we have to do is to change the name of the solute molecule (residue, in ParmEd) and ParmEd will automatically make it a new molecule on write.
            # To do this, first build a list of the residue names we want, by molecule
            resnames = [ ]
            for i in range(len(self.n_monomers)):
                # If this is not the solute, just keep what we had
                if i!=solute_index:
                    resnames += [ self.labels[i] ] * self.n_monomers[i] 
                    # If it is the solute, make the first residue be named solute and the rest what they were already
                else: 
                    resnames += [ 'solute' ] + [ self.labels[i]] * (self.n_monomers[i]-1)  

            # Make sure we didn't botch this
            assert len(resnames) == len( gromacs_topology.residues ), "Must have the same number of residues named as defined in the topology file."


            # Now we just go through and rename all the residues and we're done
            for i in range(len(resnames)):
                gromacs_topology.residues[i].name = resnames[i] 


            # Write GROMACS topology/coordinate files
            gromacs_topology.save(self.top_filename, format='gromacs')
            gromacs_topology.save(self.gro_filename)

            return
       

        # Create monomers and packmol directories
        make_path(os.path.join(self.data_path_monomers))
        make_path(os.path.join(self.data_path_packmol))

        # Call the monomers creation and packmol systems
        build_monomers(self)
        build_boxes(self)
        
        # Create amber files
        if amber:
            self.data_path_amber = os.path.join(self.data_path,'amber')
            make_path(os.path.join(self.data_path_amber))
            self.prmtop_filename = os.path.join(self.data_path_amber, self.mix_fname+'.prmtop')
            self.inpcrd_filename = os.path.join(self.data_path_amber, self.mix_fname+'.inpcrd')
            tleap_cmd = openmoltools.amber.build_mixture_prmtop(self.gaff_mol2_filenames, self.frcmod_filenames, self.pdb_filename, self.prmtop_filename, self.inpcrd_filename)

            # Create gromacs files
            if gromacs:
                self.data_path_gromacs = os.path.join(self.data_path,'gromacs')
                make_path(os.path.join(self.data_path_gromacs))
                self.top_filename = os.path.join(self.data_path_gromacs, self.mix_fname+'.top')
                self.gro_filename = os.path.join(self.data_path_gromacs, self.mix_fname+'.gro')
                convert_to_gromacs(self,solute_index)

                

#*************************
# Component Class
#*************************

class Component(object):
    """
    This Class is used to save the component parameters

    """

    
    def __init__(self, name=None, label=None, smiles=None, number=None, mole_fraction=None):
        """
        Initialization class function 
        
        Parameters
        ----------
        name : str
           the molecule name
        label : str
           the molecule label used to generates files
        smiles : str
           the molecule SMILES string
        number : int
           the number of molecule 
        mole_fraction : float
           molecular mole fraction

        REQUIRED: A name and/or label. If no SMILES is provided, the name will be used to generate SMILES and if no name is provided, the label will be used to attempt to generate SMILES. 
        """


        # Checking name and label

        ref_str = ''
        
        if not name and not label:
            raise ValueError("Error: No component parameters name or label"+
" have been provided for the component")

        if label:
            if not isinstance(label, str):
                raise ValueError("Error: The component label %s is not a string" % label)
            ref_str = label

        if name:
            if not isinstance(name, str):
                raise ValueError("Error: The component name %s is not a string" % name)
            ref_str = name    
          
        if label and not name:
            print('\nWARNING: component name not provided; label will be used as component name\n')

        # Checking smiles, molecule number and mole fraction    
            
        if smiles:
            if not isinstance(smiles, str):
                raise ValueError("Error: The SMILES % for the component %s is not a string" % (smiles, ref_str))
            #Check this is a valid SMILES string
            mol = OEMol()
            status = OEParseSmiles(mol, smiles)
            if not status:
                raise ValueError("Error: The SMILES %s for the component %s"
                 " cannot be processed by OEChem." % (smiles, ref_str) )

        if number is not None:
            if not isinstance(number, int):
                raise ValueError("Error: The molecule quantity %s for the component %s is not an integer" % (number, ref_str))
            if number < 1:
                raise ValueError("Error: The molecule quantity %s for the component %s must be a positive integer" % (number, ref_str))
            
        if mole_fraction:
            if not isinstance(mole_fraction, float):
                raise ValueError("Error: The mole fraction %s for the component %s is not a float number" % (mole_fraction, ref_str))
            if mole_fraction < 0.0:
                raise ValueError("Error: The mole fraction %s for the component %s must be positive" % (mole_fraction, ref_str))
            if mole_fraction > 1.0:
                raise ValueError("Error: The mole fraction %s for the component %s is greater than one" % (mole_fraction, ref_str))
        
        if number and mole_fraction:
            raise ValueError("Error: molecule number and mole fraction for the compound %s cannot be both specified" % ref_str)
                    
        if not smiles:
            mol = OEMol()
            if name:
                try:
                    OEParseIUPACName(mol, name)
                    smiles = OECreateIsoSmiString(mol)
                    #If smiles is empty, didn't parse correctly
                    if smiles == '':
                        raise ValueError("Error: The supplied name '%s' could not be parsed" % name)
                except:
                    raise ValueError("Error: The supplied name '%s' could not be parsed" % name)
            elif label:
                try:
                    OEParseIUPACName(mol, label)
                    smiles = OECreateIsoSmiString(mol)
                    if smiles == '':
                        raise ValueError("Error: The supplied name '%s' could not be parsed" % name)
                except:
                    raise ValueError("Error: The supplied label '%s' could not be parsed" % label)

        self.name = name
        self.label = label
        self.smiles = smiles
        self.number = number
        self.mole_fraction = mole_fraction
        
            
        return


    def __str__(self):
        """
        Printing object function
        """
        
        return "\nname = %s\nlabel =  %s\nsmiles = %s\nnumber = %s\nmole_frac = %s\n" \
        %(self.name, self.label, self.smiles, self.number, self.mole_fraction)


    
