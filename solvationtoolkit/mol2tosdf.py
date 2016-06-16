######################################################################
# SolvationToolkit: A toolkit for setting up molecular simulations of mixtures
# Copyright 2011-2015 UC Irvine and the Authors
#
# Authors: David Mobley, Gaetano Calabro, and Camila Zanette
# With thanks to Kyle Beauchamp, whose liquid_tools.py provided an initial basis for this code
# (https://github.com/choderalab/LiquidBenchmark/blob/master/src/simulation/liquid_tools.py) in April 2015
# Thanks to Christopher Bayly, whose helped with reading .mol2 files with GAFF atom types.
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


import sys
from openmoltools.utils import import_

def writeSDF(mol2_filename, sdf_filename, mol_name):
    """For generating .sdf file format from .mol2 file, using OEIFlavor (OpenEye). Creates three tags (partial_charges, partial_bond_orders, and atom_types) in the .sdf file using values from the reference (.mol2 file).
    
    Parameters
    ----------
    mol2_filename: str
        Mol2 filename used to write the .sdf file. The .mol2 file is preserved.
    sdf_filename: str
        SDF filename (path) used to save the sdf file generated.
    mol_name: str
        Name of molecule used to write the title for the .sdf file generated.
    
    Notes
    -----------
    In partial_bonds_orders tag, .mol2 bond orders are translated to .sdf bond orders using OEIFlavor, aromatic bonds ('ar') are transformed to 1 or 2.

    Limitations
    -----------
    Creates only three tags in .sdf file. The tags are partial_charges, partial_bond_orders, and atom_types (GAFF). Their values are set according to the correspondent .mol2 file properties.
    """
    oechem = import_("openeye.oechem")
    
    ifs = oechem.oemolistream(mol2_filename)
    MOL2flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor(oechem.OEFormat_MOL2, MOL2flavor)

    ofs = oechem.oemolostream(sdf_filename)
    tag_names = ['partial_charges', 'partial_bond_orders', 'atom_types']
    
    # Set the title and assign partial charges for the .mol2 file
    for mol in ifs.GetOEGraphMols():
        mol.SetTitle(mol_name)
        molToCharge = oechem.OEMol(mol)
        
        # Get partial charges and atom types from .mol2 file and store them to put into the .sdf file as tags
        charges = []
        atom_types = []
        for atom, atomCharged in zip(mol.GetAtoms(), molToCharge.GetAtoms()):
            atom.SetPartialCharge( atomCharged.GetPartialCharge() )
            charges += [atom.GetPartialCharge()]
            atom_types += [atom.GetType()]

        #print("partial charges " + str(charges))
        #print("atom types: " + str(atom_types))
        # Create the tags for the sdf file
        mol = createTag(tag_names, mol, charges, atom_types)
        oechem.OEWriteMolecule(ofs, mol)
    ifs.close()
    ofs.close()


def createTag(tag_names, mol, charges, atom_types):
    """Create the tags for sdf file."""
    oechem = import_("openeye.oechem")
    mol_id = oechem.OEGetSDData(mol, 'Mol_Index')
    for tag in tag_names:
        if tag == 'partial_charges':
            value = manipulatePartialChargesTag(charges)
        elif tag == 'partial_bond_orders':
            value = manipulateBondOrdersTag(mol)
        elif tag == 'atom_types':
            value = manipulateAtomTypes(atom_types)
        else:
            #Tag Name Error
            raise Exception('Wrong tag name')
            break
        OESetSDData(mol, OESDDataPair(tag, value))
    return mol


def manipulatePartialChargesTag(charges):
    """Transform charges from float to string format."""
    value = ''
    for charge in charges:
        value += str(charge) + '\n'
    return value


def manipulateBondOrdersTag(mol):
    """Get bonds and store their order for partial bond orders Tag"""
    partial_bonds_order = ''
    for bond in mol.GetBonds():
        partial_bonds_order += str(float(bond.GetOrder())) + '\n'
    return partial_bonds_order

def manipulateAtomTypes(atom_types):
    """Transform atom types to the right string format."""
    value = ''
    for atm in atom_types:
        value += str(atm) + '\n'
    return value
