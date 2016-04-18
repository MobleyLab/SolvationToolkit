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
from openeye.oechem import *
import openeye.oequacpac as oequacpac


def writeSDF(mol2_filename, sdf_filename, mol_name):
    """For generating .sdf file format from .mol2 file
    
    Parameters
    ----------
    mol2_filename: str
    sdf_filename: str
    mol_name: str
        
    Limitations
    -----------
    Creates only two tags in .sdf file. The tags are partial_charges and partial_bond_orders. Their values are set according to the correspondent .mol2 file properties.
    """
    
    ifs = oemolistream(mol2_filename)
    MOL2flavor = OEIFlavor_Generic_Default | OEIFlavor_MOL2_Default | OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor(OEFormat_MOL2, MOL2flavor)

    ofs = oemolostream(sdf_filename)
    tag_names = ['partial_charges', 'partial_bond_orders']
    
    # Set the title and assign partial charges for the .mol2 file
    for mol in ifs.GetOEGraphMols():
        mol.SetTitle(mol_name)
        molToCharge = OEMol(mol)
        oequacpac.OEAssignPartialCharges( molToCharge, oequacpac.OECharges_AM1BCCSym )
        
        # Set partial charges to .mol2 file and store them to put into the .sdf file
        charges = []
        for atom, atomCharged in zip(mol.GetAtoms(), molToCharge.GetAtoms()):
            atom.SetPartialCharge( atomCharged.GetPartialCharge() )
            charges += [atom.GetPartialCharge()]
        
        # Create the tags for the sdf file
        mol = createTag(tag_names, mol, charges)
        OEWriteMolecule(ofs, mol)
    ifs.close()
    ofs.close()


def createTag(tag_names, mol, charges):
    """Create the tags for sdf file."""
    mol_id = OEGetSDData(mol, 'Mol_Index')
    for tag in tag_names:
        if tag == 'partial_charges':
            value = manipulatePartialChargesTag(charges)
        elif tag == 'partial_bond_orders':
            value = manipulateBondOrdersTag(mol)
        else:
            #Tag Name Error
            raise Exception('Wrong tag name')
            break
        OESetSDData(mol, OESDDataPair(tag, value))
    return mol


def manipulatePartialChargesTag(charges):
    """Transform charges from float to string format."""
    value = ''
    print "Charges: %s" % str(charges)
    for charge in charges:
        value += str(charge) + '\n'
    return value


def manipulateBondOrdersTag(mol):
    """Get bonds and store their order for partial bond orders Tag"""
    partial_bonds_order = ''
    for bond in mol.GetBonds():
        partial_bonds_order += str(float(bond.GetOrder())) + '\n'
    print 'Partial_bonds_order: ' + partial_bonds_order
    return partial_bonds_order
