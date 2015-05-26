import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
from solvated_mixtures import MixtureSystem
import unittest
from openmoltools import utils

class TestMixtureSystem(unittest.TestCase):
    def setUp(self):
        with utils.enter_temp_directory(): 
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'test/data/',solute_index=2)
    #Test class Inizialitazion
    def test_InsufficientArgs(self):
        with utils.enter_temp_directory():
            #Check wrong number of argumetns
            self.assertRaises(TypeError,self.inst.__init__)
            self.assertRaises(TypeError,self.inst.__init__,['toluene','benzene'],['Cc1ccccc1','c1ccccc1'])
    def test_TypeArgs(self):
        #Check passed input Types
        with utils.enter_temp_directory():
            self.assertRaises(TypeError,self.inst.__init__,[1,'benzene'],['Cc1ccccc1','c1ccccc1'],[3,5],'data/', solute_index=0)
            self.assertRaises(TypeError,self.inst.__init__,['toluene','benzene'],['Cc1ccccc1','c1ccccc1'],[3.5,5],'data/', solute_index=0)
            self.assertRaises(AssertionError,self.inst.__init__,['toluene','benzene'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'test/data/', solute_index=2)
            self.assertRaises(TypeError,self.inst.__init__,['toluene','benzene'],['Cc1ccccc1','c1ccccc1'],[3,5],-3.9, solute_index=0)
    #Test class methods
    def test_build_monomers(self):
        with utils.enter_temp_directory():
            #Check smile strings
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['smile_error','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'data/', solute_index=2)
            self.assertRaises(ValueError,self.inst.build_monomers)
            #Check IO Errors
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'/', solute_index=2)
            self.assertRaises(IOError,self.inst.build_monomers)
    def test_convert_to_gromacs(self):
            #Check convert_via_acpype by using a wrong directory
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'/', solute_index=2)
            self.assertRaises(AttributeError,self.inst.convert_to_gromacs)
            #Check merge_topologies by using wrong filenames (solute_index=None case)
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'test/data/', solute_index=None)
            self.inst.top_filenames = ['test/data/Error_filenames']
            self.inst.top_filename = 'test/data/Error_filename'
            self.assertRaises(AssertionError,self.inst.convert_to_gromacs)
            #Check solute_index is in the correct range
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'test/data/', solute_index=6)
            self.assertRaises(AssertionError,self.inst.convert_to_gromacs)
            #Check merge_topologies by using wrong filenames (solute_index='auto' case)
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,1,80,7],'test/data/', solute_index='auto')
            self.inst.top_filenames = ['test/data/Error_filenames']
            self.inst.top_filename = 'test/data/Error_filename'
            self.assertRaises(AssertionError,self.inst.convert_to_gromacs)
    def test_build_boxes(self):
        #Check IO Errors
        with utils.enter_temp_directory():
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'/',solute_index=2)
            self.assertRaises(IOError,self.inst.build_boxes)
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'test/data/',solute_index=2)
            self.inst.inpcrd_filename = 'inpcrd_filename' 
            self.inst.prmtop_filename = 'prmtop_filename'
            self.assertRaises(IOError,self.inst.build_boxes)
    #Test run
    def test_run(self):
        with utils.enter_temp_directory():
            self.inst = MixtureSystem(['toluene','benzene','cyclohexane','ethane'],['Cc1ccccc1','c1ccccc1','C1CCCCC1','CC'],[3,5,80,7],'test/data/',solute_index=2)
            self.inst.run( just_build = True )
            self.inst.convert_to_gromacs()
if __name__ =='__main__':
    unittest.main()
    