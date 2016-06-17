import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
from openmoltools import utils
import unittest
from unittest import skipIf
from solvationtoolkit.solvated_mixtures import MixtureSystem


try:
    oechem = utils.import_("openeye.oechem")
    if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
    oequacpac = utils.import_("openeye.oequacpac")
    if not oequacpac.OEQuacPacIsLicensed(): raise(ImportError("Need License for oequacpac!"))
    oeiupac = utils.import_("openeye.oeiupac")
    if not oeiupac.OEIUPACIsLicensed(): raise(ImportError("Need License for OEOmega!"))        
    oeomega = utils.import_("openeye.oeomega")
    if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))    
    HAVE_OE = True
except:
    HAVE_OE = False


class TestMixtureSystem(unittest.TestCase):
    @skipIf(not HAVE_OE, "Cannot test core functionality without OpenEye tools.")
    def setUp(self):
        with utils.enter_temp_directory(): 
            self.inst = MixtureSystem()

    #Test class Initialization
    #These should be able to run without OpenEye tools at least with suitable re-architecting
    @skipIf(not HAVE_OE, "Cannot test core functionality without OpenEye tools.")
    def test_InsufficientInit(self):
        with utils.enter_temp_directory():
            #Check wrong number of arguments for adding a component - it requires a name or label at least.
            self.assertRaises(ValueError, self.inst.addComponent, smiles="CC")

            #Check what happens if we don't actually add components
            self.assertRaises(TypeError, self.inst.build )
            

    #Should be able to run without OE tools with suitable re-architecting    
    @skipIf(not HAVE_OE, "Cannot test core functionality without OpenEye tools.")
    def test_TypeArgs(self):
        #Check passed input Types
        with utils.enter_temp_directory():

            #Add a component with an integer name to ensure we catch
            self.inst = MixtureSystem()
            self.assertRaises(ValueError, self.inst.addComponent, name=1 )

            #Add a non-integer number of molecules to ensure we catch
            self.assertRaises(ValueError, self.inst.addComponent, name='phenol', number=3.5)

            #Add an invalid mole fraction
            self.assertRaises(ValueError, self.inst.addComponent, name='phenol', mole_fraction=1.2 )            
            self.assertRaises(ValueError, self.inst.addComponent, name='phenol', mole_fraction=-0.2 )            
            
            #Add mole fractions totaling greater than 1, check that we catch
            self.inst.addComponent('phenol', mole_fraction = 0.9)
            self.inst.addComponent('toluene', mole_fraction = 0.2)
            self.assertRaises( ValueError, self.inst.build )

    @skipIf(not HAVE_OE, "Cannot check smiles or name handling without OpenEye tools.")
    def test_ChemParsing(self):
        with utils.enter_temp_directory():
            self.inst = MixtureSystem()
            #Try passing invalid SMILES
            self.assertRaises( ValueError, self.inst.addComponent, name='phenol', smiles='smiles')
            #Try building with an invalid solute_index
            self.inst = MixtureSystem()
            self.inst.addComponent('toluene')
            self.assertRaises( AssertionError, self.inst.build, gromacs = True,
                solute_index = 2)
    

    #Test a bunch of different run cases which actually ought to work and ensure that they do
    @skipIf(not HAVE_OE, "Cannot test core functionality without OpenEye tools.")
    def test_run(self):
        with utils.enter_temp_directory():

            #Set up some names, labels, components
            names = ['toluene','benzene','cyclohexane','water', 'ethane']
            labels = ['toluene','benzene','cyclohexane','water', 'ethane']
            smiles = ['Cc1ccccc1','c1ccccc1','C1CCCCC1','O', 'CC'] 
            datapath = 'test' #Use non-default
            numbers = [3, 5, 80, 11, 7 ]
            mole_fractions = [ 0.1, 0.1, 0.1, 0.0, 0.7 ]
            n_components = len(names) 

            #Build using name
            self.inst = MixtureSystem( datapath )
            for n in range(n_components):
                self.inst.addComponent( name = names[n], mole_fraction = mole_fractions[n])
            self.inst.build()

            #Build using label
            self.inst = MixtureSystem( 'data' )
            for n in range(n_components):
                self.inst.addComponent( label = labels[n], mole_fraction = mole_fractions[n])
            self.inst.build()

            #Build using label and name
            self.inst = MixtureSystem( 'data2' )
            for n in range(n_components):
                self.inst.addComponent( name = names[n], label = labels[n], mole_fraction = mole_fractions[n])
            self.inst.build()

            #Build using label and smiles
            self.inst = MixtureSystem( 'data3' )
            for n in range(n_components):
                self.inst.addComponent( label = labels[n], smiles = smiles[n], mole_fraction = mole_fractions[n])
            self.inst.build()

            #Build using number rather than mole fraction
            self.inst = MixtureSystem( 'data4' )
            for n in range(n_components):
                self.inst.addComponent( label = labels[n], smiles = smiles[n], number = numbers[n])
            self.inst.build()
            

            #Build using solute_index and GROMACS (implying also AMBER)
            self.inst = MixtureSystem( 'data5' )
            for n in range(n_components):
                self.inst.addComponent( label = labels[n], smiles = smiles[n], number = numbers[n])
            self.inst.build( gromacs = True, solute_index = 2)

            #Build for AMBER only
            self.inst = MixtureSystem( 'data6' )
            for n in range(n_components):
                self.inst.addComponent( label = labels[n], smiles = smiles[n], number = numbers[n])
            self.inst.build( amber = True)


            #Make sure filling compound works
            self.inst = MixtureSystem( 'data7' )
            for n in range(n_components)[:-1]:
                self.inst.addComponent( label = labels[n], smiles = smiles[n], 
                mole_fraction = mole_fractions[n])
            self.inst.addComponent('methane')
            self.inst.build( gromacs = True )


            #We are already testing the infinite dilution case for one of the compounds


if __name__ =='__main__':
    unittest.main()
    
