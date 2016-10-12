"""
Created on 27/04/2013

@author: thom
"""
import unittest

from rdkit.Chem import AllChem as Chem

from kinetic_molecule import KineticMolecule
from semi_realistic_chemistry import SemiRealisticChemistry


class TestSemiRealisticChemistry(unittest.TestCase):

    def setUp(self):
        
        bond_energies = {
            'H1H': 104.2,
            'C1C': 83,
            'N1N': 38.4,
            'O1O': 35,
            'H1C': 99, 'C1H': 99,
            'H1N': 93, 'N1H': 93,
            'H1O': 111, 'O1H': 111,
            'C1N': 73, 'N1C': 73,
            'C1O': 85.5, 'O1C': 85.5,
            'N1O': 55, 'O1N': 55,
            'C2O': 185, 'O2C': 185,  # rough average of range
            'C2C': 146,
            'N2N': 149,
            'O2O': 119,
            'C2N': 147, 'N2C': 147,
            'N2O': 143, 'O2N': 143,
            'C3O': 258, 'O3C': 258,
            'C3C': 200,
            'N3N': 226,
            'C3N': 213, 'N3C': 213,
            'C4C': 200  # theoretically possible from valences, but in nature forms a C2C bond instead
        }
        
        self.chem = SemiRealisticChemistry(bond_energies=bond_energies)

    def tearDown(self):
        self.chem = None
        
    def testBondEnergy(self):
        # energy is energy REQUIRED => - means releases energy
        self.assertEqual(-38.4, self.chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), to_bond_type=1))  # create single bond
        self.assertEqual(-38.4, self.chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), from_bond_type=0, to_bond_type=1))  # create single bond
        self.assertEqual(38.4, self.chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), from_bond_type=1, to_bond_type=0))  # destroy single bond
        self.assertEqual(149 - 38.4, self.chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), from_bond_type=2, to_bond_type=1))  # from double to single
        self.assertEqual(38.4 - 149, self.chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), from_bond_type=1, to_bond_type=2))  # from single to double
        self.assertEqual(38.4, self.chem.get_bond_energy(Chem.Atom('N'), Chem.Atom('N'), from_bond_type=1))  # delete single bond

    def testGetBondPotential(self):
        mol = KineticMolecule('[CH2-2].[CH2-2]')
        self.assertEqual(4, self.chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = KineticMolecule('O')  # H2O
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(0, self.chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = KineticMolecule('O')  # H2O, implicit Hs
        self.assertEqual(3, mol.GetNumAtoms())  # implicit no longer...
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(0, self.chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = KineticMolecule('[H]')
        self.assertEqual(1, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(1, self.chem.get_bond_potential(mol.GetAtoms()[0]))

        mol = KineticMolecule('O=C=O')  # CO2 - bond potentials of zero all round (full octets)
        for atom in mol.GetAtoms():
            self.assertEqual(0, self.chem.get_bond_potential(atom))

        mol = KineticMolecule('[OH-]')  # = [H][O-] Hydroxl anion
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(1, mol.GetAtoms()[1].GetAtomicNum())
        self.assertEqual(2, self.chem.get_bond_potential(mol.GetAtoms()[0]))
        self.assertEqual(0, self.chem.get_bond_potential(mol.GetAtoms()[1]))

        mol = KineticMolecule('[H].[OH-]')
        self.assertEqual(1, mol.GetAtoms()[0].GetAtomicNum())  # the H in [OH-]
        self.assertEqual(1, self.chem.get_bond_potential(mol.GetAtoms()[0]))
        self.assertEqual(1, mol.GetAtoms()[2].GetAtomicNum())  # the H in [OH-]
        self.assertEqual(0, self.chem.get_bond_potential(mol.GetAtoms()[2]))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
