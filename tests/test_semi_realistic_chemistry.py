import unittest

from rdkit.Chem import AllChem as Chem
from semi_realistic_chemistry import SemiRealisticChemistry
from molecule import Molecule
from reaction import Reaction


class TestSemiRealisticChemistry(unittest.TestCase):

    def setUp(self):
        
        self.bond_energies = {
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

    def test_get_bond_energy(self):
        # energy is energy REQUIRED => - means releases energy
        chem = SemiRealisticChemistry(bond_energies=self.bond_energies)
        self.assertEqual(-38.4, chem._get_bond_energy('N', 'N', to_bond_type=1))  # create single bond
        self.assertEqual(-38.4, chem._get_bond_energy('N', 'N', from_bond_type=0, to_bond_type=1))  # create single bond
        self.assertEqual(38.4, chem._get_bond_energy('N', 'N', from_bond_type=1, to_bond_type=0))  # destroy single bond
        self.assertEqual(149 - 38.4, chem._get_bond_energy('N', 'N', from_bond_type=2, to_bond_type=1))  # from double to single
        self.assertEqual(38.4 - 149, chem._get_bond_energy('N', 'N', from_bond_type=1, to_bond_type=2))  # from single to double
        self.assertEqual(38.4, chem._get_bond_energy('N', 'N', from_bond_type=1))  # delete single bond

    def test_change_options(self):
        chem = SemiRealisticChemistry(bond_energies=self.bond_energies)
        l = chem.get_change_options(Chem.MolFromSmiles('C'))
        self.assertIsInstance(l, list)
        self.assertEqual(0, len(l))

        l = chem.get_change_options(Chem.MolFromSmiles('O=C=O'))
        for r in l:
            self.assertIsInstance(r, Reaction)
            self.assertIsInstance(r.get_reactants()[0], Molecule)
            if len(r.get_products()) > 0:
                self.assertIsInstance(r.get_products()[0], Molecule)

    def test_split(self):
        r = SemiRealisticChemistry.split(Chem.MolFromSmiles('O'))
        self.assertEqual('O', r[0].get_symbol())
        r = SemiRealisticChemistry.split(Chem.MolFromSmiles('[CH2-2].[CH2-2]'))
        self.assertEqual(2, len(r))
        self.assertEqual('[CH2-2]', r[0].get_symbol())
        self.assertEqual('[CH2-2]', r[1].get_symbol())

    def test_get_bond_potential(self):
        chem = SemiRealisticChemistry(bond_energies=self.bond_energies)

        mol = Chem.AddHs(Chem.MolFromSmiles('[CH2-2].[CH2-2]'))
        self.assertEqual(4, chem._get_bond_potential(mol.GetAtoms()[0]))

        mol = Chem.AddHs(Chem.MolFromSmiles('O'))  # H2O
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(0, chem._get_bond_potential(mol.GetAtoms()[0]))

        mol = Chem.AddHs(Chem.MolFromSmiles('O'))  # H2O, implicit Hs
        self.assertEqual(3, mol.GetNumAtoms())  # implicit no longer...
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(0, chem._get_bond_potential(mol.GetAtoms()[0]))

        mol = Chem.AddHs(Chem.MolFromSmiles('[H]'))
        self.assertEqual(1, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(1, chem._get_bond_potential(mol.GetAtoms()[0]))

        mol = Chem.AddHs(Chem.MolFromSmiles('O=C=O'))  # CO2 - bond potentials of zero all round (full octets)
        for atom in mol.GetAtoms():
            self.assertEqual(0, chem._get_bond_potential(atom))

        mol = Chem.AddHs(Chem.MolFromSmiles('[OH-]'))  # = [H][O-] Hydroxl anion
        self.assertEqual(8, mol.GetAtoms()[0].GetAtomicNum())
        self.assertEqual(1, mol.GetAtoms()[1].GetAtomicNum())
        self.assertEqual(2, chem._get_bond_potential(mol.GetAtoms()[0]))
        self.assertEqual(0, chem._get_bond_potential(mol.GetAtoms()[1]))

        mol = Chem.AddHs(Chem.MolFromSmiles('[H].[OH-]'))
        self.assertEqual(1, mol.GetAtoms()[0].GetAtomicNum())  # the H in [OH-]
        self.assertEqual(1, chem._get_bond_potential(mol.GetAtoms()[0]))
        self.assertEqual(1, mol.GetAtoms()[2].GetAtomicNum())  # the H in [OH-]
        self.assertEqual(0, chem._get_bond_potential(mol.GetAtoms()[2]))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
