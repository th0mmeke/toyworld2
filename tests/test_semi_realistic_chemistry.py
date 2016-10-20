import unittest

from rdkit.Chem import AllChem as Chem
from semi_realistic_chemistry import SemiRealisticChemistry
from molecule import Molecule
from reaction import Reaction
import bond_energies


class TestSemiRealisticChemistry(unittest.TestCase):

    def test_get_bond_energy(self):
        # energy is energy REQUIRED => - means releases energy
        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)
        self.assertEqual(-38.4, chem._get_bond_energy('N', 'N', to_bond_type=1))  # create single bond
        self.assertEqual(-38.4, chem._get_bond_energy('N', 'N', from_bond_type=0, to_bond_type=1))  # create single bond
        self.assertEqual(38.4, chem._get_bond_energy('N', 'N', from_bond_type=1, to_bond_type=0))  # destroy single bond
        self.assertEqual(149 - 38.4, chem._get_bond_energy('N', 'N', from_bond_type=2, to_bond_type=1))  # from double to single
        self.assertEqual(38.4 - 149, chem._get_bond_energy('N', 'N', from_bond_type=1, to_bond_type=2))  # from single to double
        self.assertEqual(38.4, chem._get_bond_energy('N', 'N', from_bond_type=1))  # delete single bond

    def test_change_options(self):
        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)
        l = chem._get_change_options(Reaction(reactants=[Molecule('C')]))
        self.assertIsInstance(l, list)
        self.assertEqual(4, len(l))

        l = chem._get_change_options(Reaction(reactants=[Molecule('O=C=O')]))
        for r in l:
            self.assertIsInstance(r, Reaction)
            self.assertIsInstance(r.get_reactants()[0], Molecule)
            if len(r.get_products()) > 0:
                self.assertIsInstance(r.get_products()[0], Molecule)

    def test_split(self):
        r = SemiRealisticChemistry._split(Chem.MolFromSmiles('O'))
        self.assertEqual('[H]O[H]', r[0].get_symbol())
        self.assertGreater(r[0].mass, 0)
        r = SemiRealisticChemistry._split(Chem.MolFromSmiles('[CH2-2].[CH2-2]'))
        self.assertEqual(2, len(r))
        self.assertEqual('[H][C-2][H]', r[0].get_symbol())
        self.assertEqual('[H][C-2][H]', r[1].get_symbol())
        self.assertGreater(r[0].mass, 0)
        self.assertGreater(r[1].mass, 0)

    def test_get_bond_potential(self):
        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)

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
