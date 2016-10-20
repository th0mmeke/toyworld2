import unittest

from rdkit.Chem import AllChem as Chem
from semi_realistic_chemistry import SemiRealisticChemistry
from molecule import Molecule
from chem_molecule import  ChemMolecule
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
        l = chem._get_change_options(Reaction(reactants=[ChemMolecule('C')]))
        self.assertIsInstance(l, list)
        self.assertEqual(4, len(l))

        l = chem._get_change_options(Reaction(reactants=[ChemMolecule('O=C=O')]))
        for r in l:
            self.assertIsInstance(r, Reaction)
            self.assertIsInstance(r.get_reactants()[0], ChemMolecule)
            if len(r.get_products()) > 0:
                self.assertIsInstance(r.get_products()[0], ChemMolecule)

    def test_enumerate(self):

        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('O=C=O')]))
        self.assertEqual(4, len(options))  # four options: break and drop to single bond from O=C bonds

        # This requires three mini-steps - first, recognition that the ion has free unbonded electrons, second that can use those to form
        # bond between O and H, and third, that O and H are in different components of the combined molecule
        # This only applies if two molecules - if only [OH-] without [H] or other molecule then don't have this option
        options = chem.enumerate(Reaction(reactants=[ChemMolecule('[H+].[OH-]')]))
        self.assertEqual(2, len(options))  # bond from H+ to O- (tricky - needs ion manipulation), break bond between O and H-; no bond possible between H+ and H-!

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('O')]))
        self.assertEqual(2, len(options))  # two options, both breaks of H bonds

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('[OH-]')]))
        self.assertEqual(1, len(options))  # break H bond

        self.assertEqual(3, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[C].[C]')]))))  # 3 types of bond formation - single, double, triple
        self.assertEqual(3, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[C].[C]')]))))  # 3 types of bond formation - single, double, triple

        self.assertEqual(6, len(chem.enumerate(Reaction(reactants=[ChemMolecule('C=C')]))))  # five complete breaks, and one drop from double to single
        self.assertEqual(2, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[O].[O]')]))))  # oxygen ions...pretty rare in nature - single and double bonds
        self.assertEqual(1, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[H].[O]')]))))  # oxygen ion and proton...pretty rare in nature

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
