import unittest

from rdkit.Chem import AllChem as Chem
from semi_realistic_chemistry import SemiRealisticChemistry
from molecule import Molecule
from chem_molecule import  ChemMolecule
from reaction import Reaction
import bond_energies


class TestSemiRealisticChemistry(unittest.TestCase):

    def allDifferent(self, l):
        for a in l:
            for b in l:
                for c in a:
                    if c not in b:
                        return True
        return False

    def test_all_different(self):
        self.assertTrue(self.allDifferent([['a', 'b'], ['c', 'a']]))
        self.assertFalse(self.allDifferent([['a', 'b'], ['b', 'a']]))
        self.assertTrue(self.allDifferent([['a', 'b'], ['b', 'a', 'c']]))

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

        l = chem._get_change_options(Reaction(reactants=[ChemMolecule('O=C=O')]))
        self.assertEqual(4, len(l))
        for r in l:
            self.assertIsInstance(r, Reaction)
            self.assertIsInstance(r.get_reactants()[0], ChemMolecule)
            if len(r.get_products()) > 0:
                self.assertIsInstance(r.get_products()[0], ChemMolecule)

    def test_ion_enumerate(self):
        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)

        # 3 possibilities: bond from H+ to O- (tricky - needs ion manipulation, break bond between O and H-; no bond possible between H+ and H-!
        # This requires three mini-steps - first, recognition that the ion has free unbonded electrons, second that can use those to form
        # bond between O and H, and third, that O and H are in different components of the combined molecule
        # This only applies if two molecules - if only [OH-] without [H] or other molecule then don't have this option

        r = Reaction(reactants=[ChemMolecule('[H+].[OH-]')])
        options = chem.enumerate(r)
        self.assertEqual(3, len(options))

    def test_enumerate(self):
        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('O=C=O')]))
        self.assertEqual(6, len(options))  # four options: break and drop to single bond from O=C bonds

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('[C]'), ChemMolecule('[C]')]))
        self.assertEqual(3, len(options))  # 3 types of bond formation - single, double, triple
        self.assertTrue(self.allDifferent([[p.get_symbol() for p in r.products] for r in options]))

        self.assertEqual(3, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[C].[C]')]))))  # 3 types of bond formation - single, double, triple
        self.assertEqual(2, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[O].[O]')]))))  # oxygen ions...pretty rare in nature - single and double bonds
        self.assertEqual(1, len(chem.enumerate(Reaction(reactants=[ChemMolecule('[H].[O]')]))))  # oxygen ion and proton...pretty rare in nature

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('O')]))
        self.assertEqual(2, len(options))  # two options, both breaks of H bonds

        options = chem.enumerate(Reaction(reactants=[ChemMolecule('[OH-]')]))
        self.assertEqual(2, len(options))

        self.assertEqual(7, len(chem.enumerate(Reaction(reactants=[ChemMolecule('C=C')]))))  # five complete breaks, and one drop from double to single

    def test_change_options(self):
        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)
        mol = Chem.AddHs(Chem.MolFromSmiles('[H]O[H]'))
        new_mol = chem._change_bond(mol, 0, 1, 0)
        self.assertEqual(Chem.MolToSmiles(new_mol), '[H+].[H][O-]')

    def test_split(self):
        """
        _split assumes that in parameter has explicit Hs.
        """

        r = SemiRealisticChemistry._split(Chem.AddHs(Chem.MolFromSmiles('O')))
        r = SemiRealisticChemistry._split(Chem.AddHs(Chem.MolFromSmiles('[CH2-2].[CH2-2]')))
        self.assertEqual(2, len(r))
        self.assertEqual('[H][C-2][H]', r[0].get_symbol())
        self.assertEqual('[H][C-2][H]', r[1].get_symbol())
        self.assertGreater(r[0].mass, 0)
        self.assertGreater(r[1].mass, 0)

        m = Chem.AddHs(Chem.MolFromSmiles('O=[N+]([O-])[N+](=O)[O-].[H]O[H]'))
        r = SemiRealisticChemistry._split(m)
        self.assertEqual(2, len(r))
        symbols = [x.get_symbol() for x in r]
        self.assertIn('O=[N+]([O-])[N+](=O)[O-]', symbols)
        self.assertIn('[H]O[H]', symbols)
        self.assertGreater(r[0].mass, 0)
        self.assertGreater(r[1].mass, 0)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
