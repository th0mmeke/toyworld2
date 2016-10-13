import unittest

from rdkit.Chem import AllChem as Chem
from rdkit.rdBase import DisableLog, EnableLog
from kinetic_molecule import KineticMolecule
from semi_realistic_chemistry import SemiRealisticChemistry


class Test(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')

    def tearDown(self):
        EnableLog('rdApp.error')

    def testInit(self):
        mol = KineticMolecule("O=C=O", internal_energy=30)
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(30, mol.get_internal_energy())

        mol = KineticMolecule(Chem.MolFromSmiles("O=C=O"), internal_energy=30)
        self.assertEqual(3, mol.GetNumAtoms())
        self.assertEqual(30, mol.get_internal_energy())
        self.assertEqual(44.009, mol.get_mass())

        with self.assertRaises(ValueError):
            KineticMolecule(Chem.MolFromSmiles("O"), -10)

        mol = KineticMolecule('[O][H-]')
        self.assertEqual('[H][O]', Chem.MolToSmiles(mol))  # RDKit converts to this...
        mol = KineticMolecule('[H+].[OH-]')
        self.assertEqual('[H+].[H][O-]', Chem.MolToSmiles(mol))  # RDKit converts to this...

    def testPotentialEnergy(self):
        bond_energies = {'O1O': 20, 'O2O': 30, 'O1H': 40}
        chem = SemiRealisticChemistry(bond_energies=bond_energies)

        self.assertEqual(0, KineticMolecule('[O].[O]').get_potential_energy(chem))
        self.assertEqual(chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('O'), to_bond_type=2), KineticMolecule('[O]=[O]').get_potential_energy(chem))
        self.assertTrue(KineticMolecule('[O]=[O]').get_potential_energy(chem) < KineticMolecule('[O].[O]').get_potential_energy(chem))  # bonds have NEGATIVE energy, so adding bonds REDUCES PE
        O1H_bond = chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('H'), to_bond_type=1)
        O2O_bond = chem.get_bond_energy(Chem.Atom('O'), Chem.Atom('O'), to_bond_type=2)
        self.assertEqual(O2O_bond, KineticMolecule('[O]=[O]').get_potential_energy(chem))
        self.assertEqual(2 * O1H_bond, KineticMolecule('[H]O[H]').get_potential_energy(chem))
        self.assertEqual(O1H_bond, KineticMolecule('[H+].[OH-]').get_potential_energy(chem))

    def testGetStronglyConnectedComponents(self):
        mol = KineticMolecule("O=C=O.O")
        self.assertEqual(6, mol.GetNumAtoms())
        self.assertListEqual([set([0, 1, 2]), set([3, 4, 5])], mol._get_strongly_connected_components())
        mol = KineticMolecule("[H]")
        self.assertEqual([[0]], mol._get_strongly_connected_components())

    def testSameComponents(self):
        mol = KineticMolecule("O=C=O.[H]O[H]")
        self.assertFalse(mol.same_component(0, 5))
        self.assertTrue(mol.same_component(0, 1))
        mol = KineticMolecule("CO.N", components=[set(range(6)), set(range(6, 10))])
        self.assertTrue(mol.same_component(0, 2))
        self.assertFalse(mol.same_component(0, 8))

    def testAssignFormalCharge(self):
        # [H]O[H] -> [OH-]+[H+] Hydroxyl plus proton
        mol = KineticMolecule('[H].[OH]')
        mol._assign_formal_charge()
        self.assertEqual('[H+].[H][O-]', Chem.MolToSmiles(mol))
