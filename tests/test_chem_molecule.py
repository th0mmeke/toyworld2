import unittest

from chem_molecule import ChemMolecule
from rdkit.rdBase import DisableLog, EnableLog


class TestMolecule(unittest.TestCase):

    def setUp(self):
        DisableLog('rdApp.error')  # Doesn't seem to disable SMILES Parse Error though...

    def tearDown(self):
        EnableLog('rdApp.error')

    def test_init(self):
        with self.assertRaises(TypeError):
            ChemMolecule()
        with self.assertRaises(TypeError):
            ChemMolecule(symbol=12)
        with self.assertRaises(ValueError):
            ChemMolecule(symbol="AB")
        self.assertIsInstance(ChemMolecule(symbol="C"), ChemMolecule)
        self.assertIsInstance(ChemMolecule(symbol="Br"), ChemMolecule)

    def test_add_Hs(self):
        self.assertEqual(ChemMolecule(symbol="O").get_symbol(), "[H]O[H]")
