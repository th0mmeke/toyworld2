import unittest

from molecule import Molecule


class TestMolecule(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(TypeError):
            Molecule()
        with self.assertRaises(TypeError):
            Molecule(symbol=12)
        self.assertIsInstance(Molecule(symbol="AB"), Molecule)

    def test_get_symbol(self):
        self.assertEqual(Molecule(symbol="Ab").get_symbol(), "Ab")

    def test_id(self):
        # Check unique ids
        self.assertNotEqual(Molecule('O').get_id(), Molecule('O').get_id())
        a = Molecule('O').get_id()
        b = Molecule('O').get_id()
        self.assertNotEqual(a, b)
