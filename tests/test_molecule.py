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

    def test_eq(self):
        self.assertEqual(Molecule('A'), Molecule('A'))
        self.assertNotEqual(Molecule('A'), Molecule('B'))
        self.assertIn(Molecule('A'), [Molecule('A'), Molecule('C')])
        self.assertListEqual([Molecule('A'), Molecule('C')], [Molecule('A'), Molecule('C')])
