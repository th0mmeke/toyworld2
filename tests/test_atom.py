import unittest

from atom import Atom


class TestAtom(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(TypeError):
            Atom()
        with self.assertRaises(ValueError):
            Atom(symbol='')
        with self.assertRaises(ValueError):
            Atom(symbol="ab")
        self.assertIsInstance(Atom(symbol="A"), Atom)

    def test_get_symbol(self):
        self.assertEqual(Atom(symbol="A").get_symbol(), "A")

    def test_eq(self):
        self.assertEqual(Atom('A'), Atom('A'))
        self.assertNotEqual(Atom('A'), Atom('B'))
        self.assertIn(Atom('A'), [Atom('A'), Atom('C')])
        self.assertListEqual([Atom('A'), Atom('C')], [Atom('A'), Atom('C')])
