import unittest

from atom import Atom


class TestAtom(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(ValueError):
            Atom(symbol='')
        with self.assertRaises(ValueError):
            Atom(symbol="ab")
        self.assertIsInstance(Atom(symbol="A"), Atom)

    def test_get_symbol(self):
        self.assertEqual(Atom(symbol="A").get_symbol(), "A")
