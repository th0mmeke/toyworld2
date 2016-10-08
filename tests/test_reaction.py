
import unittest

from reaction import Reaction


class Test(unittest.TestCase):

    def test_init(self):
        element = type('test', (), {})()
        with self.assertRaises(ValueError):
            Reaction(reactants=element, products=element)
        with self.assertRaises(ValueError):
            Reaction(reactants=element, products=[element])
        with self.assertRaises(ValueError):
            Reaction(reactants=[element], products=element)
        with self.assertRaises(ValueError):
            Reaction(reactants='string', products='string')
        r = Reaction(reactants=[element], products=[element])
        self.assertIsInstance(r, Reaction)

    def test_get_attributes(self):
        products = [1,2]
        reactants = [3]
        r = Reaction(reactants=reactants, products=products)
        self.assertEqual(r.products, products)
        self.assertEqual(r.reactants, reactants)
