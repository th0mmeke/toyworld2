import unittest
import json

from reaction import Reaction


class TestReaction(unittest.TestCase):

    def test_init(self):
        element = type('test', (), {})()
        with self.assertRaises(TypeError):
            Reaction(reactants=element, products=element)
        with self.assertRaises(TypeError):
            Reaction(reactants=element, products=[element])
        with self.assertRaises(TypeError):
            Reaction(reactants=[element], products=element)
        with self.assertRaises(TypeError):
            Reaction(reactants='string', products='string')
        with self.assertRaises(ValueError):
            Reaction(reactants=[], products=[element])
        r = Reaction(reactants=[element], products=[element])
        self.assertIsInstance(r, Reaction)

    def test_get_molecules(self):
        reactants = [3]
        r = Reaction(reactants=reactants)
        self.assertEqual(r.get_reactants(), reactants)
        products = [1, 2]
        r = Reaction(reactants=reactants, products=products)
        self.assertEqual(r.get_products(), products)
        self.assertEqual(r.get_reactants(), reactants)

    def test_values(self):
        self.assertEqual(Reaction(reactants=[1], products=[1]).get_reactant_value(), 0)
        self.assertEqual(Reaction(reactants=[1], products=[1], reactant_value=-1).get_reactant_value(), -1)
        self.assertEqual(Reaction(reactants=[1], products=[1], reactant_value=10).get_reactant_value(), 10)
        self.assertEqual(Reaction(reactants=[1], products=[1]).get_product_value(), 0)
        self.assertEqual(Reaction(reactants=[1], products=[1], product_value=-1).get_product_value(), -1)
        self.assertEqual(Reaction(reactants=[1], products=[1], product_value=11).get_product_value(), 11)

    def test_as_dict(self):
        products = [1, 2]
        reactants = [3]
        r = Reaction(reactants=reactants, products=products)
        actual = r.as_dict()
        self.assertEquals(actual['reactants'], [3])
        self.assertEquals(actual['products'], [1, 2])
