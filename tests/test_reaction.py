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
        r = Reaction(reactants=[element], products=[element])
        self.assertIsInstance(r, Reaction)

    def test_get_attributes(self):
        products = [1,2]
        reactants = [3]
        r = Reaction(reactants=reactants, products=products)
        self.assertEqual(r.products, products)
        self.assertEqual(r.reactants, reactants)

    def test_as_dict(self):
        products = [1, 2]
        reactants = [3]
        r = Reaction(reactants=reactants, products=products)
        actual = r.as_dict()
        self.assertEquals(actual['reactants'], [3])
        self.assertEquals(actual['products'], [1, 2])
