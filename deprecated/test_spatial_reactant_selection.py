import unittest

import pymunk as pm

from chem_molecule import ChemMolecule
from deprecated.spatial_reactant_selection import SpatialReactantSelection
from reaction import Reaction


class TestSpatialReactantSelection(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(TypeError):
            SpatialReactantSelection()
        r = SpatialReactantSelection(population=[ChemMolecule('O')])
        self.assertIsInstance(r, SpatialReactantSelection)

    def test_get_reactants(self):
        self.assertIsNone(SpatialReactantSelection(population=[ChemMolecule('O')]).get_reactants(), None)

        defn = {"[H][H]": 100, "FO": 100, "O": 200, "[O-][N+](=O)[N+]([O-])=O": 100, "N(=O)[O]": 100, "O=C=O": 200}
        population = []
        for symbol, quantity in defn.iteritems():
            for i in range(quantity):
                population.append(ChemMolecule(symbol))

        r = SpatialReactantSelection(population=population).get_reactants()
        self.assertIsInstance(r, Reaction)
        self.assertGreaterEqual(2, len(r.get_reactants()))
        self.assertIsNone(r.get_products())

    def test_react(self):
        """
        Given a SpatialReactantSelection that contains <Population> elements
        When react is called with a Reaction where reactants=<Reactants> and products=<Products>
        Then the SpatialReactantSelection contains <Expected> elements

        Rather than checking private variables we look for the implications of the change through existing methods.
        """

        a = ChemMolecule('O')
        b = ChemMolecule('C')
        c = ChemMolecule('[H]C([H])([H])[H]')
        d = ChemMolecule('[H]O[H]')

        r = SpatialReactantSelection(population=[a, b])
        sr = r.react(Reaction(reactants=[a, b], products=[c, d]))
        self.assertItemsEqual(r.get_population(), [c, d])
        self.assertIsInstance(sr, dict)

    def test_react_errors(self):
        """
        Throw a ValueError and leaves the population unchanged:
        """

        a = ChemMolecule('O')
        b = ChemMolecule('C')
        c = ChemMolecule('[H]C([H])([H])[H]')
        d = ChemMolecule('[H]O[H]')

        r = SpatialReactantSelection(population=[c])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=[a], products=[c, d]))

        r = SpatialReactantSelection(population=[a])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=[a, b], products=[c, d]))

    def test_clamp_to_vessel(self):
        x = 1.1
        y = 2.2
        location = pm.Vec2d(x, y)
        bound = SpatialReactantSelection.REACTION_VESSEL_SIZE
        self.assertEqual(pm.Vec2d(x, y), SpatialReactantSelection._clamp_to_vessel(location))

        location = pm.Vec2d(-bound-1, -bound-1)
        self.assertEqual(pm.Vec2d(-bound, -bound), SpatialReactantSelection._clamp_to_vessel(location))

        location = pm.Vec2d(bound + 1, bound + 1)
        self.assertEqual(pm.Vec2d(bound, bound), SpatialReactantSelection._clamp_to_vessel(location))

