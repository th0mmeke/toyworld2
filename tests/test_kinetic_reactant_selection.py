import unittest

import pymunk as pm

from chem_molecule import ChemMolecule
from kinetic_reactant_selection import KineticReactantSelection
from reaction import Reaction


class TestKineticReactantSelection(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(TypeError):
            KineticReactantSelection()
        r = KineticReactantSelection(population=[ChemMolecule('O')])
        self.assertIsInstance(r, KineticReactantSelection)
        self.assertEqual(100, r.ke)

        r = KineticReactantSelection(population=[ChemMolecule('O')], ke=200)
        self.assertEqual(200, r.ke)
        pre_ke = sum([r.mol2body[mol].kinetic_energy for mol in r.get_population()]) / len(r.get_population())
        self.assertAlmostEqual(200, pre_ke)

    def test_adjust_ke(self):
        a = ChemMolecule('O')
        b = ChemMolecule('C')
        c = ChemMolecule('[H]C([H])([H])[H]')
        d = ChemMolecule('[H]O[H]')

        r = KineticReactantSelection(population=[a, b], ke=200)
        self.assertEqual(200, r.ke)

        r.adjust_ke(10)

        post_ke = sum([r.mol2body[mol].kinetic_energy for mol in r.get_population()])/len(r.get_population())
        self.assertAlmostEqual(210, post_ke)

    def test_react(self):
        """
        Given a KineticReactantSelection that contains <Population> elements
        When react is called with a Reaction where reactants=<Reactants> and products=<Products>
        Then the KineticReactantSelection contains <Expected> elements

        Rather than checking private variables we look for the implications of the change through existing methods.
        """

        a = ChemMolecule('O')
        b = ChemMolecule('C')
        c = ChemMolecule('[H]C([H])([H])[H]')
        d = ChemMolecule('[H]O[H]')

        r = KineticReactantSelection(population=[a, b])
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

        r = KineticReactantSelection(population=[c])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=[a], products=[c, d]))

        r = KineticReactantSelection(population=[a])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=[a, b], products=[c, d]))

    def test_clamp_to_vessel(self):
        x = 1.1
        y = 2.2
        location = pm.Vec2d(x, y)
        bound = KineticReactantSelection.REACTION_VESSEL_SIZE
        self.assertEqual(pm.Vec2d(x, y), KineticReactantSelection._clamp_to_vessel(location))

        location = pm.Vec2d(-bound-1, -bound-1)
        self.assertEqual(pm.Vec2d(-bound, -bound), KineticReactantSelection._clamp_to_vessel(location))

        location = pm.Vec2d(bound + 1, bound + 1)
        self.assertEqual(pm.Vec2d(bound, bound), KineticReactantSelection._clamp_to_vessel(location))

    def test_pymunk_ke(self):
        # PyMunk kinetic energy is mvv, not 0.5mvv, so must adjust correctly in our code
        # From the code: vsq = self.velocity.dot(self.velocity)

        import math

        ke = 100
        mol = ChemMolecule('O')
        speed = math.sqrt(2.0 * ke / mol.mass)
        velocity = pm.Vec2d(speed, 0)
        # velocity.angle = random.uniform(-math.pi, math.pi)
        self.assertAlmostEqual(ke, 0.5 * mol.mass * speed * speed)
        self.assertAlmostEqual(speed, velocity.get_length())
        inertia = 1
        body = pm.Body(mol.mass, inertia)
        body.position = pm.Vec2d(0, 0)
        body.velocity = velocity
        self.assertAlmostEqual(ke, body.kinetic_energy/2.0)

        body.velocity = pm.Vec2d(math.sqrt(ke / mol.mass), 0)
        self.assertAlmostEqual(ke, body.kinetic_energy)

    def test_ke(self):

        import math

        ke = 100
        mol = ChemMolecule('O')

        r = KineticReactantSelection([mol], ke=ke)

        r._add_molecule(mol, location=pm.Vec2d(0, 0), velocity=pm.Vec2d(math.sqrt(ke / mol.mass), 0), collision_type=KineticReactantSelection.REACTANT)
        self.assertEqual(2, len(r.get_population()))
        for body in r.mol2body.values():
            self.assertAlmostEqual(ke, body.kinetic_energy)


