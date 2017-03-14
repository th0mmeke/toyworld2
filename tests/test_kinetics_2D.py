import random
import unittest
import math

import pymunk as pm

from kinetics_2D import Kinetics2D


class TestKinetics2D(unittest.TestCase):

    def setUp(self):
        # logging.basicConfig(level=logging.DEBUG)
        self._kinetics = Kinetics2D

    def test_get_ke(self):
        self.assertAlmostEqual(12.5, self._kinetics.get_ke(1, 3, 4))

    def test_reaction_energy(self):
        mol0 = pm.Body(mass=1, moment=1)
        mol1 = pm.Body(mass=1, moment=1)
        v = pm.Vec2d(random.uniform(0, 10), 0)
        v.angle = random.uniform(-math.pi, math.pi)
        mol0.velocity = v
        mol1.velocity = -v
#       self.assertAlmostEqual(mol0.kinetic_energy + mol0.kinetic_energy, self._kinetics.get_cm_energy([mol0, mol0]))
        self.assertAlmostEqual(0, self._kinetics.get_cm_energy([mol0, mol1]))

    def test_get_cm_velocity(self):
        mol0 = pm.Body(mass=1, moment=1)
        mol1 = pm.Body(mass=1, moment=1)
        mol0.velocity = pm.Vec2d(random.uniform(0, 10), 0)
        mol0.velocity.angle = random.uniform(-math.pi, math.pi)
        mol1.velocity = -mol0.velocity
        cm_v = self._kinetics.get_cm_velocity([mol0, mol1])
        self.assertAlmostEqual(0, cm_v.length)  # velocity = 0 given rounding

    def test_get_cm_energy(self):
        mol0 = pm.Body(mass=1, moment=1)
        mol1 = pm.Body(mass=1, moment=1)
        mol0.velocity = pm.Vec2d(random.uniform(0, 10), 0)
        mol0.velocity.angle = random.uniform(-math.pi, math.pi)
        mol1.velocity = -mol0.velocity
        self.assertAlmostEqual(0, self._kinetics.get_cm_energy([mol0, mol1]))  # equal and opposite motion -> cm stationary

    def test_inelastic_collision_angles(self):

        mol0 = pm.Body(mass=10, moment=1)
        mol1 = pm.Body(mass=2, moment=1)
        mol0.velocity = pm.Vec2d(random.uniform(0, 10), 0)
        mol0.velocity.angle = random.uniform(-math.pi, math.pi)
        mol1.velocity = pm.Vec2d(random.uniform(0, 10), 0)
        mol1.velocity.angle = random.uniform(-math.pi, math.pi)

        reactants = [mol0, mol1]
        products = [pm.Body(mass=4, moment=1), pm.Body(mass=6, moment=1), pm.Body(mass=2, moment=1)]

        delta_ke, out_v = self._kinetics.inelastic_collision(reactants, [p.mass for p in products])

        for mol, v in zip(products, out_v):
            mol.velocity = v

        # test that momentum conserved
        in_CM = self._kinetics.get_cm_velocity(reactants)
        out_CM = self._kinetics.get_cm_velocity(products)
        for in_, out_ in zip(in_CM, out_CM):
            self.assertAlmostEqual(in_, out_)

        # test that the out molecules do not all have the same velocity... possible to have correct CM and energy if all along CoM, but boring!
        # max angle - min angle > delta
        angles = [i.angle for i in out_v]
        self.assertGreater(max(angles) - min(angles), math.pi/100)

    def test_inelastic_collision_double(self):

        mol0 = pm.Body(mass=36.005, moment=1)
        mol1 = pm.Body(mass=18.015, moment=1)
        mol0.velocity = pm.Vec2d(-6.58286077125, -24.2936918394)
        mol1.velocity = pm.Vec2d(-32.4722726889, -0.437424241464)

        reactants = [mol0, mol1]
        products = [pm.Body(mass=36.005, moment=1), pm.Body(mass=18.015, moment=1)]

        failed = 0
        for i in range(0, 1000):
            try:
                out_v = self._kinetics.inelastic_collision(reactants, [p.mass for p in products])
            except ValueError:
                failed += 1
            self.assertEqual(0, failed)

    def test_inelastic_collision_single(self):

        mol0 = pm.Body(mass=1.008, moment=1)
        mol1 = pm.Body(mass=15.999, moment=1)
        mol0.velocity = pm.Vec2d(-66.6650078772, 25.300112506)
        mol1.velocity = pm.Vec2d(-0.715752846158, -0.420109446346)

        reactants = [mol0, mol1]
        products = [pm.Body(mass=17.007, moment=1)]

        out_v = self._kinetics.inelastic_collision(reactants, [p.mass for p in products])



if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testInit']
    unittest.main()
