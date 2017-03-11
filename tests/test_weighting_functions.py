import random
import unittest

import weighting_functions
from reaction import Reaction
from molecule import Molecule


class TestWeightingFunctions(unittest.TestCase):
    def setUp(self):
        self.reactions = [Reaction(reactants=[Molecule('ABC'), Molecule('DE')],
                                   products=[Molecule('ABC'), Molecule('XY')],
                                   reactant_value=0, product_value=10),
                          Reaction(reactants=[Molecule('ABC'), Molecule('FG')],
                                   products=[Molecule('AB'), Molecule('ABC')],
                                   reactant_value=30, product_value=20),
                          Reaction(reactants=[Molecule('ABCDEF'), Molecule('FG')],
                                   products=[Molecule('AB'), Molecule('BCDE'), Molecule('F')],
                                   reactant_value=0, product_value=-20)]

    def testUniformWeighting(self):
        self.assertListEqual([1, 1, 1], [weighting_functions.uniform_weighting(r, 0) for r in self.reactions])

    def testLeastEnergyWeighting(self):
        # -product_value, if product_value < 0
        # 0 if reactant_value < product_value
        # reactant_value + product_value otherwise
        self.assertListEqual([0, 10, 20], [weighting_functions.least_energy_weighting(r, 0) for r in self.reactions])
