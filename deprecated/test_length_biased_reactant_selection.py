import unittest

from deprecated.length_biased_reactant_selection import LengthBiasedReactantSelection
from molecule import Molecule
from reaction import Reaction


class TestLengthBiasedReactantSelection(unittest.TestCase):

    def test_get_reactants(self):
        x = Molecule('ABCDE')
        y = Molecule('ABDEF')
        population = [Molecule('A'),
                      Molecule('A'),
                      Molecule('A'),
                      Molecule('A'),
                      Molecule('A'),
                      Molecule('A'),
                      x,
                      y]
        self.assertIsInstance(LengthBiasedReactantSelection(population=population).get_reactants(), Reaction)
        self.assertItemsEqual(LengthBiasedReactantSelection(population=population).get_reactants().get_reactants(), [x, y])
