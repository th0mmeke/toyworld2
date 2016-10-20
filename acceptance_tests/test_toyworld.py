import unittest
from collections import Counter

from toyworld2 import ToyWorld
from chem_molecule import ChemMolecule
from uniform_reactant_selection import UniformReactantSelection
import uniform_product_selection
from semi_realistic_chemistry import SemiRealisticChemistry
from state import State
import bond_energies

class TestToyWorld(unittest.TestCase):

    def test_basic_scenario(self):

        """
        Feature: Iteratively apply reactions to generate a final population
          # Enter feature description here

          Scenario: Reactions change population under basic scenario
            Given UniformReactantSelection
            And UniformProductSelection
            And population = ['O','C']
            When ToyWorld runs for 5 generations
            Then population does not equal ['O','C']
        """

        chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)

        defn = {'O': 1, 'C': 1}
        population = []
        for symbol, quantity in defn.iteritems():
            for i in range(quantity):
                population.append(ChemMolecule(symbol))

        initial_population = Counter([str(x) for x in population])
        tw = ToyWorld(reactor=UniformReactantSelection(population=population, ke=100),
                      chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies),
                      product_selection=uniform_product_selection.product_selection)
        tw.run(generations=5, state=State(persistence=lambda x: x))

        self.assertNotEqual(initial_population, Counter([str(x) for x in tw.reactor.get_population()]))


