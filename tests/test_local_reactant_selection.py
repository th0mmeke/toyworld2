import unittest

from reaction import Reaction
from local_reactant_selection import LocalReactantSelection
from molecule import Molecule


class TestLocalReactantSelection(unittest.TestCase):

    def test_init(self):
        element = type('element', (), {})()
        with self.assertRaises(TypeError):
            LocalReactantSelection(population=element)
        r = LocalReactantSelection(population=[element])
        self.assertIsInstance(r, LocalReactantSelection)

    def test_get_reactants(self):
        self.assertIsInstance(LocalReactantSelection(population=['A']).get_reactants(), Reaction)
        self.assertEqual(LocalReactantSelection(population=['A']).get_reactants().get_reactants(), ['A'])

    def test_react(self):
        """
        Given a SimpleReactor that contains <Population> elements
        When react is called with a Reaction where reactants=<Reactants> and products=<Products>
        Then the SimpleReactor contains <Expected> elements

        | Population  | Reactants | Products  | Expected  |
        | []        | []        | []        | []        |
        | ['A']     | ['A']     | ['B','C'] | ['B,'C']  |

        """

        m = Molecule('M')
        r = LocalReactantSelection(population=[m])
        sr = r.react(Reaction(reactants=[m], products=[]))
        self.assertItemsEqual(r.population, [m])

        x = Molecule('X')
        y = Molecule('Y')
        r = LocalReactantSelection(population=[m])
        sr = r.react(Reaction(reactants=[m], products=[x, y]))
        self.assertItemsEqual(r.population, [m, x, y])

    def test_react_errors(self):
        """
        Throw a ValueError and leaves the population unchanged:

        | Population    | Reactants | Products  | Expected  |
        | ['A']         | ['D']     | -         | -         |
        | ['A']         | ['A','D'] | -         | -         |
        """

        r = LocalReactantSelection(population=['A'])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=['D'], products=['B', 'C']))

        r = LocalReactantSelection(population=['A'])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=['A', 'D'], products=['B', 'C']))





