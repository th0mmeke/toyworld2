import unittest

from reaction import Reaction
from state_record import StateRecord
from reactor import Reactor


class TestReactor(unittest.TestCase):

    def test_init(self):
        element = type('element', (), {})()
        with self.assertRaises(ValueError):
            Reactor(population=element, reactant_selection=lambda x: x)
        r = Reactor(population=[element], reactant_selection=lambda x: x)
        self.assertIsInstance(r, Reactor)

    def test_get_reactants(self):
        """
        Test idempotence of get_reactants()
        """

        p = ['A', 'B', 'C']
        r = Reactor(population=p, reactant_selection=lambda x: x)
        self.assertEqual(r.get_reactants(), p)
        self.assertEqual(r.get_reactants(), p)

    def test_get_reactants_with_reactant_selection(self):
        """
        Given a Reactor that contains only <Elements> elements
        And there are at least 2 elements
        And the reactant_selection method is to return the first two elements
        When get_reactants is called
        Then the response contains the first two elements of the population
        """

        r = Reactor(population=['A', 'B', 'C', 'D'], reactant_selection=lambda x: x[:2])
        self.assertEqual(r.get_reactants(), ['A', 'B'])

    def test_react(self):
        """
        Given a Reactor that contains <Population> elements
        When react is called with a Reaction where reactants=<Reactants> and products=<Products>
        Then the Reactor contains <Expected> elements

        | Population  | Reactants | Products  | Expected  |
        | []        | []        | []        | []        |
        | ['A']     | ['A']     | ['B','C'] | ['B,'C']  |

        Rather than checking private variables we look for the implications of the change through existing methods.
        """

        r = Reactor(population=[], reactant_selection=lambda x: x )
        sr = r.react(Reaction(reactants=[], products=[]))
        self.assertItemsEqual(r.get_reactants(), [])
        self.assertIsInstance(sr, StateRecord)

        r = Reactor(population=['A'], reactant_selection=lambda x: x)
        sr = r.react(Reaction(reactants=['A'], products=['B', 'C']))
        self.assertItemsEqual(r.get_reactants(), ['B', 'C'])
        self.assertIsInstance(sr, StateRecord)

    def test_react_errors(self):
        """
        Throw a ValueError and leaves the population unchanged:

        | Population    | Reactants | Products  | Expected  |
        | ['A']         | ['D']     | -         | -         |
        | ['A']         | ['A','D'] | -         | -         |
        """

        r = Reactor(population=['A'], reactant_selection=lambda x: x)
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=['D'], products=['B', 'C']))

        r = Reactor(population=['A'], reactant_selection=lambda x: x)
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=['A', 'D'], products=['B', 'C']))





