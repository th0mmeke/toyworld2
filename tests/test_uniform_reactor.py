import unittest

from reaction import Reaction
from state_record import StateRecord
from uniform_reactor import UniformReactor


class TestUniformReactor(unittest.TestCase):

    def test_init(self):
        element = type('element', (), {})()
        with self.assertRaises(ValueError):
            UniformReactor(population=element)
        r = UniformReactor(population=[element])
        self.assertIsInstance(r, UniformReactor)

    def test_get_reactants(self):
        '''
        Given a UniformReactor that contains only <Elements> elements
        When get_reactants() is called
        Then the response equals <Reactants>

        | Elements  | Reactants |
        | []        | []        |
        | ['A']     | []        |
        | ['A','B'] | ['A','B'] |
        '''

        self.assertEqual(UniformReactor(population=[]).get_reactants(), [])
        self.assertEqual(UniformReactor(population=[]).get_reactants(), [])
        self.assertItemsEqual(UniformReactor(population=['A', 'B']).get_reactants(), ['B', 'A'])  # order isn't important

        '''
        Given a UniformReactor that contains only <Elements> elements
        And there are at least 2 elements
        When get_reactants is called
        Then the length of the response is always 2

        | Elements |
        | ['A','B']         |
        | ['A','B','C']     |
        | ['A','B','C','D'] |
        '''

        self.assertEqual(len(UniformReactor(population=['A', 'B']).get_reactants()), 2)
        self.assertEqual(len(UniformReactor(population=['A', 'B', 'C']).get_reactants()), 2)
        self.assertEqual(len(UniformReactor(population=['A', 'B', 'C', 'D']).get_reactants()), 2)


    def test_react(self):
        '''
        Given a UniformReactor that contains <Elements> elements
        When react is called with a Reaction where reactants=<Reactants> and products=<Products>
        Then the UniformReactor contains <Expected> elements

        | Elements  | Reactants | Products  | Expected  |
        | []        | []        | []        | []        |
        | ['A']     | ['A']     | ['B','C'] | ['B,'C']  |

        Rather than checking private variables we look for the implications of the change through existing methods.
        '''

        r = UniformReactor(population=[])
        sr = r.react(Reaction(reactants=[], products=[]))
        self.assertItemsEqual(r.get_reactants(), [])
        self.assertIsInstance(sr, StateRecord)

        r = UniformReactor(population=['A'])
        sr = r.react(Reaction(reactants=['A'], products=['B', 'C']))
        self.assertItemsEqual(r.get_reactants(), ['B', 'C'])
        self.assertIsInstance(sr, StateRecord)

        '''
        Throw a ValueError and leaves the population unchanged:
        | ['A']     | ['D']     | ['B','C'] | ['A']  |
        | ['A']     | ['A','D'] | ['B','C'] | ['A']  |
        '''

        r = UniformReactor(population=['A'])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=['D'], products=['B', 'C']))
        self.assertItemsEqual(r.get_reactants(), [])

        r = UniformReactor(population=['A'])
        with self.assertRaises(ValueError):
            r.react(Reaction(reactants=['A', 'D'], products=['B', 'C']))
        self.assertItemsEqual(r.get_reactants(), [])





