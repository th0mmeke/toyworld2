import unittest

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
        Given a UniformReactor that contains the only <Elements> elements
        When get_reactants is called
        Then the response equals <Reactants>

        | Elements  | Reactants |
        | []        | []        |
        | ['A']     | []        |
        | ['A','B'] | ['A','B'  |

        :return:
        '''
        self.assertEqual(UniformReactor(population=[]).get_reactants(), [])
        self.assertEqual(UniformReactor(population=[]).get_reactants(), [])
        self.assertItemsEqual(UniformReactor(population=['A', 'B']).get_reactants(), ['B', 'A'])  # order isn't important

    @unittest.skip("Still to be done...")
    def test_react(self):
        pass
