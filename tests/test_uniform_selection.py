import unittest

from uniform_selection import UniformSelection
from atom import Atom


class TestUniformSelection(unittest.TestCase):

    def test_get_products(self):

        '''
        When get_products(<Elements>) is called
        Then the response is given by <Products>

        | Elements                              | Products                  | Notes
        | Atom('A')                             | -                         | Raises ValueError
        | []                                    | []                        |
        | [Atom('A')]                           | [Atom('A')]               | Not necessarily same Atom('A') object
        | [Atom('A'),Atom('B')]                 | [Atom('A'),Atom('B')]     | Not necessarily same objects
        | [Atom('A'),Atom('B'),Atom('C')]       | Any combination of 2      | Not necessarily same objects

        '''

        with self.assertRaises(ValueError):
            UniformSelection.get_products(Atom('A'))

        self.assertEqual(UniformSelection.get_products(reactants=[]), [])
        self.assertListEqual(UniformSelection.get_products(reactants=[Atom('A')]), [Atom('A')])

        actual = UniformSelection.get_products(reactants=[Atom('A'), Atom('B')])
        expected = [Atom('A'), Atom('B')]
        self.assertEqual(len(actual), len(expected))
        for a, e in zip(sorted(actual), sorted(expected)):
            self.assertEqual(a, e)

        # Following isn't a correct test as it compares object identities rather than symbol values:
        # self.assertItemsEqual(UniformSelection.get_products(reactants=[Atom('A'), Atom('B')]), [Atom('A'), Atom('B')])

        actual = UniformSelection.get_products(reactants=[Atom('A'), Atom('B'), Atom('C')])
        self.assertEqual(len(actual), 2)

