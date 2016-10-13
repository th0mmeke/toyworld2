import unittest

from reaction_chemistry import ReactionChemistry

from atom import Atom
from molecule import Molecule
from tests.stub_chemistry import StubChemistry


class TestReactionChemistry(unittest.TestCase):

    @unittest.skip("While designing Chemistry")
    def test_base_get_products(self):

        """
        Given the use of the StubChemistry
        When get_products(<Reactants>) is called
        Then the response is given by <Products>

        | Reactants         | Products      | Notes
        | Atom('A')         | -             | Raises TypeError
        | []                | []            | Raises ValueError

        """

        rc = ReactionChemistry(StubChemistry, lambda x: x)
        with self.assertRaises(TypeError):
            rc.get_products(Atom('A'))
        with self.assertRaises(ValueError):
            rc.get_products([])

    @unittest.skip("While designing Chemistry")
    def test_get_products(self):

        """
        Given the use of the StubChemistry
        And using <Selection> for the value of the selection parameter
        When get_products(<Reactants>) is called
        Then expect <Products>

        | Selection         | Reactants               | Products                            |
        | lambda x: x       | [Molecule('A')]         | [Molecule('A')]                     |

        """

        rc = ReactionChemistry(StubChemistry, lambda x: x)

        self.assertListEqual(rc.get_products(reactants=[Molecule('F')]), [Molecule('A'), Molecule('B')])


