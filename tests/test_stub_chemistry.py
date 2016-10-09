import unittest

from atom import Atom
from molecule import Molecule
from stub_chemistry import StubChemistry


class TestStubChemistry(unittest.TestCase):

    @unittest.skip("Not implemented")
    def test_join(self):
        self.assertEqual(StubChemistry.join([Atom('A'), Atom('B')]), Molecule('AB'))

    @unittest.skip("Not implemented")
    def test_split(self):
        self.assertEqual(StubChemistry.split(Molecule('AB')), [Atom('A'), Atom('B')])

