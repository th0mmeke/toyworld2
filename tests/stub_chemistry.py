from i_chemistry import IChemistry
from molecule import Molecule

class StubChemistry(IChemistry):

    @classmethod
    def join(cls, elements):
        if not isinstance(elements, list):
            raise ValueError
        return [Molecule('C')]

    @classmethod
    def split(cls, element):
        return [Molecule('A'), Molecule('B')]
