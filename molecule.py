from i_element import IElement


class Molecule(IElement):

    """
    Molecules are conceptually just combinations of Atoms, described by some character string symbol.
    The symbol has no particular meaning except to an IChemistry. For example, an IChemistry based on RDKit would expect symbols in SMILES form.

    Any transformations of Molecules is the responsibility of a IChemistry.
    """

    def __init__(self, symbol):
        """
        Initialise a new Molecule represented by character string.

        :param symbol: string
        """

        if not isinstance(symbol, basestring):
            raise TypeError

        self.symbol = symbol

    def get_symbol(self):
        return self.symbol

    def __str__(self):
        return self.symbol
