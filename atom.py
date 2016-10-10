from i_element import IElement

class Atom(IElement):

    """
    An Atom is the smallest element in ToyWorld2.
    """

    def __init__(self, symbol):
        """
        Initialise a new Atom represented by single character symbol
        :param symbol: Single character symbol
        """

        if not isinstance(symbol, basestring) or len(symbol) != 1:
            raise ValueError

        self.symbol = symbol

    def get_symbol(self):
        return self.symbol

    def __cmp__(self, other):
        return cmp(self.symbol, other.symbol)

    def __str__(self):
        return self.symbol