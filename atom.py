class Atom(object):

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
