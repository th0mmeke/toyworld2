from i_element import IElement


class Molecule(IElement):

    def __init__(self, symbol):
        """
        Initialise a new Molecule represented by character string
        :param symbol: string
        """

        if not isinstance(symbol, basestring):
            raise TypeError

        self.symbol = symbol

    def get_symbol(self):
        return self.symbol

    @classmethod
    def join(cls, elements):
        pass

    @classmethod
    def split(cls, element):
        pass

    def __cmp__(self, other):
        return cmp(self.symbol, other.symbol)
