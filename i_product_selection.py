class IProductSelection(object):

    @classmethod
    def get_products(cls, reactants, chemistry):

        """
        Return the products of a reaction for a given set of reactants, assuming some particular Chemistry.

        Reactants can be Atoms, Molecules or subclasses of either.

        :param reactants: list of IElement
        :param chemistry: IChemistry
        :return: list of IElement
        """

        pass
