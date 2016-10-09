import random


class UniformSelection(object):

    @classmethod
    def get_products(cls, reactants):

        """
        Return a list of products given a list of reactants.

        Reactants can be Atoms, Molecules or subclasses of either. All Molecules must share the same Chemistry.

        :param reactants: list of IElement
        :return: list of IElement
        """

        if not isinstance(reactants, list):
            raise ValueError

        if len(reactants) < 2:
            return reactants

        return random.sample(reactants, 2)
