from abc import ABCMeta, abstractmethod


class IChemistry(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def enumerate(self, partial_reaction):

        """
        Discover all possible options for reactions between the separate components of this Molecule.
        If partial_reaction is None, then also return None.
        If no options can be found, then return an empty list, [].

        :param partial_reaction: Reaction containing reactants and possibly reactant_value
        :return: [Reaction]
        """

        pass
