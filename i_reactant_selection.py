from abc import ABCMeta, abstractmethod


class IReactantSelection(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_reactants(self):

        """
        Idempotent method to identify the next reactants from the reactor population.

        :return: Reaction, for the reactants only, or None if no reactants can be found.
        """

        pass

    @abstractmethod
    def react(self, reaction):

        """
        Action the provided Reaction by replacing the reactants in the population with the products.
        Raises ValueError if any of the IElements provided as reactants doesn't match an item in the population.

        :param reaction: Reaction
        :return: StateRecord

        """

        pass

    @abstractmethod
    def get_population(self):

        """
        Return a list of Molecules in the population.
        :return: [Molecule]
        """

        pass