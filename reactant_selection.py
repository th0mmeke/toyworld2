from abc import ABCMeta, abstractmethod


class ReactantSelection(object):
    __metaclass__ = ABCMeta

    def __init__(self, population):

        if not isinstance(population, list):
            raise TypeError
        self.population = population

    @abstractmethod
    def get_reactants(self):

        """
        Idempotent method to identify the next reactants from the reactor population.

        :return: [IElement]
        """

        pass

    def react(self, reaction):
        '''
        Action the provided Reaction by replacing the reactants in the population with the products.
        Raises ValueError if any of the IElements provided as reactants doesn't match an item in the population.

        :param reaction: Reaction
        :return: StateRecord

        '''

        for x in reaction.get_reactants():
            self.population.remove(x)  # raises ValueError if this reactant is not in the population
        self.population.extend(reaction.get_products())

        return reaction.as_dict()

    def __str__(self):
        return ",".join([str(x) for x in self.population])
