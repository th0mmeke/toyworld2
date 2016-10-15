import random
from reactant_selection import ReactantSelection


class UniformReactantSelection(ReactantSelection):

    def get_reactants(self):

        """
        Idempotent method to identify the next reactants from the reactor population.

        :return: [IElement]
        """

        if len(self.population) > 1:
            reactants = random.sample(self.population, 2)
        else:
            reactants = self.population
        assert type(reactants) == list
        return reactants
