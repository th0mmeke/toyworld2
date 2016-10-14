import random
from reactant_selection import ReactantSelection


class UniformReactantSelection(ReactantSelection):

    def get_reactants(self):

        """
        Idempotent method to identify the next reactants from the reactor population.

        :return: [IElement]
        """

        if len(self.population) > 1:
            return random.sample(self.population, 2)
        else:
            return self.population
