import random
from i_reactant_selection import IReactantSelection


class UniformReactantSelection(IReactantSelection):

    def __init__(self, population):

        if not isinstance(population, list):
            raise TypeError
        self.population = population

    def get_reactants(self):

        if len(self.population) > 1:
            reactants = random.sample(self.population, 2)
        else:
            reactants = self.population
        assert type(reactants) == list
        return reactants

    def react(self, reaction):

        for x in reaction.get_reactants():
            self.population.remove(x)  # raises ValueError if this reactant is not in the population
        self.population.extend(reaction.get_products())

        return reaction.as_dict()

    def __str__(self):
        return ",".join([str(x) for x in self.population])
