import random
from i_reactant_selection import IReactantSelection
from reaction import Reaction


class UniformReactantSelection(IReactantSelection):

    def __init__(self, population, **kwargs):

        if not isinstance(population, list):
            raise TypeError
        self.population = population

    def update_environment(self, target, value):
        pass

    def get_reactants(self):

        if len(self.population) > 1:
            reactants = random.sample(self.population, 2)
        else:
            reactants = self.population
        assert type(reactants) == list
        return Reaction(reactants=reactants)

    def react(self, reaction):

        for x in reaction.get_reactants():
            self.population.remove(x)  # raises ValueError if this reactant is not in the population
        self.population.extend(reaction.get_products())

        return reaction.as_dict()

    def get_population(self):
        return self.population

    def __str__(self):
        return ",".join([str(x) for x in self.population])
