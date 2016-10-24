import random
import logging
from i_reactant_selection import IReactantSelection
from reaction import Reaction


class LengthBiasedReactantSelection(IReactantSelection):

    def __init__(self, population, **kwargs):

        if not isinstance(population, list):
            raise TypeError
        self.population = sorted(population, key=lambda x: len(x.get_symbol()), reverse=True)

    def get_reactants(self):

        if len(self.population) > 1:
            reactants = random.sample(self.population[:len(self.population)/4], 2)
        else:
            reactants = self.population
        assert type(reactants) == list
        return Reaction(reactants=reactants)

    def react(self, reaction):

        for x in reaction.get_reactants():
            self.population.remove(x)  # raises ValueError if this reactant is not in the population
        self.population.extend(reaction.get_products())
        self.population = sorted(self.population, key=lambda x: len(x.get_symbol()), reverse=True)
        return reaction.as_dict()

    def get_population(self):
        return self.population

    def __str__(self):
        return ",".join([str(x) for x in self.population])
