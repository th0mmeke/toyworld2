from reactant_selection import ReactantSelection

import random

class UniformReactor(ReactantSelection):

    def __init__(self, population):

        if not isinstance(population, list):
            raise ValueError
        self.population = population

    def get_reactants(self):
        '''

        :return: list of Element
        '''
        if len(self.population) > 1:
            return random.sample(self.population, 2)
        else:
            return []

    def react(self, reaction):
        '''

        :return: StateRecord
        '''