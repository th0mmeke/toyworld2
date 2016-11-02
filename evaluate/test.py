import json
import numpy as np


class Evaluate(object):

    def __init__(self, filename="../data/toyworld2.json"):

        with open(filename) as f:
            reactions = json.load(f)

        elements = set(reactions['initial_population'])  # COPY initial items
        for reaction in reactions['reactions']:
            elements.update(reaction['products'])
        number_of_unique_elements = len(elements)

        # Population is an array (number_of_generations + 1) x number_of_unique_elements_observed_in_the_population
        self.population = np.zeros((len(reactions['reactions'])+1, number_of_unique_elements), dtype=np.dtype(int))

        self.unique_species = list(elements)
        for x in reactions['initial_population']:
            self.population[0, self.unique_species.index(x)] += 1

        i = 1
        for reaction in reactions['reactions']:
            self.population[i, :] = self.population[i-1, :]
            for x in reaction['reactants']:
                self.population[i, self.unique_species.index(x)] -= 1
            for y in reaction['products']:
                self.population[i, self.unique_species.index(y)] += 1
            i += 1

        self.final = np.zeros(number_of_unique_elements, dtype=np.dtype(int))
        for x in reactions['final_population']:
            self.final[self.unique_species.index(x)] += 1

        assert list(self.final) == list(self.population[len(reactions['reactions']), :])


e = Evaluate()

delta_c = e.population[:, np.array(map(lambda x: len(x) > 10, e.unique_species), dtype=bool)]

presence = np.apply_along_axis(lambda x: x > 0, 0, delta_c)
delta_p = np.apply_along_axis(sum, 0, presence)
print(delta_p)
max = np.apply_along_axis(max, 0, longer_elements)
print(max)


