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

# Droop2012 Delta_C activity increment - the count of a component at a time-step
delta_c = e.population

# Bedau1998 Delta (Delta_P in Droop) activity increment - presence of a component at a time step
delta_p = np.apply_along_axis(lambda x: x > 0, 0, e.population)

evolutionary_activity_p = np.cumsum(delta_p, axis=0, dtype=float)
evolutionary_activity_c = np.cumsum(delta_c, axis=0, dtype=float)

# Bedau1998 diversity
diversity = np.apply_along_axis(sum, 1, delta_p)  # diversity for each t

# divide sum of each row by the equivalent row of diversity
# Result is A_cum by time
mean_cumulative_evolutionary_activity_c = np.apply_along_axis(sum, 1, evolutionary_activity_c) / diversity
mean_cumulative_evolutionary_activity_p = np.apply_along_axis(sum, 1, evolutionary_activity_p) / diversity


#delta_c = e.population[:, np.array(map(lambda x: len(x) > 10, e.unique_species), dtype=bool)]


