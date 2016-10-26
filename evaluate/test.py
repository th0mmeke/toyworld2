import json
import numpy as np
import unittest

# state = State(filename="toyworld2.json")
# for i in range(5000000):
#     state.add(Reaction([Molecule('A'), Molecule('B')], products=[Molecule('C'), Molecule('D')], reactant_value=i, product_value=i).as_dict())
# del state

with open('../data/toyworld2.json') as f:
    reactions = json.load(f)
    print(reactions['initial_population'])

    # 1. Long molecules
    # 2. That are similar to each other
    # 3. And that in aggregate show a period of exponential growth (assumes A+X->2A+Y rather than A+X->A+Y...)
    # 3a. n_t+1 approx 2n_t

    elements = set(reactions['initial_population'])  # COPY initial items
    for reaction in reactions['reactions']:
        elements.update(reaction['products'])
    number_of_unique_elements = len(elements)
    population = np.zeros((len(reactions['reactions'])+1, number_of_unique_elements), dtype=np.dtype(int))

    print(len(reactions['reactions']))
    unique_species = list(elements)
    for x in reactions['initial_population']:
        population[0, unique_species.index(x)] += 1

    i = 1
    for reaction in reactions['reactions']:
        population[i, :] = population[i-1, :]
        for x in reaction['reactants']:
            population[i, unique_species.index(x)] -= 1
        for y in reaction['products']:
            population[i, unique_species.index(y)] += 1
        i += 1

    final = np.zeros(number_of_unique_elements, dtype=np.dtype(int))
    for x in reactions['final_population']:
        final[unique_species.index(x)] += 1

    assert list(final) == list(population[len(reactions['reactions']), :])

    print([x for x in unique_species if len(x) > 10])
    print(population[:, np.array(map(lambda x: len(x) > 10, unique_species), dtype=bool)])

