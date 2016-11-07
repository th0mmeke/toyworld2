import json
import numpy as np
import os
import argparse


def load_json(filename):
    with open(filename) as f:
        reactions = json.load(f)

    elements = set(reactions['initial_population'])  # COPY initial items
    for reaction in reactions['reactions']:
        elements.update(reaction['products'])
    number_of_unique_elements = len(elements)

    # Population is an array (number_of_generations + 1) x number_of_unique_elements_observed_in_the_population
    _population = np.zeros((len(reactions['reactions']) + 1, number_of_unique_elements), dtype=np.dtype(int))

    _unique_species = list(elements)
    for x in reactions['initial_population']:
        _population[0, _unique_species.index(x)] += 1

    i = 1
    for reaction in reactions['reactions']:
        _population[i, :] = _population[i - 1, :]
        for x in reaction['reactants']:
            _population[i, _unique_species.index(x)] -= 1
        for y in reaction['products']:
            _population[i, _unique_species.index(y)] += 1
        i += 1

    final = np.zeros(number_of_unique_elements, dtype=np.dtype(int))
    for x in reactions['final_population']:
        final[_unique_species.index(x)] += 1

    assert list(final) == list(_population[len(reactions['reactions']), :])
    
    return _population, _unique_species

if __name__ == "__main__":

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-i', '--input', default='toyworld2.json', help="Input filename. Can be a filename or a full pathname")
    parent_parser.add_argument('-o', '--output', default='toyworld2.csv', help="Output filename. Will be placed in same directory as input file")
    parser = argparse.ArgumentParser(parents=[parent_parser])

    args = parser.parse_args()

    dirname = os.path.dirname(args.input)
    ext = os.path.splitext(args.output)[1] or '.csv'
    basename = os.path.basename(os.path.splitext(args.output)[0])
    csv_species_filename = os.path.join(dirname, 'species-'+basename+ext)
    csv_population_filename = os.path.join(dirname, 'population-'+basename+ext)

    population, unique_species = load_json(args.input)
    
    np.savetxt(csv_species_filename, unique_species, delimiter=',', fmt='%s')
    np.savetxt(csv_population_filename, population, delimiter=',', fmt='%.3f', comments="", header=",".join(unique_species))
