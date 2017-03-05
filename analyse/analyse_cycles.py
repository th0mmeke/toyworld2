import os
import json
import random
import glob

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(filename, datadir):
    basename, ext = os.path.splitext(filename)

    evaluator = 'molecules'
    evaluator_filename = os.path.join(datadir, '{}-{}sample.json'.format(basename, evaluator))
    data_filename = os.path.join(datadir, filename)
    counts = []

    with open(data_filename) as f:
        state = json.load(f)
    try:
        e = IdentifyMoleculeCycles(reactions=state['reactions'])
    except:
        pass
    else:
        print(len(e.reactants), len(e.reactant_molecules))

        for count in range(0, 10):
            population_stoichiometry = []
            sample = random.sample(e.reactant_molecules, len(e.reactant_molecules)/10)
            for mol in sample:
                for item in e.sample_reactant_stoichiometry(mol):
                    population_stoichiometry.append(item)
            counts.append(str(len(population_stoichiometry)))

        with open(evaluator_filename, mode='w') as f:
            json.dump(population_stoichiometry, f)

    return counts

datadir = '/Users/Thom/Workspace/toyworld2'
with open('sample_count.csv', mode='w') as f_counts:
    for filename in glob.glob(os.path.join(datadir, '1488600998*.json')):
        print(filename)
        counts = evaluate(filename, datadir)
        print(counts)
        f_counts.write(",".join(counts) + "\n")
