import os
import json
import random
import glob
import time

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(data_filepath, evaluator_filepath):

    p = 0.2

    with open(data_filepath) as f:
        state = json.load(f)
    try:
        e = IdentifyMoleculeCycles(reactions=state['reactions'])
    except:
        pass
    else:
        print(len(e.reactants), len(e.reactant_molecules))

        population_stoichiometry = []
        sample = random.sample(e.reactant_molecules, int(len(e.reactant_molecules)*p))
        for mol in sample:
            for item in e.sample_reactant_stoichiometry(mol, minimum_stoichiometry=2, max_depth=30):
                population_stoichiometry.append(item)

        with open(evaluator_filepath, mode='w') as f:
            json.dump(population_stoichiometry, f)


datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
#datadir = 'C:\Users\Thom\Dropbox/Experiments'

files = sorted(glob.glob(os.path.join(datadir, '1489554358*.json')))
for data_filepath in random.sample(files, 3):
    evaluator_filename = os.path.splitext(os.path.basename(data_filepath))[0] + "-cycles.json"
    evaluator_filepath = os.path.join(datadir, evaluator_filename)
    if not os.path.exists(evaluator_filepath):
        print(data_filepath, evaluator_filepath)
        evaluate(data_filepath, evaluator_filepath)
