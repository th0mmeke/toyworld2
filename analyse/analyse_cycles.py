import os
import json
import random
import glob
import time

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(data_filepath, evaluator_filepath):

    p = 0.2

    with open(data_filepath) as f:
        try:
            state = json.load(f)
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

            with open(evaluator_filepath, mode='w') as f2:
                json.dump(population_stoichiometry, f2)


datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
evaldir = '/home/cosc/guest/tjy17/Dropbox/Experiments/Sampled-at-0.20'
# datadir = 'C:\Users\Thom\Dropbox/Experiments' /home/cosc/guest/tjy17/Dropbox/Experiments/1489554358-0-3-0-cycles.json

files = sorted(glob.glob(os.path.join(datadir, '*.json')))
for data_filepath in files:
    evaluator_filename = os.path.splitext(os.path.basename(data_filepath))[0] + "-cycles.json"
    evaluator_filepath = os.path.join(evaldir, evaluator_filename)
    if not os.path.exists(evaluator_filepath):
        print(data_filepath, evaluator_filepath)
        evaluate(data_filepath, evaluator_filepath)
