import os
import json
import random
import glob
import time

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(data_filepath, evaluator_filepath):

    with open(data_filepath) as f:
        state = json.load(f)
    try:
        e = IdentifyMoleculeCycles(reactions=state['reactions'])
    except:
        pass
    else:
        print(len(e.reactants), len(e.reactant_molecules))

        population_stoichiometry = []
        sample = random.sample(e.reactant_molecules, len(e.reactant_molecules)/20)
        for mol in sample:
            for item in e.sample_reactant_stoichiometry(mol, minimum_stoichiometry=2, max_depth=30):
                population_stoichiometry.append(item)

        with open(evaluator_filepath, mode='w') as f:
            json.dump(population_stoichiometry, f)


datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
datadir = 'C:\Users\Thom\Dropbox/Experiments'

while True:
    for data_filepath in sorted(glob.glob(os.path.join(datadir, '1489311243*.json'))):
        evaluator_filename = os.path.splitext(os.path.basename(data_filepath))[0] + "-cycles.json"
        evaluator_filepath = os.path.join(datadir, evaluator_filename)
        if not os.path.exists(evaluator_filepath):
            print(data_filepath, evaluator_filepath)
            evaluate(data_filepath, evaluator_filepath)
    time.sleep(60)
