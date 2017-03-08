import os
import json
import random
import glob
import time

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(data_filepath, evaluator_filepath):

    counts = []

    with open(data_filepath) as f:
        state = json.load(f)
    try:
        e = IdentifyMoleculeCycles(reactions=state['reactions'])
    except:
        pass
    else:
        print(len(e.reactants), len(e.reactant_molecules))

        population_stoichiometry = []
        sample = random.sample(e.reactant_molecules, len(e.reactant_molecules))
        for mol in sample:
            for item in e.sample_reactant_stoichiometry(mol):
                population_stoichiometry.append(item)
        counts.append(str(len(population_stoichiometry)))

        with open(evaluator_filepath, mode='w') as f:
            json.dump(population_stoichiometry, f)

    return counts

datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

while True:
    for data_filepath in sorted(glob.glob(os.path.join(datadir, '1488846568*.json'))):
        evaluator_filename = os.path.splitext(os.path.basename(data_filepath))[0] + "-cycles.json"
        evaluator_filepath = os.path.join(datadir, evaluator_filename)
        if not os.path.exists(evaluator_filepath):
            print(data_filepath, evaluator_filepath)
            evaluate(data_filepath, evaluator_filepath)
    time.sleep(60)
