import os
import json
import random

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(filename, datadir):
    basename, ext = os.path.splitext(filename)

    evaluator = 'molecules'
    evaluator_filename = os.path.join(datadir, '{}-{}sample.json'.format(basename, evaluator))
    data_filename = os.path.join(datadir, filename)
    population_stoichiometry = []
    test = []

    with open(data_filename) as f:
        state = json.load(f)
    try:
        e = IdentifyMoleculeCycles(reactions=state['reactions'])
    except:
        pass
    else:
        print(len(e.reactants), len(e.reactant_molecules))
        sample = random.sample(e.reactant_molecules, len(e.reactant_molecules)/10)
        for mol in sample:
            for item in e.sample_reactant_stoichiometry(mol):
                [p.append({'stoichiometry': item['stoichiometry'], 'cycle': item['cycle']})
                test.append(item)
        print(len(test), len(population_stoichiometry))

    # with open(evaluator_filename, mode='w') as f:
    #     json.dump(population_stoichiometry, f)


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
datadir='..'
# analyse('1481398302-0-19-0.json', datadir)

import glob
x = os.path.join('/Users/tjy17/Workspace/toyworld2', '1488595072*.json')
print(x)
for filename in glob.glob(x):
    print(filename)
    evaluate(filename, datadir)
