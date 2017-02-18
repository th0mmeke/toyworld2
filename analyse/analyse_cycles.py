import os
import json

from identify_molecule_cycles import IdentifyMoleculeCycles


def evaluate(filename, datadir):
    basename, ext = os.path.splitext(filename)

    for evaluator in ['molecules']:
        evaluator_filename = os.path.join(datadir, '{}-{}.json'.format(basename, evaluator))
        if not os.path.exists(evaluator_filename):

            with open(evaluator_filename, 'a'):  # fairly racy way to stop another instance evaluating same data
                os.utime(evaluator_filename, None)

            data_filename = os.path.join(datadir, filename)

            with open(data_filename) as f:
                state = json.load(f)
            e = IdentifyMoleculeCycles(reactions=state['reactions'])

            population_stoichiometry = []
            count = 0
            for reactant in e.reactants:
                count += 1
                if evaluator != 'species' or len(reactant) >= 10:
                    print("{}/{}: {}".format(count, len(e.reactants), reactant))
                    s = e.get_reactant_stoichiometry(reactant, minimum_stoichiometry=2, max_depth=7)
                    for item in s:  # item = {'cycle':cycle, 'stoichiometry': stoichiometry}
                        if evaluator == 'species':  # replace id with smiles
                            cycle = []
                            for step in item['cycle']:
                                if '+' not in step and '>' not in step and '<' not in step:
                                    cycle.append(e.get_smiles(step))
                            item['cycle'] = cycle
                        population_stoichiometry.append({'stoichiometry': item['stoichiometry'], 'cycle': item['cycle']})

            with open(evaluator_filename, mode='w') as f:
                json.dump(population_stoichiometry, f)

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

# analyse('1481398302-0-19-0.json', datadir)

import glob
for filename in glob.glob(os.path.join(datadir, '1487376025-*.json')):
    print(filename)
    evaluate(filename, datadir)
