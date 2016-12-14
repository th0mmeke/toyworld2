import os
import json
import string

from evaluator_actual_cycles import EvaluatorActualCycles


def evaluate(filename, datadir):
    basename, ext = os.path.splitext(filename)

    if ext == '.json' and basename[-1:] in string.digits:

        for evaluator in ['actual', 'potential']:
            evaluator_filename = os.path.join(datadir, '{}-{}.json'.format(basename, evaluator))
            if not os.path.exists(evaluator_filename):

                data_filename = os.path.join(datadir, filename)
                print(data_filename)
                with open(data_filename) as f:
                    state = json.load(f)
                e = EvaluatorActualCycles(reactions=state['reactions'])

                population_stoichiometry = []
                count = 0
                for reactant in e.reactants:
                    count += 1
                    if evaluator != 'potential' or len(reactant) >= 10:
                        print("{}/{}: {}".format(count, len(e.reactants), reactant))
                        s = e.get_reactant_stoichiometry(reactant, minimum_stoichiometry=2, max_depth=10)
                        for item in s:
                            if evaluator == 'potential':  # replace id with smiles
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

for filename in sorted(os.listdir(datadir), reverse=True):
    evaluate(filename, datadir)
