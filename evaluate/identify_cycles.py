import os
import json
import string

from evaluator_actual_cycles import EvaluatorActualCycles


def evaluate(filename, datadir):
    basename, ext = os.path.splitext(filename)

    if ext == '.json' and basename[-1:] in string.digits:

        for evaluator in ['actual', 'potential']:
            evaluator_filename = os.path.join(datadir, '{}-{}-sc.json'.format(basename, evaluator))
            if not os.path.exists(evaluator_filename):

                data_filename = os.path.join(datadir, filename)
                print(data_filename, evaluator_filename)
                open(evaluator_filename, mode='w').close()  # 'touch' the file to mark as in-progress to another instance - not totally safe but hopefully good enough...

                with open(data_filename) as f:
                    state = json.load(f)
                e = EvaluatorActualCycles(reactions=state['reactions'])

                # t = {x: len(e.smiles[x]) for x in e.reactants}
                # import operator
                # print(sorted(t.items(), key=operator.itemgetter(1),reverse=True))

                if evaluator != 'potential':
                    population_stoichiometry = e.get_population_stoichiometry()
                else:
                    population_stoichiometry = []
                    for item in e.get_population_stoichiometry():  # item = {'cycle':cycle, 'stoichiometry': stoichiometry}
                        # replace id with smiles
                        cycle = [e.get_smiles(step) for step in item['cycle'] if '+' not in step and '>' not in step and '<' not in step]
                        population_stoichiometry.append({'stoichiometry': item['stoichiometry'], 'cycle': cycle})

                with open(evaluator_filename, mode='w') as f:
                    json.dump(population_stoichiometry, f)

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

evaluate('1481398302-0-19-0.json', datadir)

# import glob
# for filename in glob.glob(os.path.join(datadir, '1481939843*.json')):
#     evaluate(filename, datadir)
