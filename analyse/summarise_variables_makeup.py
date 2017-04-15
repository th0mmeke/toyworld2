import json
import glob
import os
import csv
import cycle_utilities
import string

symbols = string.ascii_uppercase + string.ascii_lowercase + string.digits
mapping = {}

def map_cycle_to_symbol(cycle):
    global symbols
    global mapping

    if tuple(cycle) not in mapping:
        mapping[tuple(cycle)] = symbols[0]
        symbols = symbols[1:]
    # print(cycle)

    return mapping[tuple(cycle)]


datadir = 'C:\Users\Thom\Dropbox\Experiments'
evaldir = 'C:\Users\Thom\Dropbox\Experiments\sampled-at-0.20'

evaluator_filename = os.path.join(evaldir, 'variable-makeup.csv')

with open(evaluator_filename, mode='w') as f1:
    f1.write("Dataset, Experiment, Environment, Replicate, Variable Multiplier\n")
    for filepath in sorted(glob.glob(os.path.join(evaldir, '*variables.json'))):

        filename = os.path.splitext(os.path.basename(filepath))[0]
        nc = filename.split('-')
        filebase = '{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3])

        experiment_filepath = os.path.join(datadir, filebase + ".json")

        with open(experiment_filepath) as f2:
            state = json.load(f2)
            smiles = cycle_utilities.load_smiles(state['reactions'])

        with open(filepath) as f3:
            try:
                variables = json.load(f3)
            except ValueError:
                pass
            else:
                if len(variables) > 0:
                    filename = os.path.splitext(os.path.basename(filepath))[0]
                    nc = filename.split('-')
                    for variable in variables:
                        symbolic_variable = [map_cycle_to_symbol(cycle_utilities.map_id_to_smiles(cycle, smiles)) for cycle in variable]
                        s = ','.join([nc[0], nc[1], nc[2], nc[3], ' '.join(symbolic_variable)])
                        # print(s)
                        f1.write(s + "\n")
