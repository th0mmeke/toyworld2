import json
import os
import glob
import cycle_utilities


def evaluate(experiment_filepath, data_filepath, evaluator_filepath):

    with open(data_filepath) as f:
        try:
            all_cycles = json.load(f)
        except ValueError:
            pass
        else:
            with open(experiment_filepath) as f:
                state = json.load(f)
                smiles = cycle_utilities.load_smiles(state['reactions'])

            species = cycle_utilities.discover_species(all_cycles, smiles)
            multipliers = cycle_utilities.discover_multipliers(species.itervalues())

            with open(evaluator_filepath, mode='w') as f:
                json.dump(multipliers, f)

datadir = 'C:\Users\Thom\Dropbox\Experiments'
evaldir = 'C:\Users\Thom\Dropbox\Experiments\sampled-at-0.05'

for filepath in sorted(glob.glob(os.path.join(evaldir, '*cycles.json'))):
    filename = os.path.splitext(os.path.basename(filepath))[0]
    nc = filename.split('-')
    filebase = '{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3])

    evaluator_filepath = os.path.join(evaldir, filebase + "-multipliers.json")
    data_filepath = os.path.join(evaldir, filebase + "-cycles.json")
    experiment_filepath = os.path.join(datadir, filebase + ".json")

    if not os.path.exists(evaluator_filepath):
        print(experiment_filepath, data_filepath, evaluator_filepath)
        evaluate(experiment_filepath, data_filepath, evaluator_filepath)