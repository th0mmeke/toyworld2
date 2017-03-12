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

datadir = '/Users/Thom/Dropbox/Experiments'

# /home/cosc/guest/tjy17/Dropbox/Experiments/1488820321-0-4-0-cycles.json
# /home/cosc/guest/tjy17/Dropbox/Experiments/1488846568-2-1-0.json

for filepath in sorted(glob.glob(os.path.join(datadir, '1488846568*cycles.json'))):
    filename = os.path.splitext(os.path.basename(filepath))[0]
    nc = filename.split('-')
    filebase = '{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3])

    evaluator_filepath = os.path.join(datadir, filebase + "-multipliers.json")
    data_filepath = os.path.join(datadir, filebase + "-cycles.json")
    experiment_filepath = os.path.join(datadir, filebase + ".json")

    if not os.path.exists(evaluator_filepath):
        print(experiment_filepath, data_filepath, evaluator_filepath)
        evaluate(experiment_filepath, data_filepath, evaluator_filepath)