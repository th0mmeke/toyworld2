import json
import os
import glob
import cycle_utilities


def evaluate(data_filepath, evaluator_filepath):

    with open(data_filepath) as f:
        try:
            multipliers = json.load(f)
        except ValueError:
            pass
        else:

            multipliers = cycle_utilities.discover_candidate_variable_multipliers(multipliers)
            print(multipliers)
            # with open(evaluator_filepath, mode='w') as f:
            #     json.dump(multipliers, f)

datadir = '/Users/Thom/Dropbox/Experiments'
datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

for filepath in sorted(glob.glob(os.path.join(datadir, '1489554358*multipliers.json'))):
    filename = os.path.splitext(os.path.basename(filepath))[0]
    nc = filename.split('-')
    filebase = '{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3])

    evaluator_filepath = os.path.join(datadir, filebase + "-variables.json")
    data_filepath = os.path.join(datadir, filebase + "-multipliers.json")

    if not os.path.exists(evaluator_filepath):
        print(data_filepath, evaluator_filepath)
        evaluate(data_filepath, evaluator_filepath)