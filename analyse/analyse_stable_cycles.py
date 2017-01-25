import json
import os
from collections import defaultdict
import glob

import cycle_utilities



datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'

# Construct list of stable cycles per environment

for filename in glob.glob(os.path.join(datadir, filebase+'*molecules.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    datetime, experiment, repeat, dummy1, dummy2 = basename.split('-')
    with open(os.path.join(datadir, filename)) as f:
        all_cycles = json.load(f)
    with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, repeat, dummy1))) as f:
        state = json.load(f)
        smiles = cycle_utilities.load_smiles(state['reactions'])

    stable_states, seeds = cycle_utilities.discover_stable_cycles(all_cycles, smiles)  # [[counts per cycle type]] for this file
    print(stable_states)
    evaluator_filename = os.path.join(datadir, '{}-states.json'.format(basename))
    with open(evaluator_filename, mode='w') as f:
        json.dump(stable_states, f)
