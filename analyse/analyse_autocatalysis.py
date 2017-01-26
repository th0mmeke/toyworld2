import json
import os
import glob
import cycle_utilities


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'
# filebase = '1484617345'

for filename in glob.glob(os.path.join(datadir, filebase+'*molecules.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    datetime, experiment, repeat, bistate, molecules = basename.split('-')

    with open(os.path.join(datadir, filename)) as f:
        all_cycles = json.load(f)
    with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, repeat, bistate))) as f:
        state = json.load(f)
        smiles = cycle_utilities.load_smiles(state['reactions'])

    species = cycle_utilities.discover_species(all_cycles, smiles)

    autocatalytic = cycle_utilities.discover_autocatalysis(species.itervalues())

    evaluator_filename = os.path.join(datadir, '{}-autocatalytic.json'.format(basename))
    with open(evaluator_filename, mode='w') as f:
        json.dump(autocatalytic, f)
