import json
import os
import glob
import cycle_utilities


datadir = '/Users/Thom/Workspace/toyworld2'
filebase = '1488600998'


for filename in glob.glob(os.path.join(datadir, filebase+'*sample.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    nc = basename.split('-')
    source_basename = '{}-{}-{}-{}'.format(filebase, nc[1], nc[2], nc[3])
    with open(os.path.join(datadir, filename)) as f:
        try:
            all_cycles = json.load(f)
        except ValueError:
            pass
        else:
            with open(os.path.join(datadir, '{}.json'.format(source_basename))) as f:
                state = json.load(f)
                smiles = cycle_utilities.load_smiles(state['reactions'])

            species = cycle_utilities.discover_species(all_cycles, smiles)

            multipliers = cycle_utilities.discover_multipliers(species.itervalues())

            evaluator_filename = os.path.join(datadir, '{}-multipliers.json'.format(source_basename))
            with open(evaluator_filename, mode='w') as f:
                json.dump(multipliers, f)
