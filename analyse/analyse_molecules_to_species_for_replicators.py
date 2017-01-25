import json
import os
import glob
import cycle_utilities

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'

# Construct list of stable cycles per environment

for filename in glob.glob(os.path.join(datadir, filebase+'*replicators.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    datetime, experiment, repeat, bistate, dummy2, dummy3 = basename.split('-')
    with open(os.path.join(datadir, filename)) as f:
        replicators = json.load(f)
    with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, repeat, bistate))) as f:
        state = json.load(f)
        smiles = cycle_utilities.load_smiles(state['reactions'])

    species_replicators = []
    for replicator in replicators:
        species_cycles = []
        for molecular_cycle in replicator:
            species_cycles.append(cycle_utilities.map_id_to_smiles(molecular_cycle, smiles))
        species_replicators.append(species_cycles)
    # print(species_replicator)
    # print([len(x) for x in replicators], [len(x) for x in species_replicators])

    evaluator_filename = os.path.join(datadir, '{}-species.json'.format(basename))
    print(evaluator_filename)
    with open(evaluator_filename, mode='w') as f:
        json.dump(species_replicators, f)