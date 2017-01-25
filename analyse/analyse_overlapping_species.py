import os
import json
import glob
from collections import defaultdict
import cycle_utilities


def discover_overlapping_cycles(cycles, smiles):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle

    :param cycles: [cycle in form of list of either cycle molecule or reactants ('r1+r2>') or products ('>p1+p2+...')]. Cycle molecules are identified by id, rather than by smiles.
    :return: {frozenset(cycle): [count]}
    '''

    # Construct lists of cycles of the same species
    species = defaultdict(list)
    for cycle in cycles:
        cycle_length = len(cycle['cycle'])
        if cycle_length > 8:
            key = cycle_utilities.map_id_to_smiles(cycle['cycle'], smiles)
            species[frozenset(key)].append(cycle['cycle'])

    stable_cycles = []
    for cycles in species.itervalues():
        clusters = cycle_utilities.identify_clusters(cycles)

        for cluster in clusters:
            if len(cluster) > 3:
                for cycle in cluster:
                    if cycle not in stable_cycles:
                        stable_cycles.append(cycle)

    clusters = cycle_utilities.identify_clusters(stable_cycles)

    print([len(x) for x in clusters])
    return clusters


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'
# filebase = '1484617345'

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

    clusters = discover_overlapping_cycles(all_cycles, smiles)  # [[counts per cycle type]] for this file
    print(json.dumps(clusters, indent=4))
    evaluator_filename = os.path.join(datadir, '{}-replicators.json'.format(basename))
    with open(evaluator_filename, mode='w') as f:
        json.dump(clusters, f)
