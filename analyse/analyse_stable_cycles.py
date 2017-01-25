import json
import os
from collections import defaultdict
import cycle_utilities


def discover_stable_cycles(cycles, smiles):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle

    :param cycles: [cycle in form of list of either cycle molecule or reactants ('r1+r2>') or products ('>p1+p2+...')]. Cycle molecules are identified by id, rather than by smiles.
    :return: {frozenset(cycle): [count]}
    '''

    grouped_cycles = defaultdict(list)

    for cycle in cycles:
        if type(cycle) == dict:
            grouped_cycles[len(cycle['cycle'])].append(cycle['cycle'])
        else:
            grouped_cycles[len(cycle)].append(cycle)

    stable_cycles = defaultdict(list)
    cycle_form = {}

    for cycle_type, cycles_of_length in grouped_cycles.iteritems():
        # cycles_of_length has all cycles of same length, but not guaranteed to be of same type
        smiles_cycles = defaultdict(list)
        for cycle in cycles_of_length:
            s = cycle_utilities.map_id_to_smiles(cycle, smiles)
            assert len(s) == cycle_type
            smiles_cycles[frozenset(s)].append(cycle_utilities.get_molecules_in_cycle(cycle))

            cycle_form[frozenset(s)] = s

        # smiles_cycles now contains lists of all identical cycles types, so just have to match up the connected ones

        for cycle_type, cycles_of_type in smiles_cycles.iteritems():
            # cycle_type is the smiles form, cycles_of_type every unique cycle with molecule ids
            clusters = [set(cycle_utilities.get_molecules_in_cycle(cycles_of_type.pop()))]
            counts = defaultdict(int)
            for cycle in cycles_of_type:
                molecules = set(cycle_utilities.get_molecules_in_cycle(cycle))
                for i in range(0, len(clusters)):
                    if clusters[i].intersection(molecules):
                        clusters[i].union(molecules)
                        counts[i] += 1
                        break
            if len(counts) > 0 and max(counts.values()) > 1:
                stable_cycles[max(counts.values())].append(cycle_form[cycle_type])  # Longest stable duration > 1

    return stable_cycles, [x[0] for x in cycle_form.itervalues()]

import glob

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

    stable_states, seeds = discover_stable_cycles(all_cycles, smiles)  # [[counts per cycle type]] for this file
    print(stable_states)
    evaluator_filename = os.path.join(datadir, '{}-states.json'.format(basename))
    with open(evaluator_filename, mode='w') as f:
        json.dump(stable_states, f)
