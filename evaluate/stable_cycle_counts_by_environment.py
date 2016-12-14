import json
import os
import collections
import csv
from collections import Counter
from collections import defaultdict
import pickle


def get_molecules(partial_reaction_string):
    return partial_reaction_string.replace('>', '').split('+')


def get_molecules_in_cycle(cycle):
    molecules = []
    for step in cycle:
        molecules.extend(get_molecules(step))
    return molecules


def map_id_to_smiles(id_cycle, smiles):
    smiles_cycle = []
    for step in id_cycle:
        new_step = step
        for id, s in smiles.iteritems():
            new_step = new_step.replace(id, s)
        smiles_cycle.append(new_step)
    assert len(id_cycle) == len(smiles_cycle)
    return smiles_cycle


def load_smiles(reactions):
    mapping = {}
    for reaction in reactions:
        for id, smiles in reaction['reactants'].iteritems():
            mapping[id] = smiles
        for id, smiles in reaction['products'].iteritems():
            mapping[id] = smiles
    return mapping


def discover_stable_cycles(cycles, smiles):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle

    :param cycles: [cycle in form of list of either cycle molecule or reactants ('r1+r2>') or products ('>p1+p2+...')]. Cycle molecules are identified by id, rather than by smiles.
    :return: {frozenset(cycle): count}
    '''

    sorted_cycles = defaultdict(list)
    for cycle in cycles:
        sorted_cycles[len(cycle)].append(cycle)

    stable_cycles = {}

    for cycle_type, cycles_of_length in sorted_cycles.iteritems():
        # cycles_of_length has all cycles of same length, but not guaranteed to be of same type
        smiles_cycles = defaultdict(list)
        for cycle in cycles_of_length:
            s = map_id_to_smiles(cycle, smiles)
            assert len(s) == cycle_type
            smiles_cycles[frozenset(s)].append(get_molecules_in_cycle(cycle))

        # smiles_cycles now contains lists of all identical cycles types, so just have to match up the connected ones

        for cycle_type, cycles_of_type in smiles_cycles.iteritems():
            # cycle_type is the smiles form, cycles_of_type every unique cycle with molecule ids
            clusters = [set(get_molecules_in_cycle(cycles_of_type.pop()))]
            counts = defaultdict(int)
            for cycle in cycles_of_type:
                molecules = set(get_molecules_in_cycle(cycle))
                for i in range(0,len(clusters)):
                    if clusters[i].intersection(molecules):
                        clusters[i].union(molecules)
                        counts[i] += 1
                        break
            if len(counts) > 0:
                stable_cycles[frozenset(cycle_type)] = counts.values()

    return stable_cycles


datadir = "C:\Users\Thom\Dropbox\Experiments"

# Load environmental metrics

metric = collections.defaultdict(lambda: collections.defaultdict(int))

experiment = environment = '0'
with open(os.path.join(datadir, '1480963448-metadata.csv'), 'rb') as csvfile:
    r = csv.reader(csvfile, delimiter=',')
    for row in r:
        metric[int(experiment)][int(environment)] = row[-3]
        environment = str(int(environment) + 1)
        if row[0] != experiment:
            experiment = row[0]
            environment = 0

# Construct list of stable cycles per environment

for filename in os.listdir(datadir):

    basename, ext = os.path.splitext(filename)
    # Load actual cycle data
    if ext == '.json' and basename[-3:] == 'ual' and basename[:10] == '1480963448':
        print(filename)
        datetime, experiment, environment, repeat, dummy = basename.split('-')

        with open(os.path.join(datadir, filename)) as f:
            all_cycles = json.load(f)
        with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, environment, repeat))) as f:
            state = json.load(f)
            smiles = load_smiles(state['reactions'])

        stable_cycles = discover_stable_cycles(all_cycles, smiles)
        print(stable_cycles)

        with open(os.path.join(datadir, '{}-{}-{}-{}-stablestates.json'.format(datetime, experiment, environment, repeat)), 'wb') as f:
            json.dump(stable_cycles.values(), f, skipkeys=True)

# with open(os.path.join(datadir, '1480963448-counts.csv'), 'wb+') as csvfile:
#     w = csv.writer(csvfile, delimiter=',')
#     for seed, values in states.iteritems():
#         if len(values) > 3:
#             for hurst, count in values.iteritems():
#                 w.writerow([seed, hurst, count])

exit()


for experiment, experiment_metrics in metric.iteritems():
    for w in sorted(states, key=lambda x: len(states[x]), reverse=True):
        # print(experiment_number, w, np.average(stable_states[w].values()), np.std(stable_states[w].values()), stable_states[w].values())
        [metric[experiment_number][environment_number]] = states[w][environment_number]

# print("File, No. of Cycles")
# for filename in os.listdir(evaldir):
#     basename, ext = os.path.splitext(filename)
#     if ext == '.json' and (len(basename) > 7 and basename[-7:] == '-actual'):
#         with open(os.path.join(evaldir, filename)) as f:
#             all_cycles = json.load(f)
#             print(filename, len(all_cycles))