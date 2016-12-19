import json
import os
import collections
import csv
from itertools import chain
from collections import defaultdict


def get_metrics(filename):
    metric = collections.defaultdict(lambda: collections.defaultdict(int))

    experiment = None
    with open(os.path.join(datadir, filename), 'rb') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        for row in r:
            if experiment is None or int(row[0]) != experiment:
                experiment = int(row[0])
                environment = 0
            metric[experiment][environment] = row[-3]
            environment += 1

    return metric


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
    :return: {frozenset(cycle): [count]}
    '''

    grouped_cycles = defaultdict(list)

    for cycle in cycles:
        if type(cycle) == dict:
            grouped_cycles[len(cycle['cycle'])].append(cycle['cycle'])
        else:
            grouped_cycles[len(cycle)].append(cycle)

    stable_cycles = {}
    seeds = {}

    for cycle_type, cycles_of_length in grouped_cycles.iteritems():
        # cycles_of_length has all cycles of same length, but not guaranteed to be of same type
        smiles_cycles = defaultdict(list)
        for cycle in cycles_of_length:
            s = map_id_to_smiles(cycle, smiles)
            assert len(s) == cycle_type
            smiles_cycles[frozenset(s)].append(get_molecules_in_cycle(cycle))
            seeds[frozenset(s)] = s[0]

        # smiles_cycles now contains lists of all identical cycles types, so just have to match up the connected ones

        for cycle_type, cycles_of_type in smiles_cycles.iteritems():
            # cycle_type is the smiles form, cycles_of_type every unique cycle with molecule ids
            clusters = [set(get_molecules_in_cycle(cycles_of_type.pop()))]
            counts = defaultdict(int)
            for cycle in cycles_of_type:
                molecules = set(get_molecules_in_cycle(cycle))
                for i in range(0, len(clusters)):
                    if clusters[i].intersection(molecules):
                        clusters[i].union(molecules)
                        counts[i] += 1
                        break
            if len(counts) > 0:
                stable_cycles[frozenset(cycle_type)] = max(counts.values())  # longest stable time!

    return stable_cycles, seeds


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1481670569'

# Load environmental metrics

metric = get_metrics('1481670569-metadata.csv')  # metric[experiment][environment]

# Construct list of stable cycles per environment

with open(os.path.join(datadir, filebase + '-stablestates.csv'), 'wb') as csvfile:
    w = csv.writer(csvfile, delimiter=',')
    w.writerow(['experiment', 'environment', 'hurst', 'seeds', 'min', 'mean', 'max'])

    for filename in os.listdir(datadir):

        basename, ext = os.path.splitext(filename)

        # Load actual cycle data
        if ext == '.json' and basename[-3:] == 'ual' and basename[:len(filebase)] == filebase:
            print(filename)
            datetime, experiment, environment, repeat, dummy = basename.split('-')
            with open(os.path.join(datadir, filename)) as f:
                all_cycles = json.load(f)
            with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, environment, repeat))) as f:
                state = json.load(f)
                smiles = load_smiles(state['reactions'])

            hurst = metric[int(experiment)][int(environment)]
            stable_states, seeds = discover_stable_cycles(all_cycles, smiles)  # [[counts per cycle type]] for this file

            values_by_seed = defaultdict(list)
            for state, count in stable_states.iteritems():
                values_by_seed[seeds[state]].append(count)  # seed:longest length of stable pathway for each state
            counts_by_seed = {k:len(v) for k, v in values_by_seed.iteritems()}
            print(values_by_seed)
            print(counts_by_seed)

            if len(counts_by_seed) == 0:
                w.writerow([experiment, environment, hurst, len(counts_by_seed), 'nan', 'nan', 'nan'])
            else:
                w.writerow([experiment, environment, hurst, len(counts_by_seed), min(counts_by_seed.values()), sum(counts_by_seed.values())*1.0/len(counts_by_seed), max(counts_by_seed.values())])