import json
import os
import collections
import csv
from itertools import chain
from collections import defaultdict


def get_metrics(filebase):
    metric = collections.defaultdict(lambda: collections.defaultdict(int))

    experiment = environment = '0'
    with open(os.path.join(datadir, filebase + '-metadata.csv'), 'rb') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        for row in r:
            # print(experiment, environment, row[-3], row[-2], row[-1])
            metric[int(experiment)][int(environment)] = row[-3]
            environment = str(int(environment) + 1)
            if row[0] != experiment:
                experiment = row[0]
                environment = 0
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

    for cycle_type, cycles_of_length in grouped_cycles.iteritems():
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

    print(stable_cycles)
    return stable_cycles


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1481670569'

# Load environmental metrics

metric = get_metrics(filebase)  # metric[experiment][environment]

# Construct list of stable cycles per environment

with open(os.path.join(datadir, filebase + '-stablecounts.csv'), 'wb') as csvfile:
    w = csv.writer(csvfile, delimiter=',')

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

            stable_states = discover_stable_cycles(all_cycles, smiles).values()  # [[counts per cycle type]] for this file
            hurst = metric[int(experiment)][int(environment)]
            values = list(chain.from_iterable(stable_states))
            if len(values) == 0:
                w.writerow([experiment, environment, hurst, 'nan', 'nan', 'nan'])
            else:
                w.writerow([experiment, environment, hurst, min(values), sum(values)*1.0/len(values), max(values)])