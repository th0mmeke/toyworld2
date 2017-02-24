import json
import os
import collections
import csv
from collections import defaultdict
import cycle_utilities

def get_metrics(filename):
    metric = collections.defaultdict(lambda: collections.defaultdict(int))

    experiment = None
    with open(os.path.join(datadir, filename), 'rb') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        for row in r:
            if experiment is None or int(row[0]) != experiment:
                experiment = int(row[0])
                environment = 0
            metric[experiment][environment] = row[-2]
            environment += 1

    return metric


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'

# Load environmental metrics

metric = get_metrics(filebase+'-metadata.csv')  # metric[experiment][environment]

# Construct list of stable cycles per environment

with open(os.path.join(datadir, filebase + '-number-stablestates.csv'), 'wb') as csvfile:
    w = csv.writer(csvfile, delimiter=',')
    w.writerow(['experiment', 'environment', 'repeat', 'dfa', 'seeds', 'min', 'mean', 'max'])

    for filename in os.listdir(datadir):

        basename, ext = os.path.splitext(filename)

        # Load actual cycle data
        if ext == '.json' and basename[-3:] == 'ule' and basename[:len(filebase)] == filebase:
            print(filename)
            datetime, experiment, environment, repeat, dummy2 = basename.split('-')
            with open(os.path.join(datadir, filename)) as f:
                all_cycles = json.load(f)
            with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, environment, repeat))) as f:
                state = json.load(f)
                smiles = cycle_utilities.load_smiles(state['reactions'])

            dfa = metric[int(experiment)][int(environment)]
            stable_states, seeds = cycle_utilities.discover_multipliers(all_cycles, smiles)  # [[counts per cycle type]] for this file

            values_by_seed = defaultdict(list)
            for state, count in stable_states.iteritems():
                values_by_seed[seeds[state]].append(count)  # seed:longest length of stable pathway for each state
            counts_by_seed = {k: len(v) for k, v in values_by_seed.iteritems()}  # counts_by_seed = seed: number of states

            # counts_by_seed is number of states for each seed
            if len(counts_by_seed) == 0:
                w.writerow([experiment, environment, repeat, dfa, len(counts_by_seed), 'nan', 'nan', 'nan'])
            else:
                w.writerow([experiment, environment, repeat, dfa, len(counts_by_seed), min(counts_by_seed.values()), sum(counts_by_seed.values()) * 1.0 / len(counts_by_seed), max(counts_by_seed.values())])