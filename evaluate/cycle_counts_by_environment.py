import json
import os
import collections
import numpy as np
import csv


datadir = "C:\Users\Thom\Dropbox\Experiments"

metric = collections.defaultdict(lambda: collections.defaultdict(int))
experiment = environment = 0
with open(os.path.join(datadir, '1480963448-metadata.csv'), 'rb') as csvfile:
    r = csv.reader(csvfile, delimiter=',')
    for row in r:
        #print(experiment, environment, row[-3], row[-2], row[-1])
        metric[int(experiment)][int(environment)] = row[-3]
        environment = str(int(environment) + 1)
        if row[0] != experiment:
            experiment = row[0]
            environment = 0

stable_states = collections.defaultdict(lambda: collections.defaultdict(float))
states_by_environment = collections.defaultdict(lambda: collections.defaultdict(int))

for filename in os.listdir(datadir):
    basename, ext = os.path.splitext(filename)
    if ext == '.json' and basename[-3:] == 'ial':
        with open(os.path.join(datadir, filename)) as f:
            datetime, experiment, environment, repeat, dummy = basename.split('-')
            all_cycles = json.load(f)
            for cycle in all_cycles:
                stable_states[cycle['cycle'][0]][basename] += 1
                x = metric[int(experiment)][int(environment)]
                states_by_environment[cycle['cycle'][0]][x] += 1

with open(os.path.join(datadir, '1480963448-counts.csv'), 'wb+') as csvfile:
    w = csv.writer(csvfile, delimiter=',')
    for seed, values in states_by_environment.iteritems():
        if len(values) > 3:
            for hurst, count in values.iteritems():
                w.writerow([seed, hurst, count])

exit()
for experiment, experiment_metrics in metric.iteritems():
    for w in sorted(states_by_environment, key=lambda x: len(states_by_environment[x]), reverse=True):
        # print(experiment_number, w, np.average(stable_states[w].values()), np.std(stable_states[w].values()), stable_states[w].values())
        [metric[experiment_number][environment_number]] = states_by_environment[w][environment_number]

# print("File, No. of Cycles")
# for filename in os.listdir(evaldir):
#     basename, ext = os.path.splitext(filename)
#     if ext == '.json' and (len(basename) > 7 and basename[-7:] == '-actual'):
#         with open(os.path.join(evaldir, filename)) as f:
#             all_cycles = json.load(f)
#             print(filename, len(all_cycles))