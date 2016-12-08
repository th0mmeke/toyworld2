import json
import os
import collections
import numpy as np


evaldir = "/home/cosc/guest/tjy17/Dropbox/Experiments"
datadir = "/home/cosc/guest/tjy17/Dropbox/Experiments"


stable_states = collections.defaultdict(lambda: collections.defaultdict(int))
for filename in os.listdir(evaldir):
    basename, ext = os.path.splitext(filename)
    if ext == '.json' and basename[-3:] == 'ial':
        with open(os.path.join(evaldir, filename)) as f:

            all_cycles = json.load(f)
            for cycle in all_cycles:
                stable_states[cycle['cycle'][0]][basename] += 1

for seed, states in stable_states.iteritems():
    print(seed, len(states.values()), np.average(states.values()), np.std(states.values()))

for w in sorted(stable_states, key=lambda x: len(stable_states[x]), reverse=True):
    print w, stable_states[w]

# print("File, No. of Cycles")
# for filename in os.listdir(evaldir):
#     basename, ext = os.path.splitext(filename)
#     if ext == '.json' and (len(basename) > 7 and basename[-7:] == '-actual'):
#         with open(os.path.join(evaldir, filename)) as f:
#             all_cycles = json.load(f)
#             print(filename, len(all_cycles))