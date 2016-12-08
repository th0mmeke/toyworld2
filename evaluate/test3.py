import json
import os
import collections


evaldir = "C:\Users\Thom\Dropbox/Experiments"
datadir = "C:\Users\Thom\Dropbox/Experiments"



print("File, No. of Cycles")
for filename in os.listdir(evaldir):
    basename, ext = os.path.splitext(filename)
    if ext == '.json' and basename[-3:] == 'ial':
        with open(os.path.join(evaldir, filename)) as f:
            stable_states = collections.defaultdict(int)
            all_cycles = json.load(f)
            print(filename, len(all_cycles))
            for cycle in all_cycles:
                #print(cycle['stoichiometry'], cycle['cycle'])
                stable_states[cycle['cycle'][0]] += 1

            print(sum(stable_states.values())/len(stable_states), stable_states)

# print("File, No. of Cycles")
# for filename in os.listdir(evaldir):
#     basename, ext = os.path.splitext(filename)
#     if ext == '.json' and (len(basename) > 7 and basename[-7:] == '-actual'):
#         with open(os.path.join(evaldir, filename)) as f:
#             all_cycles = json.load(f)
#             print(filename, len(all_cycles))