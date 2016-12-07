import json
import os


evaldir = "C:\Users\Thom\Dropbox/Experiments"
datadir = "C:\Users\Thom\Dropbox/Experiments"

print("File, No. of reactions, Size of final population")
for filename in os.listdir(datadir):
    basename, ext = os.path.splitext(filename)
    if ext == '.json' and  (len(basename) > 7 and not basename[-7:] == '-actual'):
        with open(os.path.join(datadir, filename)) as f:
            reactions = json.load(f)
            print(filename, len(reactions['reactions']), len(reactions['final_population']))

print("File, No. of Cycles")
for filename in os.listdir(evaldir):
    basename, ext = os.path.splitext(filename)
    if ext == '.json' and (len(basename) > 7 and basename[-7:] == '-actual'):
        with open(os.path.join(evaldir, filename)) as f:
            all_cycles = json.load(f)
            print(filename, len(all_cycles))