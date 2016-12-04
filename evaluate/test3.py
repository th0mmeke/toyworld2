import json
import os


evaldir = "."
datadir = "..\data"

for filename in os.listdir(datadir):
    basename, ext = os.path.splitext(filename)
    if ext == '.json':  # in working directory
        with open(os.path.join(datadir, filename)) as f:
            reactions = json.load(f)
            print(filename, len(reactions['reactions']), len(reactions['final_population']))

for filename in os.listdir(evaldir):
    basename, ext = os.path.splitext(filename)
    if ext == ".json":
        with open(filename) as f:
            all_cycles = json.load(f)
            print(filename, len(all_cycles))