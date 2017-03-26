import json
import os
import glob
import cycle_utilities

datadir = 'C:\Users\Thom\Dropbox\Experiments'
evaldir = 'C:\Users\Thom\Dropbox\Experiments\sampled-at-0.20'

for filepath in sorted(glob.glob(os.path.join(evaldir, '*multipliers.json'))):
    filename = os.path.splitext(os.path.basename(filepath))[0]
    nc = filename.split('-')
    filebase = '{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3])

    experiment_filepath = os.path.join(datadir, filebase + ".json")

    print(experiment_filepath, filepath)

    with open(experiment_filepath) as f:
        state = json.load(f)
        smiles = cycle_utilities.load_smiles(state['reactions'])

    with open(filepath) as f2:
        try:
            clusters_by_species = json.load(f2)
        except ValueError:
            pass
        else:
            species = len(clusters_by_species)
            for clusters in clusters_by_species:
                for cluster in clusters:
                    assert len(set([len(cycle) for cycle in cluster])) == 1
                    for cycle in cluster:
                        print(cycle_utilities.map_id_to_smiles(cycle, smiles))
