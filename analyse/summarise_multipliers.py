import json
import glob
import os
import csv
import cycle_utilities


def get_metrics(filename):
    metadata = {}

    experiment = -1
    with open(filename, 'rb') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        next(r)  # skip header
        for row in r:
            if row[0] != experiment:
                experiment = row[0]
                environment = 0
                metadata[experiment] = {}
            metadata[experiment][str(environment)] = {'s_reactant': row[1], 's_product': row[2], 'target': row[3], 'shape': row[4], 'dfa': row[7], 'sampen': row[8]}
            environment += 1

    return metadata


datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
datadir = '/Users/Thom/Dropbox/Experiments'
filebase = '1488846568'

metadata_filename = os.path.join(datadir, filebase + "-metadata.csv")
evaluator_filename = os.path.join(datadir, 'multipliers.csv')

metadata = get_metrics(metadata_filename)

with open(evaluator_filename, mode='w') as f:

    f.write("Dataset, Experiment, Environment, Replicate, S_Reactant, S_Product, Target, Shape, DFA, Sample Entropy, Multiplier Species, Lineages, Average Lineage Size\n")

    for data_filepath in sorted(glob.glob(os.path.join(datadir, filebase+'*multipliers.json'))):

        with open(data_filepath) as f2:
            try:
                clusters_by_species = json.load(f2)
            except ValueError:
                pass
            else:
                species = len(clusters_by_species)
                clusters = cycle_utilities.flatten(clusters_by_species)
                if len(clusters) == 0:
                    average_entities = 0
                else:
                    average_entities = sum([len(cluster) for cluster in clusters])/len(clusters)
                filename = os.path.splitext(os.path.basename(data_filepath))[0]
                nc = filename.split('-')

                data = metadata[nc[1]][nc[2]]
                s = ','.join([filebase, nc[1], nc[2], nc[3], data['s_reactant'], data['s_product'], data['target'], data['shape'], data['dfa'], data['sampen'], str(species), str(len(clusters)), str(average_entities)])
                print(s)
                f.write(s + "\n")