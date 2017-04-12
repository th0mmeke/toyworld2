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
            if len(row) <= 7:
                row.extend([str(0), str(0)])

            metadata[experiment][str(environment)] = {'s_reactant': row[1], 's_product': row[2], 'target': row[3], 'shape': row[4], 'dfa': row[7], 'sampen': row[8]}
            environment += 1

    return metadata


datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
evaldir = '/home/cosc/guest/tjy17/Dropbox/Experiments/Sampled-at-0.20'

#datadir = '/Users/Thom/Dropbox/Experiments'

evaluator_filename = os.path.join(evaldir, 'multipliers.csv')

current_metadata_filename = None
with open(evaluator_filename, mode='w') as f:

    f.write("Dataset, Experiment, Environment, Replicate, S_Reactant, S_Product, Target, Shape, DFA, Sample Entropy, Multiplier Species, Lineages, Average Lineage Size\n")

    for data_filepath in sorted(glob.glob(os.path.join(evaldir, '*multipliers.json'))):

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

                metadata_filename = os.path.join(datadir, nc[0] + "-metadata.csv")
                if metadata_filename != current_metadata_filename:
                    metadata = get_metrics(metadata_filename)
                    current_metadata_filename = metadata_filename

                data = metadata[nc[1]][nc[2]]
                s = ','.join([nc[0], nc[1], nc[2], nc[3], data['s_reactant'], data['s_product'], data['target'], data['shape'], data['dfa'], data['sampen'], str(species), str(len(clusters)), str(average_entities)])
                print(s)
                f.write(s + "\n")