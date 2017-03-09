import json
import glob
import os
import csv


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
#datadir = '/Users/Thom/Dropbox/Experiments'

evaluator_filename = os.path.join(datadir, 'metadata.csv')

with open(evaluator_filename, mode='w') as f:

    f.write("Dataset, Experiment, Environment, S_Reactant, S_Product, Target, Shape, DFA, Sample Entropy\n")

    for data_filepath in sorted(glob.glob(os.path.join(datadir, '*-metadata.csv'))):

        with open(data_filepath) as f2:
            metadata = get_metrics(data_filepath)
            filename = os.path.splitext(os.path.basename(data_filepath))[0]
            datetime = filename.split('-')[0]
            for experiment in sorted(metadata.keys()):
                for environment in sorted(metadata[experiment].keys()):
                    data = metadata[experiment][environment]
                    s = ','.join([datetime, experiment, environment, data['s_reactant'], data['s_product'], data['target'], data['shape'], data['dfa'], data['sampen']])
                    print(s)
                    f.write(s + "\n")