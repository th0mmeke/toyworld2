import json
import glob
import os
import csv


def get_metrics(filename):
    s_reactant = {}
    s_product = {}
    target = {}
    shape = {}

    experiment = None
    with open(filename, 'rb') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        for row in r:
            if experiment is None or int(row[0]) != experiment:
                experiment = row[0]
                s_reactant[experiment] = row[1]
                s_product[experiment] = row[2]
                target[experiment] = row[3]
                shape[experiment] = row[4]

    return s_reactant, s_product, target, shape


datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1488846568'
metadata_filename = os.path.join(datadir, filebase + "-metadata.csv")
evaluator_filename = os.path.join(datadir, 'multipliers.csv')

s_reactant, s_product, target, shape = get_metrics(metadata_filename)

with open(evaluator_filename, mode='w') as f:

    f.write("Datetime, Experiment, Environment, Replicate, S_Reactant, S_Product, Target, Shape, Multipliers\n")

    for data_filepath in sorted(glob.glob(os.path.join(datadir, filebase+'*multipliers.json'))):

        with open(data_filepath) as f2:
            try:
                clusters = json.load(f2)
            except ValueError:
                pass
            else:
                filename = os.path.splitext(os.path.basename(data_filepath))[0]
                nc = filename.split('-')
                s = ','.join([filebase, nc[1], nc[2], nc[3], s_reactant[nc[1]], s_product[nc[1]], target[nc[1]], shape[nc[1]], str(len(clusters))])
                f.write(s + "\n")