import json
import glob
import os
import collections
import csv


def get_metrics(filename):
    s_reactant = {}
    s_product = {}
    experiment = None
    with open(filename, 'rb') as csvfile:
        r = csv.reader(csvfile, delimiter=',')
        for row in r:
            if experiment is None or int(row[0]) != experiment:
                experiment = row[0]
                s_reactant[experiment] = row[1]
                s_product[experiment] = row[2]

    return s_reactant, s_product


def get_data(datadir, file_regexp):

    data_strings = []
    for filename in glob.glob(os.path.join(datadir, file_regexp)):

        basename, ext = os.path.splitext(filename)
        nc = os.path.basename(basename).split('-')

        dataset = nc[0]
        if dataset == '1481939843':
            experiment = nc[2]
            environment = nc[3]
            replicate = nc[4]
            basedata_filename = '-'.join([dataset, 'energy', experiment, environment, replicate]) + '.json'
        else:
            experiment = nc[1]
            environment = nc[2]
            replicate = nc[3]
            basedata_filename = '-'.join([dataset, experiment, environment, replicate]) + '.json'

        metadata_filename = os.path.join(datadir, nc[0] + "-metadata.csv")
        s_reactant, s_product = get_metrics(metadata_filename)

        with open(os.path.join(datadir, basedata_filename)) as f2:
            try:
                clusters = json.load(f2)
            except ValueError:
                pass
            else:
                s = ','.join([dataset, experiment, environment, replicate, s_reactant[experiment], s_product[experiment], str(len(clusters))])
                data_strings.append(s)

    return data_strings

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

evaluator_filename = os.path.join(datadir, 'multipliers.csv')

with open(evaluator_filename, mode='w') as f:
    f.write("Datetime, Experiment, Environment, Replicate, S_Reactant, S_Product, Multipliers\n")

    for s in get_data(datadir, '*multipliers.json'):
        # f.write(s + '\n')
        print(s)

exit()


# filename = '../data/toyworld2-500000.json'
#
# with open(filename) as f:
#    reactions = json.load(f)
# e = IdentifySpeciesCycles(reactions=reactions['reactions'])
#
#
# seed = '[H]c1n[c][c]n[c-][n-][c]n1'
# cycles = e.get_reactant_stoichiometry(seed, minimum_stoichiometry=4, max_depth=10)
#
#
# e.write_subset_to_graphml(seed, [x['cycle'] for x in cycles])


