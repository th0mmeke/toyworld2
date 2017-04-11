import json
import glob
import os
import cycle_utilities


evaldir = 'C:\Users\Thom\Dropbox\Experiments\samplerate'
evaluator_filename = os.path.join(evaldir, 'cycles_and_multipliers.csv')

with open(evaluator_filename, mode='w') as f:

    f.write("Dataset, Experiment, Environment, Replicate, Number of Cycles, Multiplier Species, Lineages, Average Lineage Size\n")

    for filepath in sorted(glob.glob(os.path.join(evaldir, '*p*.json'))):

        filename = os.path.splitext(os.path.basename(filepath))[0]
        nc = filename.split('-')
        filebase = '{}-{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3], nc[4])
        print(filepath, filebase)

        with open(os.path.join(evaldir, filebase + '.json')) as f1:
            number_of_cycles = len(json.load(f1))

        with open(os.path.join(evaldir, filebase + '-multipliers.json')) as f2:
            clusters_by_species = json.load(f2)

            species = len(clusters_by_species)
            clusters = cycle_utilities.flatten(clusters_by_species)
            if len(clusters) == 0:
                average_entities = 0
            else:
                average_entities = sum([len(cluster) for cluster in clusters])/len(clusters)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            nc = filename.split('-')

            s = ','.join([nc[0], nc[1], nc[2], nc[3], nc[4][1:], str(number_of_cycles), str(species), str(len(clusters)), str(average_entities)])
            print(s)
            f.write(s + "\n")