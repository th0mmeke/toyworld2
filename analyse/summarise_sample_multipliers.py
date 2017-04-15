import json
import glob
import os
import cycle_utilities
import ijson


evaldir = 'C:\Users\Thom\Dropbox\Experiments\samplerate'
#evaldir = '/home/cosc/guest/tjy17/Dropbox/Experiments/samplerate'
evaluator_filename = os.path.join(evaldir, 'cycles_and_multipliers.csv')

with open(evaluator_filename, mode='w') as f:

    f.write("Dataset, Experiment, Environment, Replicate, Sample Rate, Number of Cycles, Multiplier Species, Lineages, Average Lineage Size\n")

    for filepath in sorted(glob.glob(os.path.join(evaldir, '*p*.json'))):

        filename = os.path.splitext(os.path.basename(filepath))[0]
        nc = filename.split('-')
        filebase = '{}-{}-{}-{}-{}'.format(nc[0], nc[1], nc[2], nc[3], nc[4])
        print(filepath, filebase)

        ok = True

        number_of_cycles = 0
        try:
            with open(os.path.join(evaldir, filebase + '.json')) as f1:
                for cycle in ijson.items(f1, "item"):
                    number_of_cycles += 1
            print("Number of cycles", number_of_cycles)
        except:
            ok = False

        if ok:
            try:
                with open(os.path.join(evaldir, filebase + '-multipliers.json')) as f2:
                    clusters_by_species = json.load(f2)
                print("Clusters", len(clusters_by_species))
            except:
                ok = False


        if ok:
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