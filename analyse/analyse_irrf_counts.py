import json
import glob
import os

datadir = 'C:\Users\Thom\Dropbox\Experiments'

filenames = glob.glob(os.path.join(datadir, '*irrraf.json'))
datetimes = set([os.path.basename(filename).split('-')[0] for filename in filenames])
print(filenames)
evaluator_filename = os.path.join(datadir, 'irrf.csv')
print(evaluator_filename)


with open(evaluator_filename, mode='w') as f:
    f.write("Datetime, Experiment, Environment, Replicate, Irrf\n")
    for datetime in datetimes:

        for filename in glob.glob(os.path.join(datadir, datetime + '*irrraf.json')):

            basename, ext = os.path.splitext(filename)
            print(filename)
            nc = os.path.basename(basename).split('-')
            if datetime == '1484540618':
                experiment = str(0)
                environment = str(0)
                replicate = nc[2]
            elif datetime == '1481939843':
                experiment = nc[2]
                environment = nc[3]
                replicate = nc[4]
            elif datetime == '1481952255':
                experiment = nc[3]
                environment = nc[4]
                replicate = nc[5]
            else:
                experiment = nc[1]
                environment = nc[2]
                replicate = nc[3]

            with open(filename) as f2:
                try:
                    irrf = json.load(f2)
                except:
                    pass
                else:
                    print(len(irrf))
            #         s = ','.join([datetime, experiment, environment, replicate, str(len(irrf))])
            #         f.write(s + "\n")


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


