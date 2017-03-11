import glob
import json
import os

import fgenerated_utilities
from identify_species_cycles import IdentifySpeciesCycles

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

N_IRRRAF = 5

for filename in glob.glob(os.path.join(datadir, '1484540618*bistate.json')):

    print(filename)
    with open(filename) as f:
        state = json.load(f)

    e = IdentifySpeciesCycles(reactions=state['reactions'])

    basename, ext = os.path.splitext(filename)
    evaluator_filename = os.path.join(datadir, '{}-irrraf.json'.format(basename))

    foodset = state['initial_population'].keys()  # foodset is 'source' nodes for graph e
    with open(evaluator_filename, mode='w', buffering=0) as f:
        f.write("[")
        for i in range(0, N_IRRRAF):
            print("{}/{}".format(i, N_IRRRAF))
            irr_f = fgenerated_utilities.get_irr_fgenerated(e.g, foodset)
            json.dump(irr_f, f)
            if i < N_IRRRAF-1:
                f.write(",")
        f.write("]")

