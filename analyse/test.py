import json
from identify_species_cycles import IdentifySpeciesCycles
#import Levenshtein
import itertools
import networkx as nx
import os
import glob
import random

datadir = 'C:\Users\Thom\Dropbox\Experiments'

# Construct list of stable cycles per environment


for filename in glob.glob(os.path.join(datadir, '*autocatalytic*.json')):

    print(filename)
    with open(filename) as f:
        try:
            state = json.load(f)
            print(len(state))
            for cluster in state:
                if len(cluster) > 2:
                    uniq = set(tuple(x) for x in cluster)
                    print(uniq)
        except:
            pass

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


