import json
import os
import glob
import cycle_utilities
from cycle_utilities import flatten


def discover_autocatalysis(species):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle
    '''

    autocatalytic = []
    for cycles in species.itervalues():
        unclustereds, clusters = cycle_utilities.identify_clusters(cycles)

        for cluster in clusters:
            # now look for cycle that is upstream of two or more other cycles
            # "two or more other cycles":
            #   linkages = find all molecules in two (or more) cycles
            #   find all cycles with two (or more) linkage molecules
            # "upstream":
            #   linkage molecule appears as product, not reactant in upstream cycle
            # That is, all cycles with two or more linkage molecules, as products, to different cycles

            products = set(flatten([cycle_utilities.get_products(cycle) for cycle in cluster]))
            reactants = set(flatten([cycle_utilities.get_reactants(cycle) for cycle in cluster]))
            linkages = products.intersection(reactants)
            if linkages:
                product_linkage_molecules_by_cycle = [cycle_utilities.get_products(cycle).intersection(reactants) for cycle in cluster]
                for linkage_molecules_in_cycle in product_linkage_molecules_by_cycle:
                    reactant_linkages = [list(cycle_utilities.get_reactants(cycle).intersection(linkage_molecules_in_cycle)) for cycle in cluster]
                    number_downstream_cycles = sum(map(int, [len(x) > 0 for x in reactant_linkages]))
                    # print(reactant_linkages)
                    if number_downstream_cycles > 1:  # two or more downstream cycles for single upstream -> autocatalytic
                        autocatalytic.append(reactant_linkages)

    return autocatalytic


datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'
# filebase = '1484617345'

for filename in glob.glob(os.path.join(datadir, filebase+'*molecules.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    datetime, experiment, repeat, bistate, dummy2 = basename.split('-')
    with open(os.path.join(datadir, filename)) as f:
        all_cycles = json.load(f)
    with open(os.path.join(datadir, '{}-{}-{}-bistate.json'.format(datetime, experiment, repeat))) as f:
        state = json.load(f)
        smiles = cycle_utilities.load_smiles(state['reactions'])
    species = cycle_utilities.discover_species(all_cycles, smiles)
    autocatalytic = discover_autocatalysis(species)
    print(len(species), len(all_cycles))
    print(len(autocatalytic))
    evaluator_filename = os.path.join(datadir, '{}-autocatalytic.json'.format(basename))
    with open(evaluator_filename, mode='w') as f:
        json.dump(autocatalytic, f)
