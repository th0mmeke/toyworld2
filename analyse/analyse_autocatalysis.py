import json
import os
import glob
from collections import defaultdict
from identify_species_cycles import IdentifySpeciesCycles


def get_molecules(partial_reaction_string):
    return partial_reaction_string.replace('>', '').split('+')


def get_molecules_in_cycle(cycle):
    molecules = []
    for step in cycle:
        molecules.extend(get_molecules(step))
    return molecules


def map_id_to_smiles(id_cycle, smiles):
    smiles_cycle = []
    for step in id_cycle:
        new_step = step
        for id, s in smiles.iteritems():
            new_step = new_step.replace(id, s)
        smiles_cycle.append(new_step)
    assert len(id_cycle) == len(smiles_cycle)
    return smiles_cycle


def load_smiles(reactions):
    mapping = {}
    for reaction in reactions:
        for id, smiles in reaction['reactants'].iteritems():
            mapping[id] = smiles
        for id, smiles in reaction['products'].iteritems():
            mapping[id] = smiles
    return mapping

def identify_clusters(cycles):

    """
    Add cycles into clusters of two or more cycles linked by one or more molecules in common.

    :param cycles:
    :param clusters:
    :return:
    """

    clusters = []
    unclustereds = []

    for cycle in cycles:

        cycle_molecules = set(get_molecules_in_cycle(cycle))
        new_clusters = []
        can_cluster = False

        # First check if part of any existing cluster
        for cluster in clusters:
            new_cluster = cluster[:]
            for cluster_molecule in cluster:
                if set(get_molecules_in_cycle(cluster_molecule)).intersection(cycle_molecules):
                    new_cluster.append(cycle)
                    can_cluster = True
                    break
            if new_cluster:
                new_clusters.append(new_cluster)

        clusters = new_clusters

        # Now check if can form new clusters by checking against all previously unclustered cycles
        if not can_cluster:
            for unclustered in unclustereds:
                if set(get_molecules_in_cycle(unclustered)).intersection(cycle_molecules):
                    clusters.append([unclustered, cycle])
                    can_cluster = True

        # If still can't cluster, then add to unclustereds and hope for later...
        if not can_cluster:
            unclustereds.append(cycle)

    return unclustereds, clusters

def get_products(cycle):
    return set(flatten([get_molecules(item) for item in cycle if IdentifySpeciesCycles.is_reaction_b(item)]))


def get_reactants(cycle):
    return set(flatten([get_molecules(item) for item in cycle if IdentifySpeciesCycles.is_reaction_a(item)]))


def discover_species(cycles, smiles):

    species = defaultdict(list)
    for cycle in cycles:
        cycle_length = len(cycle['cycle'])
        if cycle_length > 8:
            key = map_id_to_smiles(cycle['cycle'], smiles)
            species[frozenset(key)].append(cycle['cycle'])

    return species


def flatten(list_of_lists):
    return [a for b in list_of_lists for a in b]


def discover_autocatalysis(species):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle
    '''

    autocatalytic = []
    for cycles in species.itervalues():
        unclustereds, clusters = identify_clusters(cycles)

        for cluster in clusters:
            # now look for cycle that is upstream of two or more other cycles
            # "two or more other cycles":
            #   linkages = find all molecules in two (or more) cycles
            #   find all cycles with two (or more) linkage molecules
            # "upstream":
            #   linkage molecule appears as product, not reactant in upstream cycle
            # That is, all cycles with two or more linkage molecules, as products, to different cycles

            products = set(flatten([get_products(cycle) for cycle in cluster]))
            reactants = set(flatten([get_reactants(cycle) for cycle in cluster]))
            linkages = products.intersection(reactants)
            if linkages:
                product_linkage_molecules_by_cycle = [get_products(cycle).intersection(reactants) for cycle in cluster]
                for linkage_molecules_in_cycle in product_linkage_molecules_by_cycle:
                    reactant_linkages = [list(get_reactants(cycle).intersection(linkage_molecules_in_cycle)) for cycle in cluster]
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
        smiles = load_smiles(state['reactions'])
    species = discover_species(all_cycles, smiles)
    autocatalytic = discover_autocatalysis(species)
    print(len(species), len(all_cycles))
    print(len(autocatalytic))
    evaluator_filename = os.path.join(datadir, '{}-autocatalytic.json'.format(basename))
    with open(evaluator_filename, mode='w') as f:
        json.dump(autocatalytic, f)
