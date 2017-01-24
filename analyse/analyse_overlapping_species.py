import json
import os
import collections
import csv
import itertools

from collections import defaultdict


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

def discover_stable_cycles(cycles, smiles):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle

    :param cycles: [cycle in form of list of either cycle molecule or reactants ('r1+r2>') or products ('>p1+p2+...')]. Cycle molecules are identified by id, rather than by smiles.
    :return: {frozenset(cycle): [count]}
    '''

    # Construct lists of cycles of the same species
    species = defaultdict(list)
    for cycle in cycles:
        cycle_length = len(cycle['cycle'])
        if cycle_length > 8:
            key = map_id_to_smiles(cycle['cycle'], smiles)
            species[frozenset(key)].append(cycle['cycle'])

    stable_cycles = []
    for cycles in species.itervalues():
        unclustereds, clusters = identify_clusters(cycles)

        for cluster in clusters:
            if len(cluster) > 3:
                for cycle in cluster:
                    if cycle not in stable_cycles:
                        stable_cycles.append(cycle)

    unclustereds, clusters = identify_clusters(stable_cycles)

    print([len(x) for x in clusters])
    return clusters

import glob

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'
# filebase = '1484617345'

# Construct list of stable cycles per environment

for filename in glob.glob(os.path.join(datadir, filebase+'*molecules.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    # datetime, experiment, repeat, dummy1, dummy2 = basename.split('-')
    # with open(os.path.join(datadir, filename)) as f:
    #     all_cycles = json.load(f)
    # with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, repeat, dummy1))) as f:
    #     state = json.load(f)
    #     smiles = load_smiles(state['reactions'])
    #
    # clusters = discover_stable_cycles(all_cycles, smiles)  # [[counts per cycle type]] for this file
    # print(json.dumps(clusters, indent=4))
    # evaluator_filename = os.path.join(datadir, '{}-replicators.json'.format(basename))
    # with open(evaluator_filename, mode='w') as f:
    #     json.dump(clusters, f)
