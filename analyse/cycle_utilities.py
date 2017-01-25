from collections import defaultdict
from identify_species_cycles import IdentifySpeciesCycles


def map_id_to_smiles(molecule_cycle, smiles):
    smiles_cycle = []
    for step in molecule_cycle:
        new_step = step
        for id, s in smiles.iteritems():
            new_step = new_step.replace(id, s)
        smiles_cycle.append(new_step)
    assert len(molecule_cycle) == len(smiles_cycle)
    return smiles_cycle


def load_smiles(reactions):
    mapping = {}
    for reaction in reactions:
        for id, smiles in reaction['reactants'].iteritems():
            mapping[id] = smiles
        for id, smiles in reaction['products'].iteritems():
            mapping[id] = smiles
    return mapping


def get_molecules(partial_reaction_string):
    return partial_reaction_string.replace('>', '').split('+')


def get_molecules_in_cycle(cycle):
    """
    :param cycle:
    :return: {all molecules in cycle}
    """
    molecules = set()
    for step in cycle:
        molecules = molecules.union(set(get_molecules(step)))
    return molecules


def get_products(cycle):
    return set(flatten([get_molecules(item) for item in cycle if IdentifySpeciesCycles.is_reaction_b(item)]))


def get_reactants(cycle):
    return set(flatten([get_molecules(item) for item in cycle if IdentifySpeciesCycles.is_reaction_a(item)]))


def discover_species(molecular_cycles, smiles, length=9):

    species = defaultdict(list)
    for cycle in molecular_cycles:
        if len(cycle['cycle']) >= length:
            key = map_id_to_smiles(cycle['cycle'], smiles)
            species[frozenset(key)].append(cycle['cycle'])

    return species


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

        cycle_molecules = get_molecules_in_cycle(cycle)
        new_clusters = []
        can_cluster = False

        # First check if part of any existing cluster
        for cluster in clusters:
            new_cluster = cluster[:]
            for cluster_molecule in cluster:
                if get_molecules_in_cycle(cluster_molecule).intersection(cycle_molecules):
                    new_cluster.append(cycle)
                    can_cluster = True
                    break
            if new_cluster:
                new_clusters.append(new_cluster)

        clusters = new_clusters

        # Now check if can form new clusters by checking against all previously unclustered cycles
        if not can_cluster:
            for unclustered in unclustereds:
                if get_molecules_in_cycle(unclustered).intersection(cycle_molecules):
                    clusters.append([unclustered, cycle])
                    can_cluster = True

        # If still can't cluster, then add to unclustereds and hope for later...
        if not can_cluster:
            unclustereds.append(cycle)

    return unclustereds, clusters


def flatten(list_of_lists):
    return [a for b in list_of_lists for a in b]

