from collections import defaultdict
from identify_species_cycles import IdentifySpeciesCycles


def flatten(list_of_lists):
    return [a for b in list_of_lists for a in b]


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

    """
    Return dictionary with species (as frozenset) as key, and value a list of each molecular cycle of that species
    :param molecular_cycles:
    :param smiles:
    :param length:
    :return:
    """
    species = defaultdict(list)
    for cycle in molecular_cycles:
        if len(cycle['cycle']) >= length:
            key = map_id_to_smiles(cycle['cycle'], smiles)
            species[frozenset(key)].append(cycle['cycle'])

    return species


def identify_clusters(molecular_cycles):

    """
    Add cycles into clusters of two or more cycles linked by one or more molecules in common.

    :param molecular_cycles:
    :param clusters:
    :return:
    """

    clusters = []
    unclustereds = []

    for cycle in molecular_cycles:

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
                    clusters.append([unclustered, cycle])  # TODO: should now removed unclustered from unclusters...
                    can_cluster = True

        # If still can't cluster, then add to unclustereds and hope for later...
        if not can_cluster:
            unclustereds.append(cycle)

    return clusters


def discover_stable_cycles(cycles, smiles):
    '''
    Stable cycles are those where there are a chain of two or more cycles...
    - linked by a product in one cycle being a reactant in another
    - and where the cycles have the same "form", or sequence of molecule species (smiles) in the cycle

    :param cycles: [cycle in form of list of either cycle molecule or reactants ('r1+r2>') or products ('>p1+p2+...')]. Cycle molecules are identified by id, rather than by smiles.
    :return: {frozenset(cycle): [count]}
    '''

    grouped_cycles = defaultdict(list)

    for cycle in cycles:
        if type(cycle) == dict:
            grouped_cycles[len(cycle['cycle'])].append(cycle['cycle'])
        else:
            grouped_cycles[len(cycle)].append(cycle)

    stable_cycles = defaultdict(list)
    cycle_form = {}

    for cycles_of_length in grouped_cycles.itervalues():
        # cycles_of_length has all cycles of same length, but not guaranteed to be of same type
        smiles_cycles = defaultdict(list)
        for cycle in cycles_of_length:
            s = cycle_utilities.map_id_to_smiles(cycle, smiles)
            smiles_cycles[frozenset(s)].append(cycle_utilities.get_molecules_in_cycle(cycle))

            cycle_form[frozenset(s)] = s

        # smiles_cycles now contains lists of all identical cycles types, so just have to match up the connected ones

        for cycle_type, cycles_of_type in smiles_cycles.iteritems():
            # cycle_type is the smiles form, cycles_of_type every unique cycle with molecule ids
            clusters = [cycles_of_type.pop()]
            counts = defaultdict(int)
            for cycle_molecules in cycles_of_type:
                for i in range(0, len(clusters)):
                    if clusters[i].intersection(cycle_molecules):
                        clusters[i].union(cycle_molecules)
                        counts[i] += 1
                        break
            if len(counts) > 0 and max(counts.values()) > 1:
                stable_cycles[max(counts.values())].append(cycle_form[cycle_type])  # Longest stable duration > 1

    return stable_cycles, [x[0] for x in cycle_form.itervalues()]