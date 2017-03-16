from collections import defaultdict


def flatten(list_of_lists):
    return [a for b in list_of_lists for a in b]


def is_reaction(node):
    return node[-1] == '>' or node[0] == '>'


def is_reaction_a(node):
    return node[-1] == '>'


def is_reaction_b(node):
    return node[0] == '>'


def map_id_to_smiles(molecule_cycle, id_to_smiles):
    smiles_cycle = []
    for step in molecule_cycle:
        new_step = step  # string so copies
        for mol in get_molecules(new_step):
            new_step = new_step.replace(mol, id_to_smiles[mol])
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
    return set(flatten([get_molecules(item) for item in cycle if is_reaction_b(item)]))


def get_reactants(cycle):
    return set(flatten([get_molecules(item) for item in cycle if is_reaction_a(item)]))


def discover_species(molecular_cycles, smiles, length=9):

    """
    Group all instances of cycles of the same species.

    :param molecular_cycles:
    :param smiles:
    :param length: Minimum length of any cycle to be considered for clustering (ignore short cycles)
    :return: dictionary with species (as frozenset) as key, and value a list of each molecular cycle of that species
    """

    species = defaultdict(list)
    for cycle in molecular_cycles:
        if len(cycle['cycle']) >= length:
            key = map_id_to_smiles(cycle['cycle'], smiles)
            species[frozenset(key)].append(cycle['cycle'])

    return species


def identify_clusters(molecular_cycles):

    """
    Add cycles into clusters of two or more cycles of any species linked by one or more molecules in common.
    Because molecules are consumed in reactions, if a molecule exists in two cycles it must be a product in one and a reactant in the other.

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
                    clusters.append([unclustered, cycle])  # TODO: should now remove unclustered from unclusters...
                    can_cluster = True

        # If still can't cluster, then add to unclustereds and hope for later...
        if not can_cluster:
            unclustereds.append(cycle)

    return clusters


def discover_multipliers(molecular_cycles_by_species):
    '''
    Stoichiometric autocatalysis occurs when a cycles has two or more linkage molecules, as products,
    to two or more different cycles.

    :param molecular_cycles_by_species: [[cycle1, cycle2]] - a list of lists of cycles for a species
    :return: List of multipliers for a species, where each multiplier is a list of cycles
    '''

    multipliers_by_species = []

    for cycles in molecular_cycles_by_species:
        multipliers_in_species = []
        clusters = identify_clusters(cycles)

        for cluster in clusters:

            products = set(flatten([get_products(cycle)-get_reactants(cycle) for cycle in cluster]))
            reactants = set(flatten([get_reactants(cycle)-get_products(cycle) for cycle in cluster]))
            linking_molecules = products.intersection(reactants)

            # Skip over clusters that are not interconnected

            if linking_molecules:

                product_cycles = [cycle for cycle in cluster if get_products(cycle).intersection(linking_molecules)]
                reactant_cycles = [cycle for cycle in cluster if get_reactants(cycle).intersection(linking_molecules)]

                print(len(reactant_cycles), len(product_cycles), len(cluster), len(linking_molecules))
                # If we have more cycles produced than consumed then autocatalytic by stoichiometry, and whole cluster is autocatalytic (as by defn of cluster, cluster is interconnected)
                if len(reactant_cycles) > len(product_cycles):
                    autocatalytic_cycles_in_cluster = [cycle for cycle in cluster if get_products(cycle).intersection(linking_molecules) or get_reactants(cycle).intersection(linking_molecules)]
                    multipliers_in_species.append(autocatalytic_cycles_in_cluster)

        if len(multipliers_in_species) > 0:
            multipliers_by_species.append(multipliers_in_species)
            print(multipliers_in_species)

    return multipliers_by_species


def discover_candidate_variable_multipliers(multipliers_by_species):

    '''
    Variable multipliers are chains of multipliers linked by common molecules, where:
    - there is some repetition of multipliers e.g., A - A - B - B - B - A - A (not tested here)
    - and the change of pattern is associated with a change in selective pressure (not tested here)

    :param multipliers: [cycles by multiplier by species]
    :return:
    '''

    candidate_variable_multiplier = []
    multipliers = flatten(multipliers_by_species)

    products = set(flatten([get_products(cycle) - get_reactants(cycle) for cycle in multipliers]))
    reactants = set(flatten([get_reactants(cycle) - get_products(cycle) for cycle in multipliers]))
    linking_molecules = products.intersection(reactants)

    if linking_molecules:

        candidate_variable_multiplier = [cycle for cycle in cluster if get_products(cycle).intersection(linking_molecules) or get_reactants(cycle).intersection(linking_molecules)]

    return candidate_variable_multiplier

