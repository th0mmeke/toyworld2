import random

from identify_species_cycles import IdentifySpeciesCycles

def compute_closure(g, foodset):
    closure = set(foodset[:])
    node_front = set(closure)
    visited = set(closure)

    while len(node_front) > 0:
        new_front = []
        for node in node_front:
            successors = g.successors(node)
            for successor in g.successors(node):
                if successor not in visited:

                    if IdentifySpeciesCycles.is_reaction(successor):
                        reaction_b = g.successors(successor)[0]
                        visited.add(successor)
                        visited.add(reaction_b)
                        # all reactants must be available
                        all_reactants_available = True
                        for predecessor in g.predecessors(successor):
                            if predecessor not in closure:
                                all_reactants_available = False
                        if all_reactants_available:
                            new_front.append(reaction_b)
                    else:
                        closure.add(successor)
                        new_front.append(successor)
                        visited.add(node)

        node_front = new_front

    # closure = molecules generated from foodset under reactions in graph e
    # closure != e.reactants as e.reactants doesn't include final 'sink' products
    return closure


def get_reactions_b(g):
    return [node for node in g.nodes_iter() if IdentifySpeciesCycles.is_reaction(node) and node[0] == '>']


def get_reactions_a(g):
    return [node for node in g.nodes_iter() if IdentifySpeciesCycles.is_reaction(node) and node[-1] == '>']


def get_fgenerated(e, foodset):

    """
    Return all molecules (including foodset molecules) and reactions for the reaction graph generated from foodset
    :param e:
    :param foodset:
    :return:
    """

    g = e.copy()
    finished = False

    while not finished:
        finished = True
        closure = compute_closure(g, foodset)

        for reaction in get_reactions_b(g):
            if len(g.predecessors(reaction)) != 1:
                print("FAIL", reaction, g.predecessors(reaction))
            reaction_a = g.predecessors(reaction)[0]  # can guarantee only one predecessor
            reactants = set(g.predecessors(reaction_a))
            if len(reactants - closure) > 0:  # some reactant which isn't in the closure set
                # Now remove reaction from g and repeat
                g.remove_node(reaction)
                g.remove_node(reaction_a)

                finished = False

    closure = compute_closure(g, foodset)
    for node in g.nodes():
        if not IdentifySpeciesCycles.is_reaction(node) and node not in closure:
            g.remove_node(node)
    return g


def get_irr_fgenerated(e, foodset):
    """
    Return list of molecules and reactions in one of the possible minimal F-generated sets
    :param e:
    :param foodset:
    :return: {'molecules': all molecules except foodset, 'reactions': reactions}
    """

    g = e.copy()  # don't modify original graph
    reactions = get_reactions_b(g)
    random.shuffle(reactions)

    print(reactions)
    for reaction in reactions:
        if g.has_node(reaction):  # case where reduced graph no longer contains this reaction from original graph
            g_copy = g.copy()
            reaction_a = g.predecessors(reaction)[0]
            g_copy.remove_node(reaction)
            g_copy.remove_node(reaction_a)
            sub = get_fgenerated(g_copy, foodset)

            if len(get_reactions_b(sub)) > 0:  # still an f-generated set, check to see if can remove further
                g = sub
                print("{}{} can be removed".format(reaction_a, reaction))
            else:
                print("{}{} cannot".format(reaction_a, reaction))

    reactions = []

    molecules = [node for node in g.nodes() if not IdentifySpeciesCycles.is_reaction(node) and node not in foodset]
    return {'molecules': molecules, 'reactions': reactions}

