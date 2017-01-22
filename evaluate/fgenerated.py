import random

from identify_species_cycles import IdentifySpeciesCycles

def compute_closure(g, foodset):
    closure = set(foodset[:])
    node_front = set(closure)
    visited = set()

    while len(node_front) > 0:
        new_front = []
        for node in node_front:
            for successor in g.successors(node):
                if successor not in visited:
                    new_front.append(successor)
                    if not IdentifySpeciesCycles.is_reaction(successor):
                        closure.add(successor)
            visited.add(node)
        node_front = new_front

    # closure = molecules generated from foodset under reactions in graph e
    # closure != e.reactants as e.reactants doesn't include final 'sink' products
    return closure


def get_reactions(g):
    return [node for node in g.nodes_iter() if IdentifySpeciesCycles.is_reaction(node) and node[0] == '>']


def get_reactions_b(g):
    return [node for node in g.nodes_iter() if IdentifySpeciesCycles.is_reaction(node) and node[-1] == '>']


def get_fgenerated(e, foodset):
    g = e.copy()
    finished = False

    while not finished:
        finished = True
        closure = compute_closure(g, foodset)
        reactions = get_reactions_b(g)

        # print("reactions", len(reactions))
        for reaction in reactions:
            reactants = set(g.predecessors(reaction))
            if len(reactants - closure) > 0:  # some reactant which isn't in the closure set
                # Now remove reaction from g and repeat
                for product in set(g.successors(reaction)) - set(foodset):  # don't remove foodset molecules
                    if len(g.predecessors(product)) == 1:  # this reaction is only source of this product
                        g.remove_node(product)
                g.remove_node(reaction)

                finished = False
    return g


def get_irr_fgenerated(e, foodset):
    #  irrRAF algorithm

    g = e.copy()  # don't modify original graph
    reactions_a = get_reactions(g)
    random.shuffle(reactions_a)

    for reaction in reactions_a:
        if g.has_node(reaction):  # case where reduced graph no longer contains this reaction from original graph
            g_copy = g.copy()
            g_copy.remove_node(reaction)
            sub_raf = get_fgenerated(g_copy, foodset)

            if len(get_reactions(sub_raf)) > 0:
                g = sub_raf

    return compute_closure(g, foodset)
