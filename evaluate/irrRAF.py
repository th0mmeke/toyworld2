import json
import os
import glob
import random

from identify_species_cycles import IdentifySpeciesCycles

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

N_IRRRAF = 50

def compute_closure(g, foodset):
    closure = set(foodset[:])
    node_front = set(closure)
    visited = set()
    finished = False
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


def RAF(e, foodset):
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
                g.remove_node(reaction)
                finished = False
    return g


def irrRAF(e, foodset):
    #  irrRAF algorithm

    g = e.copy()  # don't modify original graph
    reactions_a = get_reactions(g)
    random.shuffle(reactions_a)

    for reaction in reactions_a:
        if g.has_node(reaction):  # case where reduced graph no longer contains this reaction from original graph
            g_copy = g.copy()
            g_copy.remove_node(reaction)
            sub_raf = RAF(g_copy, foodset)
            if sub_raf.number_of_nodes() > 0:
                g = sub_raf

    return [node for node in g.nodes() if not IdentifySpeciesCycles.is_reaction(node)]


#for filename in glob.glob(os.path.join(datadir, '1484540618-0-*-selection.json')):
for filename in glob.glob(os.path.join(datadir, '1484617345-0-0-selection.json')):

    print(filename)
    with open(filename) as f:
        state = json.load(f)

    e = IdentifySpeciesCycles(reactions=state['reactions'])

    basename, ext = os.path.splitext(filename)
    evaluator_filename = os.path.join(datadir, '{}-irrraf.json'.format(basename))

    foodset = state['initial_population'].keys()  # foodset is 'source' nodes for graph e
    with open(evaluator_filename, mode='w') as f:
        for i in range(0, N_IRRRAF):
            print("{}/{}".format(i, N_IRRRAF))
            mols_in_irrRAF = irrRAF(e.g, foodset)
            json.dump(mols_in_irrRAF, f)
