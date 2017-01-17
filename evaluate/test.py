import json
from evaluator_cycles import EvaluatorCycles
#import Levenshtein
import itertools
import networkx as nx
import os
import glob
import random

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'

# Construct list of stable cycles per environment


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
                    if not EvaluatorCycles.is_reaction(successor):
                        closure.add(successor)
            visited.add(node)
        node_front = new_front

    # closure = molecules generated from foodset under reactions in graph e
    # closure != e.reactants as e.reactants doesn't include final 'sink' products
    return closure


def get_reactions(g):
    return [node for node in g.nodes_iter() if EvaluatorCycles.is_reaction(node) and node[0] == '>']

def get_reactions_b(g):
    return [node for node in g.nodes_iter() if EvaluatorCycles.is_reaction(node) and node[-1] == '>']

def RAF(e, foodset):
    g = e.copy()
    finished = False

    while not finished:
        finished = True
        closure = compute_closure(g, foodset)
        reactions = get_reactions_b(g)
        # print("reactions", len(reactions))
        for reaction in reactions:
            #reactants = reaction.replace('>', '').split('+')  # must remove stoichiometry...
            reactants = set(g.predecessors(reaction))
            if len(reactants - closure) > 0:  # some reactant which isn't in the closure set
                # Now remove reaction from g and repeat
                g.remove_node(reaction)
                finished = False
    return g


def irrRAF(e, foodset):
    #  irrRAF algorithm
    reactions_a = get_reactions(e)

    print("Original", len([node for node in e.nodes() if not EvaluatorCycles.is_reaction(node)]))

    for i in range(0, 10):
        g = e.copy()
        random.shuffle(reactions_a)
        for reaction in reactions_a:
            if g.has_node(reaction):  # case where reduced graph no longer contains this reaction from original graph
                g_copy = g.copy()
                g_copy.remove_node(reaction)
                sub_raf = RAF(g_copy, foodset)
                if sub_raf.number_of_nodes() > 0:
                    g = sub_raf
        irrRAF = [node for node in g.nodes() if not EvaluatorCycles.is_reaction(node)]
        print("irrRAF", len(irrRAF))

for filename in glob.glob(os.path.join(datadir, '1484540618-0-*-selection.json')):
#for filename in glob.glob(os.path.join(datadir, '1484617345-0-0-selection.json')):

    print(filename)
    with open(filename) as f:
        state = json.load(f)
    foodset = state['initial_population'].keys()  # foodset is 'source' nodes for graph e

    e = EvaluatorCycles(reactions=state['reactions'])

    irrRAF(e.g, foodset)



exit()

population_stoichiometry = []
count = 0
for reactant in e.reactants:
    count += 1
    if len(reactant) >= 10:
        print("{}/{}: {}".format(count, len(e.reactants), reactant))
        population_stoichiometry.extend(e.get_reactant_stoichiometry(reactant, minimum_stoichiometry=3, max_depth=10))

# e.get_population_stoichiometry(max_depth=10, minimum_length=10, minimum_stoichiometry=3)
all_cycles = [x['cycle'] for x in population_stoichiometry]
with open('1480307656-0-0-0-actual.json', mode='w') as f:
    json.dump(all_cycles, f)

# with open('all_cycles.json', mode='r') as f:
#     all_cycles = json.load(f)

for cycle in all_cycles:
    print(cycle[0], e.get_food_set(cycle))

exit()


with open(filename) as f:
    reactions = json.load(f)
e = EvaluatorCycles(reactions=reactions['reactions'])

#e.get_population_stoichiometry(minimum_stoichiometry=2, max_depth=8, minimum_length=10)

all_cycles = [x['cycle'] for x in e.get_population_stoichiometry(max_depth=8, minimum_length=10, minimum_stoichiometry=3)]
with open('all_cycles.json', mode='w') as f:
    json.dump(all_cycles, f)

for cycle in all_cycles:
    print(e.get_food_set(cycle))

exit()

print(nx.dfs_tree(e.g, '[H]c1n[c][c]n[c-][n-][c]n1').nodes())
g = nx.Graph(e.g)
for seed in g.nodes():
    for node in nx.dfs_tree(g, seed):
        if node[0] != '>' and node[-1] != '<':  # molecule
            if node not in closure_set:
                # trace back to producing reaction, and add reactants
                pass


# closure_set =  all products of reactions in cycle + food set
# reactant_set = all reactants in cycle
# while can pop molecule from reactant_set
#    reactant_set.extend(get_reactants(molecule))

def get_reactants(molecule, closure_set):
    reactant_set = set()
    if molecule not in closure_set:
        # find producing reaction - what if more than one?
        for reactant in reactions:
            reactant_set.extend(get_reactants(reactant))
    return reactant_set


g = nx.DiGraph(e.g)
print(nx.number_strongly_connected_components(g))
exit()
# population_stoichiometry = [{'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-]=[N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-]=[N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-]=[N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[H]c1n[c][c]nc([H])n([H])[c]n1+1[c]1[c]n[c-][n-][c]n[c-]n1>', '>1[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-][N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-][N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-][N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-]=[N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-]=[N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-]=[N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[H]c1n[c][c]n[c-]n[c][n-]1>', '>1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-][N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-][N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-][N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}]
all_cycles = [x['cycle'] for x in e.get_population_stoichiometry(max_depth=10, minimum_length=10, minimum_stoichiometry=3)]

molecules = []
for cycle in all_cycles:
    molecules.extend([unicode(x) for x in cycle if x[0] != '>' and x[-1] != '>'])
print(molecules)

for a, b in itertools.combinations(molecules, 2):
    if len(a) > 10 and len(b) > 10:
        d = Levenshtein.distance(a, b)
        if 0 < d < 5:
            print(d, a, b)

exit()


filename = '../data/toyworld2-500000.json'

with open(filename) as f:
   reactions = json.load(f)
e = EvaluatorCycles(reactions=reactions['reactions'])


seed = '[H]c1n[c][c]n[c-][n-][c]n1'
cycles = e.get_reactant_stoichiometry(seed, minimum_stoichiometry=4, max_depth=10)


e.write_subset_to_graphml(seed, [x['cycle'] for x in cycles])


