import json
from evaluator_cycles import EvaluatorCycles
import Levenshtein
import itertools
import networkx as nx

filename = '../data/toyworld2-500000.json'

with open(filename) as f:
    reactions = json.load(f)
e = EvaluatorCycles(reactions=reactions['reactions'])

all_cycles = [x['cycle'] for x in e.get_population_stoichiometry(max_depth=10, minimum_length=10, minimum_stoichiometry=3)]
with open('all_cycles-500000.json', mode='w') as f:
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


