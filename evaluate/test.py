import json
from evaluator_cycles import EvaluatorCycles
import networkx as nx
import string
import re

filename = '../data/toyworld2-500000.json'

with open(filename) as f:
   reactions = json.load(f)
e = EvaluatorCycles(reactions=reactions['reactions'])

#print(e.g.number_of_nodes(), e.g.number_of_edges())

# e.get_population_stoichiometry(max_depth=10, minimum_length=10, minimum_stoichiometry=3)

seed = '[H]c1n[c][c]n[c-][n-][c]n1'
rs = e.get_reactant_stoichiometry(seed, minimum_stoichiometry=4, max_depth=10)
# rs = [{'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-]=[N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-]=[N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-]=[N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[H]c1n[c][c]nc([H])n([H])[c]n1+1[c]1[c]n[c-][n-][c]n[c-]n1>', '>1[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-][N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-][N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+2[c]1[c]n[c-][n-][c]n[c-]n1', u'[c]1[c]n[c-][n-][c]n[c-]n1', '1[c]1[c]n[c-][n-][c]n[c-]n1+1[N-][N-]>', '>2[N-]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H][N][N][H]+1[H+]>', '>1[H+]+2[H][N]', u'[H][N]', '1[H][N]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-]=[N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-]=[N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-]=[N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[H]c1n[c][c]n[c-]n[c][n-]1>', '>1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-][N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-][N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>1[N-]+1[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[N-][N-]>', '>1[H+]+2[N-]', u'[N-]', '1[N-]+1[H]c1n[c][c]n[c-]n([H])[c]n1>', '>1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1+1[N-]', '[H]c1n[c][c]n[c-][n-][c]n1']}, {'stoichiometry': 4, 'cycle': ['[H]c1n[c][c]n[c-][n-][c]n1', '1[H+]+1[H]c1n[c][c]n[c-][n-][c]n1>', '>2[H+]+1[c]1[c]n[c-][n-][c]n[c-]n1', u'[H+]', '1[H+]+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+1[H]c1n[c][c]nc([H])[n-][c]n1', u'[H]c1n[c][c]nc([H])[n-][c]n1', '1[H]c1n[c][c]n[c-][n-][c]n1+1[H]c1n[c][c]nc([H])[n-][c]n1>', '>1[H+]+2[H]c1n[c][c]n[c-][n-][c]n1', '[H]c1n[c][c]n[c-][n-][c]n1']}]

re_exp = r'(?<=\d).*?(?=\+\d|>|$)'  # without $: bug >1[H]c1n[c][c]nc([H])[n-][c]n1 should produce [H]c1n[c][c]nc([H])[n-][c]n1, instead produces []
nodes = set()
for cycle in rs:
    for node in cycle['cycle']:
        nodes.add(node)
        if node[-1] == '>' or node[0] == '>':
            nodes.update(re.findall(re_exp, node))  # nodes that appear in reactions only

try:
    nodes.remove('')  # artifact of the re_exp used
except KeyError:
    pass

subg = e.g.subgraph(nodes)

if len(nodes) < len(string.ascii_letters):
    mapping = {}
    reaction_mapping = {}
    letters = list(string.ascii_letters)

    nodes.remove(seed)
    mapping[seed] = letters.pop(0)

    for node in nodes:
        if node[-1] != '>' and node[0] != '>':
            new_label = letters.pop(0)
            mapping[node] = new_label

    for node in nodes:
        if node[-1] == '>' or node[0] == '>':
            # substitute new labels for old in reaction description, and add to mapping
            new_label = node

            for old in re.findall(re_exp, new_label):
                if old != '':  # discard these - artifact of the lookahead in the re expression consuming > before matching a molecule that ends in a digit
                    if old not in mapping.keys():  # reactant from outside the cycle entering part-way through
                        new_label = letters.pop(0)
                        mapping[old] = new_label
                    new_label = string.replace(new_label, old, mapping[old])

            reaction_mapping[node] = new_label

    mapping.update(reaction_mapping)
    print(mapping)

    nx.relabel_nodes(subg, mapping, copy=False)  # remap labels in-place

nx.write_graphml(subg, "subgraph.graphml")

