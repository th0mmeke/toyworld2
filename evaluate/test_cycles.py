import networkx as nx
import collections

def make_canonical(reactants):
    rle = collections.Counter(reactants)
    rle_reactants = ["{}x{}".format(count, reactant) for reactant, count in rle.iteritems()]
    return "+".join(rle_reactants)

reactions = [
    {'reactants': ['a', 'b', 'a'], 'products': ['c', 'd']},
    {'reactants': ['a', 'b', 'a'], 'products': ['c', 'c', 'f']},
    {'reactants': ['c', 'f'], 'products': ['a', 'd']}
]

# build graph

g = nx.MultiDiGraph()

for reaction in reactions:

    canonical_reactants = make_canonical(reaction['reactants'])
    canonical_products = make_canonical(reaction['products'])
    g.add_edge(canonical_reactants, canonical_products)
    for reactant in reaction['reactants']:
        g.add_edge(reactant, canonical_reactants)
    for product, count in collections.Counter(reaction['products']).iteritems():
        g.add_edge(canonical_products, product, stoichiometry=count)

# search graph

print(list(nx.simple_cycles(g)))

