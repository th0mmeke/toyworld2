"""
Created on 6/05/2013

@author: thom
"""

import networkx as nx
import itertools
import json
import copy


class EvaluatorCycles(object):

    def __init__(self, filename="../data/toyworld2.json"):

        with open(filename) as f:
            reactions = json.load(f)

        g = nx.Graph()

        for reaction in reactions['reactions']:
            g.add_edges_from(itertools.product(reaction['reactants'], reaction['products']))
        print(g.number_of_nodes(), g.number_of_edges())

        cycle_count = {x: 0 for x in g.nodes()}
        for acs_candidate in g.nodes():
            for reactant, product in nx.bfs_edges(g, acs_candidate):
                if product == acs_candidate:
                    print(acs_candidate)  # found
                    cycle_count[acs_candidate] += 1


e = EvaluatorCycles()

