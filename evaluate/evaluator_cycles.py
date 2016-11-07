import networkx as nx
import itertools
import collections


class EvaluatorCycles(object):

    def __init__(self, reactions):

        self.cycle_stoichiometry = []

        # Build graph

        g = nx.MultiDiGraph()

        for reaction in reactions:

            canonical_reactants = EvaluatorCycles.make_canonical(reaction['reactants'])
            canonical_products = EvaluatorCycles.make_canonical(reaction['products'])
            g.add_edge(canonical_reactants, canonical_products)
            for reactant in reaction['reactants']:
                g.add_edge(reactant, canonical_reactants)
            for product, count in collections.Counter(reaction['products']).iteritems():
                g.add_edge(canonical_products, product, stoichiometry=count)

        # Search graph

        for cycle in nx.simple_cycles(g):
            stoichiometry = 1
            for start_node, end_node in zip(cycle[0:-1], cycle[1:]):
                # get the max stoichiometry - looking for ACS of high stoichiometry, assume take the highest reaction pathway
                stoichiometries = [edge['stoichiometry'] for edge in g.edge[start_node][end_node].itervalues() if 'stoichiometry' in edge.keys()]
                if len(stoichiometries) > 0:
                    stoichiometry *= max(stoichiometries)

            self.cycle_stoichiometry.append({'cycle': cycle, 'stoichiometry': stoichiometry})

    @classmethod
    def make_canonical(cls, reactants):
        rle = collections.Counter(sorted(reactants))  # Cannot rely on dict(a) == dict(b) if items in a == items in b, but in different order, so sort them first.
        rle_reactants = ["{}{}".format(count, reactant) for reactant, count in rle.iteritems()]
        return "+".join(rle_reactants)

