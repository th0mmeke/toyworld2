import networkx as nx
import itertools
import collections


class EvaluatorCycles(object):

    def __init__(self, reactions):
        """
        Construct initial reaction graph.
        :param reactions:
        """

        # Build graph

        self.g = nx.MultiDiGraph()
        self.reactants = set()

        for reaction in reactions:

            # Distinct treatment of reactant and product sets to avoid edges being added from reactants back to molecules
            canonical_reactants = EvaluatorCycles.make_canonical(reaction['reactants']) + '>'
            canonical_products = '>' + EvaluatorCycles.make_canonical(reaction['products'])
            self.reactants.update(reaction['reactants'])

            self.g.add_edge(canonical_reactants, canonical_products)

            for reactant in reaction['reactants']:
                # Only add if not already present
                if not self.g.has_edge(reactant, canonical_reactants):
                    self.g.add_edge(reactant, canonical_reactants)
            for product, count in collections.Counter(reaction['products']).iteritems():
                add_edge = True
                if self.g.has_edge(canonical_products, product):
                    existing_stoichiometries = set([edge['stoichiometry'] for edge in self.g.edge[canonical_products][product].itervalues()])
                    if count in existing_stoichiometries:
                        add_edge = False
                if add_edge:
                    self.g.add_edge(canonical_products, product, stoichiometry=count)

        print(self.g.number_of_nodes(), self.g.number_of_edges())

    def get_cycle_stoichiometry(self, minimum_length=0, minimum_stoichiometry=1, max_depth=5):
        """
        Return dictionary of stoichiometry for each unique cycle in the reaction set
        :param minimum_length:
        :param minimum_stoichiometry:
        :param max_depth:
        :return:
        """

        cycle_stoichiometry = []
        for candidate_acs_seed in self.reactants:

            if len(candidate_acs_seed) >= minimum_length:
                cycles = list(nx.all_simple_paths(self.g, candidate_acs_seed, candidate_acs_seed, cutoff=max_depth))
                for cycle in cycles:
                    stoichiometry = 1
                    for start_node, end_node in zip(cycle[0:-1], cycle[1:]):
                        # get the max stoichiometry - looking for ACS of high stoichiometry, assume take the highest reaction pathway
                        stoichiometries = [edge['stoichiometry'] for edge in self.g.edge[start_node][end_node].itervalues() if 'stoichiometry' in edge.keys()]
                        if len(stoichiometries) > 0:
                            stoichiometry *= max(stoichiometries)

                    if stoichiometry >= minimum_stoichiometry:
                        cycle_stoichiometry.append({'cycle': cycle, 'stoichiometry': stoichiometry})
                        print({'cycle': cycle, 'stoichiometry': stoichiometry})
        return cycle_stoichiometry

    @classmethod
    def make_canonical(cls, reactants):
        rle = collections.Counter(sorted(reactants))  # Cannot rely on dict(a) == dict(b) if items in a == items in b, but in different order, so sort them first.
        rle_reactants = ["{}{}".format(count, reactant) for reactant, count in rle.iteritems()]
        return "+".join(rle_reactants)

[u'[H][N-]c1nc2c(nc([H])n2[H])c(=O)n1[H]', '1[H][N-]c1nc2c(nc([H])n2[H])c(=O)n1[H]+1[H]c1nc2c(=O)n([H])c(N([H])[H])nc2n1[H]>', '>1[H+]+2[H][N-]c1nc2c(nc([H])n2[H])c(=O)n1[H]', u'[H][N-]c1nc2c(nc([H])n2[H])c(=O)n1[H]']}
[u'[H]c1nc2c(=O)n([H])[c]nc2[n-]1', '1[H]c1nc2c(=O)n([H])[c]nc2[n-]1+1[H]c1nc2c(=O)n([H])[c]nc2n1[H]>', '>1[H+]+2[H]c1nc2c(=O)n([H])[c]nc2[n-]1', u'[H]c1nc2c(=O)n([H])[c]nc2[n-]1']}