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

            if not self.g.has_edge(canonical_reactants, canonical_products):
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

    def get_population_stoichiometry(self, minimum_length=0, minimum_stoichiometry=1, max_depth=5):
        """
        Return dictionary of stoichiometry for each unique cycle in the reaction set
        :param minimum_length:
        :param minimum_stoichiometry:
        :param max_depth:
        :return:
        """

        population_stoichiometry = []
        for reactant in self.reactants:
            if len(reactant) >= minimum_length:
                population_stoichiometry.extend(self.get_reactant_stoichiometry(reactant, minimum_stoichiometry, max_depth))

        return population_stoichiometry

    def get_reactant_stoichiometry(self, acs_seed, minimum_stoichiometry=1, max_depth=5):

        # print("Seed: {}".format(acs_seed))
        reactant_stoichiometry = []

        cycles = list(nx.all_simple_paths(self.g, acs_seed, acs_seed, cutoff=max_depth))

        for cycle in cycles:
            stoichiometry = 1
            for start_node, end_node in zip(cycle[0:-1], cycle[1:]):

                # get the max stoichiometry - as looking for ACS of high stoichiometry, take the highest
                stoichiometries = [edge['stoichiometry'] for edge in self.g.edge[start_node][end_node].itervalues() if 'stoichiometry' in edge.keys()]
                if len(stoichiometries) > 0:
                    stoichiometry *= max(stoichiometries)

            if stoichiometry >= minimum_stoichiometry:
                reactant_stoichiometry.append({'cycle': cycle, 'stoichiometry': stoichiometry})
                # print({'cycle': cycle, 'stoichiometry': stoichiometry})

        return reactant_stoichiometry

    @classmethod
    def make_canonical(cls, reactants):
        rle = collections.Counter(sorted(reactants))  # Cannot rely on dict(a) == dict(b) if items in a == items in b, but in different order, so sort them first.
        rle_reactants = ["{}{}".format(count, reactant) for reactant, count in rle.iteritems()]
        return "+".join(rle_reactants)