import networkx as nx
import re
import collections
import string


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

    def write_subset_to_graphml(self, seed, cycles, path='subgraph.graphml'):
        """
        Write a subset of the graph containing all cycles that include the molecule seed as reactant or product, replacing full labels by single letter ones that are better suited
        to graph visualisation. Includes nodes for each molecule that is a reactant or a product in any of these cycles.
        :param seed: Label of node to mark as node 'a' in the subgraph
        :param cycles: [[Labels of graph nodes]]
        :param path: Output filepath
        :return:
        """

        re_exp = r'(?<=\d).*?(?=\+\d|>|$)'  # without $: bug >1[H]c1n[c][c]nc([H])[n-][c]n1 should produce [H]c1n[c][c]nc([H])[n-][c]n1, instead produces []

        # Find all nodes to be included in the subgraph - nodes in cycles plus reactant or product molecules
        nodes = set()
        for cycle in cycles:
            for node in cycle:
                nodes.add(node)
                if node[-1] == '>' or node[0] == '>':
                    nodes.update(re.findall(re_exp, node))  # nodes that appear in reactions only

        try:
            nodes.remove('')  # artifact of the re_exp used
        except KeyError:
            pass

        subgraph = self.g.subgraph(nodes)  # Warning: subgraph attributes will change in self.g; nodes/edges however are independent (nx.subgraph() documentation)

        if len(nodes) > len(string.ascii_letters):
            raise ValueError("Too many nodes for visualisation")
        else:
            mapping = {}
            reaction_mapping = {}
            letters = list(string.ascii_letters)

            # Label the seed node as 'a'
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

            nx.relabel_nodes(subgraph, mapping, copy=False)  # remap labels in-place

        nx.write_graphml(subgraph, path)

    @classmethod
    def make_canonical(cls, reactants):
        rle = collections.Counter(sorted(reactants))  # Cannot rely on dict(a) == dict(b) if items in a == items in b, but in different order, so sort them first.
        rle_reactants = ["{}{}".format(count, reactant) for reactant, count in rle.iteritems()]
        return "+".join(rle_reactants)