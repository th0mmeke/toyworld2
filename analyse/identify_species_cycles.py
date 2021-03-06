import networkx as nx
import re
from collections import Counter
import string
import cycle_utilities


class IdentifySpeciesCycles(object):

    def __init__(self, reactions):
        """
        Construct initial reaction graph.
        :param reactions:
        """

        # Build graph

        self.g = nx.DiGraph()
        self.reactants = set()

        for reaction in reactions:

            # Distinct treatment of reactant and product sets to avoid edges being added from reactants back to molecules
            canonical_reactants = self.make_canonical(reaction['reactants'].values()) + '>'
            canonical_products = '>' + self.make_canonical(reaction['products'].values())
            self.reactants.update(reaction['reactants'].values())

            if not self.g.has_edge(canonical_reactants, canonical_products):
                self.g.add_edge(canonical_reactants, canonical_products)

            for reactant, count in Counter(reaction['reactants'].values()).iteritems():
                # Only add if not already present
                if not self.g.has_edge(reactant, canonical_reactants):
                    self.g.add_edge(reactant, canonical_reactants, stoichiometry=count)

            for product, count in Counter(reaction['products'].values()).iteritems():
                if not self.g.has_edge(canonical_products, product):
                    self.g.add_edge(canonical_products, product, stoichiometry=count)

    def get_population_stoichiometry(self, minimum_stoichiometry=1, max_depth=5):
        """
        Return dictionary of stoichiometry for each unique cycle in the reaction set
        :param minimum_length:
        :param minimum_stoichiometry:
        :param max_depth:
        :return:
        """

        population_stoichiometry = []
        for reactant in self.reactants:
            population_stoichiometry.extend(self.get_reactant_stoichiometry(reactant, minimum_stoichiometry, max_depth))

        return population_stoichiometry

    def get_reactant_stoichiometry(self, acs_seed, minimum_stoichiometry=0, max_depth=5):

        reactant_stoichiometry = []

        cycles = list(nx.all_simple_paths(self.g, acs_seed, acs_seed, cutoff=max_depth))

        for cycle in cycles:
            stoichiometry = 1
            for start_node, end_node in zip(cycle[0:-1], cycle[1:]):
                if cycle_utilities.is_reaction(start_node) and not cycle_utilities.is_reaction(end_node):  # end_node is product
                    stoichiometry *= float(self.g[start_node][end_node]['stoichiometry'])
                if cycle_utilities.is_reaction(end_node) and not cycle_utilities.is_reaction(start_node):  # start_node must be reactant
                    stoichiometry /= float(self.g[start_node][end_node]['stoichiometry'])

            if stoichiometry >= minimum_stoichiometry:
                reactant_stoichiometry.append({'cycle': cycle, 'stoichiometry': stoichiometry})

        return reactant_stoichiometry

    # def get_foodset(self,  cycle):
    #
    #     """
    #     Food set is all the reactants required by the cycle less those it can generate itself.
    #     :param cycle:
    #     :return:
    #     """
    #
    #     cycle = set(cycle)
    #     foodset = cycle
    #     new_foodset = foodset
    #
    #     while new_foodset != set():
    #         new_foodset = self.get_reactant_set(foodset).difference(self.get_product_set(foodset))
    #         foodset = foodset.union(new_foodset)
    #     return foodset - cycle
    #
    # def get_reactant_set(self, cycle):
    #     """
    #     Return the complete set of reactants required for this cycle - that is, incoming reactants plus the reactants generated by the cycle
    #     :param cycle:
    #     :return:
    #     """
    #     reactants = set([node for node in cycle if not cycle_utilities.is_reaction(node)])
    #     node_front = reactants
    #     for i in range(0, 3):
    #         new_front = []
    #         for node in node_front:
    #             predecessors = self.g.predecessors(node)
    #             for predecessor in predecessors:
    #                 if cycle_utilities.is_reaction(predecessor):
    #                     new_front.append(predecessor)
    #                 else:
    #                     reactants.add(predecessor)
    #         node_front = new_front
    #     return reactants
    #
    # def get_product_set(self, cycle):
    #     """
    #     Return all possible products from this cycle. A product is anything which appears on the right-hand side of a reaction in the cycle.
    #     :param cycle:
    #     :return:
    #     """
    #     products = set([node for node in cycle if not cycle_utilities.is_reaction(node)])
    #     node_front = products
    #     for i in range(0, 3):
    #         new_front = []
    #         for node in node_front:
    #             successors = self.g.successors(node)
    #             for successor in successors:
    #                 if cycle_utilities.is_reaction(successor):
    #                     new_front.append(successor)
    #                 else:
    #                     products.add(successor)
    #         node_front = new_front
    #     return products

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
                if cycle_utilities.is_reaction(node):
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
                if cycle_utilities.is_reaction(node):
                    new_label = letters.pop(0)
                    mapping[node] = new_label

            for node in nodes:
                if cycle_utilities.is_reaction(node):
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

    def make_canonical(self, reactants):
        rle = Counter(sorted(reactants))  # Cannot rely on dict(a) == dict(b) if items in a == items in b, but in different order, so sort them first.
        rle_reactants = ["{}{}".format(count, reactant) for reactant, count in rle.iteritems()]
        return "+".join(rle_reactants)
