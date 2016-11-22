import networkx as nx
from collections import Counter
from collections import defaultdict


from evaluator_cycles import EvaluatorCycles


class EvaluatorActualCycles(EvaluatorCycles):

    def __init__(self, reactions, minimum_length=0):

        """
        Identify all reaction cycles from the reaction network graph. A cycle is a sequential set of reactions from a molecule
        returning back to a molecule of the same type, and where at least one product molecule (not just type) from one step is a reactant molecule of the
        next. This differs from a potential cycle where the requirement is simply that products and reactants are of the same *type*.

        Note that the final step in a cycle must be to return to a molecule of the same type, not the same molecule, as we started from. That earlier molecule
        won't still be around...

        :param reactions:
        """

        self.g = nx.DiGraph()
        self.reactants = set()
        self.smiles = defaultdict(list)  # smiles[x] = [all ids with smiles of x]

        for reaction in reactions:

            self.reactants.update(reaction['reactants'].values())

            canonical_reactants = EvaluatorCycles.make_canonical(reaction['reactants'].keys()) + '>'
            canonical_products = '>' + EvaluatorCycles.make_canonical(reaction['products'].keys())

            if not self.g.has_edge(canonical_reactants, canonical_products):
                self.g.add_edge(canonical_reactants, canonical_products)

            reactant_stoichiometry = Counter(reaction['reactants'].values())
            product_stoichiometry = Counter(reaction['products'].values())

            for id, smiles in reaction['reactants'].iteritems():
                if not self.g.has_node(id):  # reactant may already exist in the network
                    self.g.add_node(id, smiles=smiles)
                self.smiles[smiles].append(id)
                self.g.add_edge(id, canonical_reactants, stoichiometry=reactant_stoichiometry[smiles])

            for id, smiles in reaction['products'].iteritems():
                self.g.add_node(id, smiles=smiles)  # products are always new
                self.smiles[smiles].append(id)
                if not self.g.has_edge(canonical_products, id):
                    self.g.add_edge(canonical_products, id, stoichiometry=product_stoichiometry[smiles])

    def get_reactant_stoichiometry(self, acs_seed, minimum_stoichiometry=0, max_depth=5):
        """
        Find the shortest cycle from each node that has the SMILES of acs_seed, if the node is in a cycle.

        :param acs_seed:
        :param minimum_stoichiometry:
        :param max_depth:
        :return:
        """
        reactant_stoichiometry = []

        for id_seed in self.smiles[acs_seed]:

            cycles = list(self.find_shortest_paths(self.g, id_seed, acs_seed, max_depth))

            # eliminate duplicate cycles
            canonical_cycles = set([self._make_canonical(cycle) for cycle in cycles])

            for cycle in canonical_cycles:
                stoichiometry = 1
                for start_node, end_node in zip(cycle[0:-1], cycle[1:]):
                    if self.is_reaction(start_node) and not self.is_reaction(end_node):  # end_node is product
                        stoichiometry *= float(self.g[start_node][end_node]['stoichiometry'])
                    if self.is_reaction(end_node) and not self.is_reaction(start_node):  # start_node must be reactant
                        stoichiometry /= float(self.g[start_node][end_node]['stoichiometry'])

                if stoichiometry >= minimum_stoichiometry:
                    reactant_stoichiometry.append({'cycle': cycle, 'stoichiometry': stoichiometry})

        return reactant_stoichiometry

    @staticmethod
    def find_shortest_paths(network, source, target, max_depth=5):  # http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
        stack = [(source, [source])]
        while stack:
            (vertex, path) = stack.pop()
            for next_node in set(network.predecessors(vertex)) - set(path):
                if not EvaluatorCycles.is_reaction(next_node) and network.node[next_node]['smiles'] == target:
                    yield list(reversed(path + [next_node]))  # reverse order as path found in reverse direction
                else:
                    if len(path) < max_depth:
                        stack.append((next_node, path + [next_node]))

    @classmethod
    def _make_canonical(cls, path):
        hashes = [hash(x) for x in path]
        start_idx = hashes.index(min(hashes))
        return tuple(path[(i + start_idx) % len(path)] for i in range(len(path)))
