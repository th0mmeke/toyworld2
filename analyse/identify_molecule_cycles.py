import networkx as nx
from collections import Counter
from collections import defaultdict


from identify_species_cycles import IdentifySpeciesCycles


class IdentifyMoleculeCycles(IdentifySpeciesCycles):

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
        self.reactants = set()  # set of smiles
        self.reactant_molecules = set()  # set of molecule ids
        self.smiles = defaultdict(set)  # smiles[x] = [all ids with smiles of x]

        print("Reactions=", len(reactions))
        for reaction in reactions:

            self.reactants.update(reaction['reactants'].values())
            self.reactant_molecules.update(reaction['reactants'].keys())

            canonical_reactants = self.make_canonical(reaction['reactants'].keys()) + '>'
            canonical_products = '>' + self.make_canonical(reaction['products'].keys())

            if not self.g.has_edge(canonical_reactants, canonical_products):
                self.g.add_edge(canonical_reactants, canonical_products)

            reactant_stoichiometry = Counter(reaction['reactants'].values())
            product_stoichiometry = Counter(reaction['products'].values())

            for id, smiles in reaction['reactants'].iteritems():
                if not self.g.has_node(id):  # reactant may already exist in the network
                    self.g.add_node(id, smiles=smiles)
                self.smiles[smiles].add(id)
                self.g.add_edge(id, canonical_reactants, stoichiometry=reactant_stoichiometry[smiles])

            for id, smiles in reaction['products'].iteritems():
                self.g.add_node(id, smiles=smiles)  # products are always new
                self.smiles[smiles].add(id)
                if not self.g.has_edge(canonical_products, id):
                    self.g.add_edge(canonical_products, id, stoichiometry=product_stoichiometry[smiles])

    def sample_reactant_stoichiometry(self, molecule, minimum_stoichiometry=0, max_depth=15):
        """
        Find all unique cycles that include each node that has the SMILES of acs_seed.

        :param acs_seed:
        :param minimum_stoichiometry:
        :param max_depth:
        :return: [{'cycle':cycle, 'stoichiometry':stoichiometry}]
        """
        reactant_stoichiometry = []

        #  find all cycles - each will be unique as from unique seed
        cycles = list(self.find_cycles_from_seed(self.g, molecule, max_depth))

        for cycle in cycles:
            stoichiometry = 1
            for start_node, end_node in zip(cycle[0:-1], cycle[1:]):
                if self.is_reaction(start_node) and not self.is_reaction(end_node):  # end_node is product
                    stoichiometry *= float(self.g[start_node][end_node]['stoichiometry'])
                if self.is_reaction(end_node) and not self.is_reaction(start_node):  # start_node must be reactant
                    stoichiometry /= float(self.g[start_node][end_node]['stoichiometry'])

            if stoichiometry >= minimum_stoichiometry:
                #  first make sure don't already have this cycle to a different end molecule of the same species
                unique = True
                for candidate_cycle in reactant_stoichiometry:
                    x = candidate_cycle['cycle']
                    if self.g.node[x[0]]['smiles'] == self.g.node[cycle[0]]['smiles'] and x[1:-1] == cycle[1:-1] and self.g.node[x[-1]]['smiles'] == self.g.node[cycle[-1]]['smiles']:
                        unique = False
                        break
                if unique:
                    reactant_stoichiometry.append({'cycle': cycle, 'stoichiometry': stoichiometry})

        return reactant_stoichiometry

    @staticmethod
    def find_cycles_from_seed(network, source, max_depth=5):  # http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
        target = network.node[source]['smiles']
        stack = [(source, [source])]
        while stack:
            (vertex, path) = stack.pop()
            for next_node in set(network.predecessors(vertex)) - set(path):
                if not IdentifySpeciesCycles.is_reaction(next_node) and network.node[next_node]['smiles'] == target:
                    yield list(reversed(path + [next_node]))  # reverse order as path found in reverse direction
                else:
                    if len(path) < max_depth:
                        stack.append((next_node, path + [next_node]))

    def get_smiles(self, id):
        return self.g.node[id]['smiles']

    def make_canonical(self, reactants):
        """Simpler for reactions involving molecule ids rather than smiles as no need for stoichiometry as each id is by definition unique."""
        return "+".join(sorted(reactants))