import networkx as nx
from collections import OrderedDict
import itertools

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

        reaction_network = nx.DiGraph()
        cycles = {}
        id_cycles = {}
        active_molecules = {}
        cycle_id = 0

        for reaction in reactions:

            reactants = [{'id': id, 'smiles': smiles} for id, smiles in reaction['reactants'].iteritems()]
            products = [{'id': id, 'smiles': smiles} for id, smiles in reaction['products'].iteritems()]
            
            for reactant, product in itertools.product(reactants, products):
                
                # First, check if this combination extends an existing cycle
                try:
                    candidate_cycles = active_molecules[reactant['id']]
                except KeyError:
                    pass
                else:
                    # yes it does, and candidate_cycles contains pointers to all of them
                    for cycle_id in candidate_cycles:
                        expected_smiles, count = cycles[cycle_id].popitem(last=False)
                        active_molecules[reactant['id']].remove(cycle_id)
                        if product['smiles'] == expected_smiles:

                            # cycle continues - so update count for this mol, and advance to next position in the cycle
                            cycles[cycle_id][expected_smiles] = count + 1
                            try:
                                active_molecules[product['id']].append(cycle_id)
                            except KeyError:
                                active_molecules[product['id']] = [cycle_id]
                        else:
                            cycles[cycle_id][expected_smiles] = count  # just restore it if cycle broken to preserve counts

                # Second, check if we have a new cycle
                if reaction_network.has_node(reactant['id']):  # first iteration won't have reactants added yet

                    # check for all paths from reactant['id'] to any node in the network with same smiles as product (product['smiles'])
                    paths = self.find_shortest_paths(reaction_network, reactant['id'], product['smiles'])
                    for id_path in paths:
                        if len(id_path) >= minimum_length:
                            new_cycle = OrderedDict()
                            for x in id_path:
                                new_cycle[reaction_network.node[x]['smiles']] = 1
                            start_smiles, count = new_cycle.popitem(last=False)
                            new_cycle[start_smiles] = count  # move the current location to last in list; first item is now the next one we'll expect
                            # add the cycle to the list of active cycles
                            cycle_id += 1
                            cycles[cycle_id] = new_cycle  # cycles recorded by smiles
                            id_cycles[cycle_id] = id_path  # cycles recorded by id
                            # now add this cycle to the list of cycles of which this molecule is a participant
                            try:
                                active_molecules[product['id']].append(cycle_id)
                            except KeyError:
                                active_molecules[product['id']] = [cycle_id]

                reaction_network.add_edge(reactant['id'], product['id'])

                reaction_network.node[reactant['id']]['smiles'] = reactant['smiles']
                reaction_network.node[product['id']]['smiles'] = product['smiles']

        # Now aggregate by canonical form
        self.compliant_cycles = {}
        for cycle_id in cycles.keys():
            cycle = EvaluatorActualCycles.make_canonical(tuple(cycles[cycle_id].keys()))
            try:
                self.compliant_cycles[cycle] += 1
            except KeyError:
                self.compliant_cycles[cycle] = 1

    def get_population_stoichiometry(self, minimum_stoichiometry=1):
        """
        Return dictionary of stoichiometry for each unique cycle in the reaction set
        :param minimum_length:
        :param minimum_stoichiometry:
        :param max_depth:
        :return:
        """

        return [{'cycle': cycle, 'stoichiometry': stoichiometry} for cycle, stoichiometry in self.compliant_cycles.iteritems() if stoichiometry >= minimum_stoichiometry]

    @classmethod
    def make_canonical(cls, path):
        hashes = [hash(x) for x in path]
        start_idx = hashes.index(min(hashes))
        return tuple(path[(i + start_idx) % len(path)] for i in range(len(path)))

    @staticmethod
    def find_shortest_paths(network, source, target):  # http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
        stack = [(source, [source])]
        while stack:
            (vertex, path) = stack.pop()
            for next_node in set(network.predecessors(vertex)) - set(path):
                if network.node[next_node]['smiles'] == target:
                    yield list(reversed(path + [next_node]))  # reverse order as path found in reverse direction
                else:
                    stack.append((next_node, path + [next_node]))