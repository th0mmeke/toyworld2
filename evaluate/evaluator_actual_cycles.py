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
        self.smiles = defaultdict(set)  # smiles[x] = [all ids with smiles of x]

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
                self.smiles[smiles].add(id)
                self.g.add_edge(id, canonical_reactants, stoichiometry=reactant_stoichiometry[smiles])

            for id, smiles in reaction['products'].iteritems():
                self.g.add_node(id, smiles=smiles)  # products are always new
                self.smiles[smiles].add(id)
                if not self.g.has_edge(canonical_products, id):
                    self.g.add_edge(canonical_products, id, stoichiometry=product_stoichiometry[smiles])

    def get_population_stoichiometry(self, minimum_stoichiometry=0, max_depth=5):
        """
        Find the shortest cycle from each node that has the SMILES of acs_seed, if the node is in a cycle.

        :param minimum_stoichiometry:
        :param max_depth:
        :return: [{'cycle':cycle, 'stoichiometry':stoichiometry}]
        """
        reactant_stoichiometry = []

        #  find all cycles - each will be unique as from unique seed
        cycles = list(self.simple_cycles())
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

    def get_smiles(self, id):
        return self.g.node[id]['smiles']

    def simple_cycles(self):
        """Reimplementation of nx.simple_cycles with check for same smiles to close cycle rather than same node.
        Source taken from https://github.com/networkx/networkx/blob/master/networkx/algorithms/cycles.py
        """

        def _unblock(thisnode, blocked, B):
            stack = set([thisnode])
            while stack:
                node = stack.pop()
                if node in blocked:
                    blocked.remove(node)
                    stack.update(B[node])
                    B[node].clear()

        # Johnson's algorithm requires some ordering of the nodes.
        # We assign the arbitrary ordering given by the strongly connected comps
        # There is no need to track the ordering as each node removed as processed.
        subG = type(self.g)(self.g.edges())  # save the actual graph so we can mutate it here
        # We only take the edges because we do not want to
        # copy edge and node attributes here.
        sccs = list(nx.strongly_connected_components(subG))
        while sccs:
            scc = sccs.pop()
            # order of scc determines ordering of nodes
            startnode = scc.pop()
            # Processing node runs "circuit" routine from recursive version
            path = [startnode]
            blocked = set()  # vertex: blocked from search?
            closed = set()  # nodes involved in a cycle
            blocked.add(startnode)
            B = defaultdict(set)  # graph portions that yield no elementary circuit
            stack = [(startnode, list(subG[startnode]))]  # subG gives component nbrs
            while stack:
                thisnode, nbrs = stack[-1]
                if nbrs:
                    nextnode = nbrs.pop()
                    if 'smiles' in self.g.node[nextnode].keys() and 'smiles' in self.g.node[startnode].keys() and self.g.node[nextnode]['smiles'] == self.g.node[startnode]['smiles']:
                        yield path[:] + [nextnode]
                        closed.update(path)
                    elif nextnode not in blocked:
                        path.append(nextnode)
                        stack.append((nextnode, list(subG[nextnode])))
                        closed.discard(nextnode)
                        blocked.add(nextnode)
                        continue
                # done with nextnode... look for more neighbors
                if not nbrs:  # no more nbrs
                    if thisnode in closed:
                        _unblock(thisnode, blocked, B)
                    else:
                        for nbr in subG[thisnode]:
                            if thisnode not in B[nbr]:
                                B[nbr].add(thisnode)
                    stack.pop()
                    path.pop()
            # done processing this node
            subG.remove_node(startnode)
            H = subG.subgraph(scc)  # make smaller to avoid work in SCC routine

            sccs.extend(list(nx.strongly_connected_components(H)))