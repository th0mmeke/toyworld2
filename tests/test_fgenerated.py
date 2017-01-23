import unittest
import networkx as nx

from identify_species_cycles import IdentifySpeciesCycles
import fgenerated

class TestFgenerated(unittest.TestCase):

    def test_get_reactions(self):
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        e = IdentifySpeciesCycles(reactions=reactions)

        self.assertItemsEqual(['>1c', '>1e', '>1a', '>2b', '>2a'], fgenerated.get_reactions_b(e.g))
        self.assertItemsEqual(['1a>', '1b>', '1c>', '1e>'], fgenerated.get_reactions_a(e.g))

    def test_closure(self):

        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        closure = fgenerated.compute_closure(e.g, ['a', 'b'])
        self.assertEqual(set(['a', 'b', 'c']), closure)

        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'4': 'b', '5': 'e'}, 'products': {'6': 'd'}},
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        closure = fgenerated.compute_closure(e.g, ['a', 'b'])
        self.assertEqual(set(['a', 'b', 'c']), closure)

        foodset = ['a', 'b']
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c', '5': 'b'}, 'products': {'6': 'd'}},
            {'reactants': {'11': 'a', '12': 'b'}, 'products': {'13': 'f'}},
            {'reactants': {'13': 'f', '15': 'b'}, 'products': {'16': 'g'}}
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        closure = fgenerated.compute_closure(e.g, foodset)
        self.assertEqual(set(['a', 'b', 'c', 'd', 'f', 'g']), closure)

    def test_get_f(self):

        # Based on Fig 2 Hordijk and Steel 2004
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c', '5': 'b'}, 'products': {'6': 'd'}}
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        f = fgenerated.get_fgenerated(e.g, ['a', 'b'])
        print("f", f.nodes())
        self.assertTrue(nx.is_isomorphic(e.g, f))  # small network has single RAF

        # Steel et al 2013 fig 2b - F-generated
        foodset = ['f1', 'f2', 'f3', 'f4', 'f5']
        reactions = [
            {'reactants': {'1': 'f1', '2': 'p3'}, 'products': {'3': 'p1', '4': 'p4'}},  # r1
            {'reactants': {'5': 'f2', '6': 'p1'}, 'products': {'7': 'p2'}},  # r2
            {'reactants': {'8': 'f3', '9': 'p2'}, 'products': {'10': 'p3'}},  # r3
            {'reactants': {'11': 'f4', '12': 'f5'}, 'products': {'13': 'p1'}},  # r4
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        expected_nodes = [node for node in e.g.nodes() if not IdentifySpeciesCycles.is_reaction(node)]
        raf = fgenerated.get_fgenerated(e.g, foodset)
        self.assertItemsEqual(expected_nodes, list(fgenerated.compute_closure(raf, foodset)))
        actual_nodes = [node for node in raf if not IdentifySpeciesCycles.is_reaction(node)]
        self.assertListEqual(expected_nodes, actual_nodes)

        # Steel et al 2013 fig 2a - not F-generated
        foodset = ['f1', 'f2']
        reactions = [
            {'reactants': {'1': 'f1', '2': 'p3'}, 'products': {'3': 'p1', '4': 'p4'}},  # r1
            {'reactants': {'5': 'f2', '6': 'p1'}, 'products': {'7': 'p2'}},  # r2
            {'reactants': {'8': 'f3', '9': 'p2'}, 'products': {'10': 'p3'}},  # r3
        ]
        e = IdentifySpeciesCycles(reactions=reactions)

        f = fgenerated.get_fgenerated(e.g, foodset)
        self.assertItemsEqual(['f1', 'f2'], f.nodes())

    def test_get_irr_f(self):
        """
        Based on Fig 2 Hordijk and Steel 2004
        """
        foodset = ['a', 'b']
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c', '5': 'b'}, 'products': {'6': 'd'}}
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        molecules = fgenerated.get_irr_fgenerated(e.g, foodset)['molecules']
        expected_nodes = ['c']

        self.assertTrue(molecules == ['c'] or molecules == ['d'])  # small network has single RAF and hence irrRAF

    def test_reduction_to_irr_f(self):
        """
        Based on fig 2 from Hordijk and Steel 2004
        """
        foodset = ['a', 'b']

        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c', '5': 'b'}, 'products': {'6': 'd'}},
            {'reactants': {'11': 'a', '12': 'b'}, 'products': {'13': 'f'}},
            {'reactants': {'13': 'f', '15': 'b'}, 'products': {'16': 'g'}}
        ]

        e = IdentifySpeciesCycles(reactions=reactions)
        for i in range(0, 40):
            irr = fgenerated.get_irr_fgenerated(e.g, foodset)['molecules']
            print(irr)
            #self.assertEqual(1, len(irr))
            #self.assertTrue(irr == ['c'] or irr == ['f'])
