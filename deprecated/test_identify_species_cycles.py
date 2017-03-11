import unittest

from deprecated.identify_species_cycles import IdentifySpeciesCycles


class TestIdentifySpeciesCycles(unittest.TestCase):

    def testSimple(self):
        # a->a->b->a => a->b
        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'a'}}
                ]
        r = IdentifySpeciesCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)
        self.assertEqual(2, len(r))

        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
                {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
                ]
        r = IdentifySpeciesCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)
        self.assertEqual(3, len(r))

        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        self.assertEqual(5, len(IdentifySpeciesCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testMultiple(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = IdentifySpeciesCycles(reactions=reactions)
        r = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(2, len(r))

    def testStoichiometry(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c'}, 'products': {'4': 'a', '5': 'a', '6': 'd'}}
        ]

        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        for cycle in cycles:
            self.assertEqual(2, cycle['stoichiometry'])

        reactions = [
            {'reactants': {'1': 'a', '2': 'a'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c'}, 'products': {'4': 'a', '5': 'a', '6': 'd'}}
        ]

        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        for cycle in cycles:
            self.assertEqual(1, cycle['stoichiometry'])

    def test_get_reactants(self):
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual({'a', 'c', 'e'}, e.get_reactant_set(cycles[1]['cycle']))

        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual({'a', 'c', 'f', 'b'}, e.get_reactant_set(cycles[1]['cycle']))

    def test_get_products(self):
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual({'a', 'b', 'e'}, e.get_product_set(cycles[1]['cycle']))

        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual({'a', 'c', 'f', 'd'}, e.get_product_set(cycles[1]['cycle']))

    def test_get_foodset(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]
        e = IdentifySpeciesCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual({'b'}, e.get_foodset(cycles[1]['cycle']))

        # Based on Fig 2 Hordijk and Steel 2004
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c', '5': 'b'}, 'products': {'6': 'd'}}
        ]
        e = IdentifySpeciesCycles(reactions=reactions)

        foodset = e.get_foodset({'c'})
        self.assertItemsEqual({'a', 'b'}, foodset)

        foodset = e.get_foodset({'d'})
        self.assertItemsEqual({'a', 'b', 'c'}, foodset)

