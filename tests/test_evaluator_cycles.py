import unittest

from evaluator_cycles import EvaluatorCycles


class TestEvaluatorCycles(unittest.TestCase):

    def testSimple(self):
        # a->a->b->a => a->b
        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'a'}}
                ]
        self.assertEqual(2, len(EvaluatorCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
                {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
                ]
        self.assertEqual(3, len(EvaluatorCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        self.assertEqual(5, len(EvaluatorCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testMultiple(self):
        # aba -> ccf -> cf -> ad has stoichiometry of 2, other two reaction cycles have stoichiometry of 1.
        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = EvaluatorCycles(reactions=reactions)
        self.assertEqual(6, len(e.get_population_stoichiometry(max_depth=10)))

    def testStoichiometry(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c'}, 'products': {'4': 'a', '5': 'a', '6': 'd'}}
        ]

        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        for cycle in cycles:
            self.assertEqual(2, cycle['stoichiometry'])

    def test_get_reactants(self):
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(set({'a', 'c', 'e'}), e.get_reactant_set(cycles[1]['cycle']))

        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(set({'a', 'c', 'f', 'b'}), e.get_reactant_set(cycles[1]['cycle']))

    def test_get_products(self):
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(set({'a', 'b', 'e'}), e.get_product_set(cycles[1]['cycle']))

        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(set({'a', 'c', 'f', 'd'}), e.get_product_set(cycles[1]['cycle']))

    def test_get_food_set(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]
        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(set({'b'}), e.get_food_set(cycles[1]['cycle']))
