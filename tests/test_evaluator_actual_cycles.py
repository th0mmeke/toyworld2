import unittest

from evaluator_actual_cycles import EvaluatorActualCycles


class TestEvaluatorActualCycles(unittest.TestCase):

    def testSimple(self):
        # a->a->b->a => a->b
        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b'}},
                {'reactants': {'2': 'b'}, 'products': {'4': 'a'}}
                ]
        self.assertEqual(1, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testSimple2(self):
        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
                {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
                ]
        self.assertEqual(1, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
                {'reactants': {'4': 'c'}, 'products': {'5': 'd'}},
                {'reactants': {'5': 'd'}, 'products': {'6': 'c'}},
                ]
        self.assertEqual(1, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testSimple3(self):
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
            {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
            {'reactants': {'4': 'c'}, 'products': {'5': 'a'}},
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        self.assertEqual(2, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testLoop(self):
        # Two loops - through 7 and through 8? Or just one? Just one - don't double count
        reactions = [
            {'reactants': {'5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        self.assertEqual(1, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
            {'reactants': {'4': 'a', '5': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a', '8': 'a'}},
        ]
        self.assertEqual(1, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
            {'reactants': {'4': 'a'}, 'products': {'6': 'e'}},
            {'reactants': {'9': 'a'}, 'products': {'10': 'b'}},
            {'reactants': {'6': 'e'}, 'products': {'7': 'a'}},
            {'reactants': {'10': 'b'}, 'products': {'11': 'a'}},
        ]
        self.assertEqual(2, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testBrokenCycles(self):
        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b', '3': 'b'}},
                {'reactants': {'3': 'b'}, 'products': {'4': 'c'}},
                {'reactants': {'6': 'c'}, 'products': {'5': 'a'}},  # Not the same 'c'
                ]
        self.assertEqual(0, len(EvaluatorActualCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testMultiple(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b', '3': 'a'}, 'products': {'4': 'c', '5': 'd'}},
            {'reactants': {'6': 'b', '7': 'a', '8': 'a'}, 'products': {'9': 'c', '10': 'c', '11': 'f'}},
            {'reactants': {'9': 'c', '11': 'f'}, 'products': {'12': 'a', '13': 'd'}}
        ]

        e = EvaluatorActualCycles(reactions=reactions)
        self.assertEqual(1, len(e.get_population_stoichiometry(max_depth=10)))  # same path

    def testStoichiometry(self):
        reactions = [
            {'reactants': {'1': 'a', '2': 'b'}, 'products': {'3': 'c'}},
            {'reactants': {'3': 'c'}, 'products': {'4': 'a', '5': 'a', '6': 'd'}}
        ]

        e = EvaluatorActualCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        for cycle in cycles:
            self.assertEqual(2, cycle['stoichiometry'])

    def testTwoCycles(self):

        # Simple case - one cycle, then another, sharing same smiles
        reactions = [
                {'reactants': {'1': 'a'}, 'products': {'2': 'b'}},
                {'reactants': {'2': 'b'}, 'products': {'3': 'a'}},
                {'reactants': {'4': 'a'}, 'products': {'5': 'b'}},
                {'reactants': {'5': 'b'}, 'products': {'6': 'a'}}
                ]
        e = EvaluatorActualCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(2, len(cycles))

        # Order of reactions shouldn't matter...
        reactions = [
            {'reactants': {'1': 'a'}, 'products': {'2': 'b'}},
            {'reactants': {'4': 'a'}, 'products': {'5': 'b'}},
            {'reactants': {'2': 'b'}, 'products': {'3': 'a'}},
            {'reactants': {'5': 'b'}, 'products': {'6': 'a'}}
        ]
        e = EvaluatorActualCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        self.assertEqual(2, len(cycles))


