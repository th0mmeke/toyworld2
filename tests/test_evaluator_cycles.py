import unittest

from evaluator_cycles import EvaluatorCycles


class TestEvaluatorCycles(unittest.TestCase):

    def testSimple(self):
        # a->a->b->a => a->b
        reactions = [
                {'reactants': ['a'], 'products': ['b']},
                {'reactants': ['b'], 'products': ['a']}
                ]
        self.assertEqual(2, len(EvaluatorCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
                {'reactants': ['a'], 'products': ['b', 'b']},
                {'reactants': ['b'], 'products': ['c']},
                {'reactants': ['c'], 'products': ['a']},
                ]
        self.assertEqual(3, len(EvaluatorCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

        reactions = [
            {'reactants': ['a'], 'products': ['b', 'b']},
            {'reactants': ['b'], 'products': ['c']},
            {'reactants': ['c'], 'products': ['a']},
            {'reactants': ['a'], 'products': ['e']},
            {'reactants': ['e'], 'products': ['a', 'a']},
        ]
        self.assertEqual(5, len(EvaluatorCycles(reactions=reactions).get_population_stoichiometry(max_depth=10)))

    def testMultiple(self):
        # aba -> ccf -> cf -> ad has stoichiometry of 2, other two reaction cycles have stoichiometry of 1.
        reactions = [
            {'reactants': ['a', 'b', 'a'], 'products': ['c', 'd']},
            {'reactants': ['b', 'a', 'a'], 'products': ['c', 'c', 'f']},
            {'reactants': ['c', 'f'], 'products': ['a', 'd']}
        ]

        e = EvaluatorCycles(reactions=reactions)
        self.assertEqual(6, len(e.get_population_stoichiometry(max_depth=10)))

    def testStoichiometry(self):
        reactions = [
            {'reactants': ['a', 'b'], 'products': ['c']},
            {'reactants': ['c'], 'products': ['a', 'a', 'd']}
        ]

        e = EvaluatorCycles(reactions=reactions)
        cycles = e.get_population_stoichiometry(max_depth=10)
        for cycle in cycles:
            self.assertEqual(2, cycle['stoichiometry'])

