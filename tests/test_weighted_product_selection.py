import random
import unittest

import weighted_product_selection


class TestWeightedProductSelection(unittest.TestCase):

    def testChoice(self):
        with self.assertRaises(ValueError):
            weighted_product_selection.choice([1], [2, 3])
        self.assertEqual(1, weighted_product_selection.choice([1], weights=[1]))
        self.assertEqual(None, weighted_product_selection.choice([1], weights=[0]))

        random.seed(123456890)
        self.assertEqual(1, weighted_product_selection.choice([1, 2], weights=[0.5, 0.5]))
        n = 10000
        actual = [weighted_product_selection.choice([1, 0], weights=[2, 4]) for i in range(n)]
        self.assertAlmostEqual(actual.count(0), n*4/6, delta=700)
        self.assertAlmostEqual(actual.count(1), n*2/6, delta=400)

    def testProductSelection(self):
        self.assertIsNone(weighted_product_selection.product_selection(None))