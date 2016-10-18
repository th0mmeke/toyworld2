import random
import unittest

import uniform_product_selection


class TestUniformProductSelection(unittest.TestCase):

    def testProductSelection(self):
        self.assertIsNone(uniform_product_selection.product_selection(None))