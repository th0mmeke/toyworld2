import unittest
from my_json import MyJSON

import json

from state import State
from atom import Atom


class TestState(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(State(persistence=lambda x: x), State)

    def test_add(self):
        expected = {'a': 1}
        actual = State(persistence=lambda x: json.dumps(x, cls=MyJSON)).add(expected)
        self.assertEqual(json.loads(actual), expected)

        expected = [Atom('A')]
        actual = State(persistence=lambda x: json.dumps(x, cls=MyJSON)).add(expected)
        self.assertEqual(json.loads(actual), [str(Atom('A'))])
        
        # Note that {Atom('A'): 1} will fail encoding - limitation of encoder that cannot use object as key
