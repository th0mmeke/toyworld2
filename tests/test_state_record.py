import unittest
import json

from state_record import StateRecord

class TestStateRecord(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(TypeError):
            StateRecord(data=[])
        self.assertIsInstance(StateRecord(data={}), StateRecord)

    def test_data_to_json(self):
        sr = StateRecord(data={'a': 1, 'b': 2})
        self.assertIsInstance(json.loads(sr.__str__()), object)