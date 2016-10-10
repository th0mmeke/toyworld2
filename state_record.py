import json


class StateRecord(object):

    def __init__(self, data):
        if not isinstance(data, dict):
            raise TypeError
        self.data = data

    def __str__(self):
        return json.dumps(self.data)

