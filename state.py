from my_json import MyJSON
import json
import logging


class State(object):

    def __init__(self, filename):
        """
        Throw IOError if cannot open file
        :param filename:
        """

        self.i = 0
        self.f = open(filename, 'w+')
        logging.info('Opened state file {}'.format(self.f.name))
        self.f.write('[\n')

    def add(self, state_entry):
        s = json.dumps(state_entry, cls=MyJSON)
        logging.info(s)
        self.f.write((',\n' if self.i != 0 else '') + '{}'.format(s))
        self.i += 1

    def __del__(self):
        logging.info('Closing state file {}'.format(self.f.name))
        self.f.write('\n]\n')
        self.f.close()
