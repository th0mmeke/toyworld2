from my_json import MyJSON
import json
import logging


class State(object):

    def __init__(self, filename, initial_population):
        """
        Throw IOError if cannot open file
        :param filename:
        """

        self.i = 0
        self.f = open(filename, 'w+')
        logging.info('Opened state file {}'.format(self.f.name))

        self.f.write('{"initial_population": ' + json.dumps(initial_population, cls=MyJSON) + ',\n"reactions": [\n')

    def add(self, state_entry):
        s = json.dumps(state_entry, cls=MyJSON)
        self.f.write((',\n' if self.i != 0 else '') + '{}'.format(s))
        self.i += 1

    def close(self, final_population):
        logging.info('Closing state file {}'.format(self.f.name))
        self.f.write('\n],\n"final_population":' + json.dumps(final_population, cls=MyJSON) + '}')
        self.f.close()
