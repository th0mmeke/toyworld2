from json import JSONEncoder
from atom import Atom
from molecule import Molecule


class MyJSON(JSONEncoder):

    def default(self, o):

        if isinstance(o, (Atom, Molecule)):
            return str(o)
        JSONEncoder.default(self, o)
