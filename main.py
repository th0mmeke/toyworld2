import random
import json
from my_json import MyJSON

from toyworld import ToyWorld
from reactor import Reactor
from atom import Atom
from state import State


def uniform_selection(population):
    if len(population) > 1:
        random.sample(population, 2)
    return population


def dummy(x):
    print(json.dumps(x, cls=MyJSON))
    return x

if __name__ == "__main__":

    population = [Atom(i) for i in ['A', 'B', 'C', 'D', 'E']]
    tw = ToyWorld(reactor=Reactor(population=population, reactant_selection=uniform_selection),
                  product_selection=uniform_selection)
    tw.run(generations=20, state=State(persistence=dummy))
