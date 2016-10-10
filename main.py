import random
from toyworld import ToyWorld
from reactor import Reactor
from atom import Atom


def uniform_selection(x):
    if len(x) > 1:
        random.sample(x, 2)
    return x

if __name__ == "__main__":

    population = [Atom(i) for i in ['A', 'B', 'C', 'D', 'E']]
    tw = ToyWorld(reactor=Reactor(population=population, reactant_selection=uniform_selection),
                  product_selection=uniform_selection)
    tw.run(generations=20)
