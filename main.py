import random
import json
from my_json import MyJSON

from toyworld import ToyWorld
from reactant_selection import Reactor
from molecule import Molecule
from state import State
from semi_realistic_chemistry import SemiRealisticChemistry


def uniform_selection(population):
    if len(population) > 1:
        random.sample(population, 2)
    return population


def dummy(x):
    print(json.dumps(x, cls=MyJSON))
    return x

if __name__ == "__main__":
    bond_energies = {
        'H1H': 104.2,
        'C1C': 83,
        'N1N': 38.4,
        'O1O': 35,
        'H1C': 99, 'C1H': 99,
        'H1N': 93, 'N1H': 93,
        'H1O': 111, 'O1H': 111,
        'C1N': 73, 'N1C': 73,
        'C1O': 85.5, 'O1C': 85.5,
        'N1O': 55, 'O1N': 55,
        'C2O': 185, 'O2C': 185,  # rough average of range
        'C2C': 146,
        'N2N': 149,
        'O2O': 119,
        'C2N': 147, 'N2C': 147,
        'N2O': 143, 'O2N': 143,
        'C3O': 258, 'O3C': 258,
        'C3C': 200,
        'N3N': 226,
        'C3N': 213, 'N3C': 213,
        'C4C': 200  # theoretically possible from valences, but in nature forms a C2C bond instead
    }

    chem = SemiRealisticChemistry(bond_energies=bond_energies)

    population = [Molecule(i) for i in ['C', 'O', 'C', 'H', 'H']]
    tw = ToyWorld(reactor=Reactor(population=population, reactant_selection=uniform_selection),
                  generate_reaction_set=chem.generate_reaction_options,
                  product_selection=uniform_selection)
    tw.run(generations=20, state=State(persistence=dummy))
