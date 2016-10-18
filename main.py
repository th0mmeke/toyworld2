import random
import os
import logging
import argparse

from toyworld2 import ToyWorld
from uniform_reactant_selection import UniformReactantSelection
from spatial_reactant_selection import SpatialReactantSelection
import uniform_product_selection
import least_energy_product_selection
from molecule import Molecule
from state import State
from semi_realistic_chemistry import SemiRealisticChemistry


def dummy(x):
    #print(json.dumps(x, cls=MyJSON))
    pass


def initialise_logging(args, basedir):
    level = getattr(logging, args.log_level.upper())
    logger = logging.getLogger()
    logger.setLevel(level)

    if not args.log_filename:
        args.log_filename = os.path.basename(basedir) + ".log"
    fh = logging.FileHandler(os.path.join(basedir, args.log_filename))
    fh.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)

    formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s:%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)

def get_args():

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level")
    parent_parser.add_argument('-f', '--log_filename', help="Filename for logging (relative to location of evaluation design file) (optional)")
    parser = argparse.ArgumentParser(parents=[parent_parser])

    return parser.parse_args()

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

    args = get_args()
    initialise_logging(args, os.getcwd())

    chem = SemiRealisticChemistry(bond_energies=bond_energies)

    defn = {"[H][H]": 100, "O=O": 100, "O": 200, "[O-][N+](=O)[N+]([O-])=O": 100, "N(=O)[O]": 100, "O=C=O": 200}
    population = []
    for symbol, quantity in defn.iteritems():
        for i in range(quantity):
            population.append(Molecule(symbol))

    tw = ToyWorld(reactor=SpatialReactantSelection(population=population, ke=10),
                  chemistry=SemiRealisticChemistry(bond_energies=bond_energies),
                  product_selection=least_energy_product_selection.product_selection)
    state = tw.run(generations=200, state=State(persistence=dummy))

