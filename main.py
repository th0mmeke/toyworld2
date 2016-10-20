import random
import os
import logging
import argparse

from toyworld2 import ToyWorld2
from uniform_reactant_selection import UniformReactantSelection
from spatial_reactant_selection import SpatialReactantSelection
import uniform_product_selection
import least_energy_product_selection
from molecule import Molecule
from chem_molecule import ChemMolecule
from state import State
from semi_realistic_chemistry import SemiRealisticChemistry
import bond_energies


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


if __name__ == "__main__":

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level")
    parent_parser.add_argument('-f', '--log_filename', help="Filename for logging (relative to location of evaluation design file) (optional)")
    parser = argparse.ArgumentParser(parents=[parent_parser])

    args = parser.parse_args()
    initialise_logging(args, os.getcwd())

    chem = SemiRealisticChemistry(bond_energies=bond_energies.bond_energies)

    defn = {"[H][H]": 100, "FO": 100, "O": 200, "[O-][N+](=O)[N+]([O-])=O": 100, "N(=O)[O]": 100, "O=C=O": 200}
    population = []
    for symbol, quantity in defn.iteritems():
        for i in range(quantity):
            population.append(ChemMolecule(symbol))

    tw = ToyWorld2(reactor=SpatialReactantSelection(population=population, ke=100),
                   chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies),
                   product_selection=uniform_product_selection.product_selection)
    state = tw.run(generations=5, state=State(persistence=dummy))

    #tw = ToyWorld(reactor=SpatialReactantSelection(population=population, ke=100),
    #                chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies),
    #                product_selection=least_energy_product_selection.product_selection)
