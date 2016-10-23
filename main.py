from collections import Counter
import os
import logging
import argparse

from toyworld2 import ToyWorld2
from spatial_reactant_selection import SpatialReactantSelection
import uniform_product_selection
from chem_molecule import ChemMolecule
from state import State
from semi_realistic_chemistry import SemiRealisticChemistry
import bond_energies


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
    parent_parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL")
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

    reactor = SpatialReactantSelection(population=population, ke=100)
    tw = ToyWorld2(reactor=reactor,
                   chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies),
                   product_selection=uniform_product_selection.product_selection)

    logging.info("Initial population: {}".format(Counter([str(x) for x in reactor.get_population()])))
    state = tw.run(generations=5000, state=State(filename="toyworld2.json"))
    logging.info("Final population: {}".format(Counter([str(x) for x in reactor.get_population()])))
