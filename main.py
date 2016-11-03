from collections import Counter
import os
import logging
import argparse

from toyworld2 import ToyWorld2
from length_biased_reactant_selection import LengthBiasedReactantSelection
from chem_molecule import ChemMolecule
from semi_realistic_chemistry import SemiRealisticChemistry
from state import State
import weighting_functions
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
    parent_parser.add_argument('-g', '--generations', type=int, help="Number of generations", default=500000)
    parser = argparse.ArgumentParser(parents=[parent_parser])

    args = parser.parse_args()
    initialise_logging(args, os.getcwd())

    #defn = {"[H][H]": 10, "FO": 10, "O": 20, "[O-][N+](=O)[N+]([O-])=O": 10, "N(=O)[O]": 10, "O=C=O": 20}
    proteinogenic_amino_acids = ["O=C(O)[C@@H](N)C",  # ALA
                   "NC(CCCNC(N)=N)C(O)=O",  # ARG
                   "O=C(N)C[C@H](N)C(=O)O",  # ASN
                   "O=C(O)CC(N)C(=O)O",  # ASP
                   "C([C@@H](C(=O)O)N)S",  # CYS
                   "C(CC(=O)O)C(C(=O)O)N",  # GLU
                   "O=C(N)CCC(N)C(=O)O",  # GLN
                   "C(C(=O)O)N",  # GLY
                   "O=C([C@H](CC1=CNC=N1)N)O",  # HIS
                   "CC[C@H](C)[C@@H](C(=O)O)N",  # ILE
                   "CC(C)C[C@@H](C(=O)O)N",  # LEU
                   "C(CCN)CC(C(=O)O)N",  # LYS
                   "CSCCC(C(=O)O)N",  # MET
                   "C1=CC=C(C=C1)CC(C(=O)O)N",  # PHE
                   "OC(=O)C1CCCN1",  # PRO
                   "C([C@@H](C(=O)O)N)O",  # SER
                   "C[C@H]([C@@H](C(=O)O)N)O",  # THR
                   "c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N",  # TRP
                   "N[C@@H](Cc1ccc(O)cc1)C(O)=O",  # TYR
                   "CC(C)[C@@H](C(=O)O)N",  # VAL
                   ]
    defn = {x: 40 for x in proteinogenic_amino_acids}

    population = []
    for symbol, quantity in defn.iteritems():
        for i in range(quantity):
            population.append(ChemMolecule(symbol))

    reactor = LengthBiasedReactantSelection(population=population, ke=100)
    tw = ToyWorld2(reactor=reactor,
                   chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies),
                   product_selection=weighting_functions.replicant_weighting)

    logging.info("Generations: {}".format(args.generations))
    logging.info("Initial population: {}".format(Counter([str(x) for x in reactor.get_population()])))
    state = State(filename="data/toyworld2.json", initial_population=reactor.get_population())
    state = tw.run(generations=args.generations, state=state)
    state.close(final_population=reactor.get_population())
    logging.info("Final population: {}".format(Counter([str(x) for x in reactor.get_population()])))

