from collections import Counter
import os
import logging
import argparse
import random
import itertools
import nolds
import time
import copy

from toyworld2 import ToyWorld2
from length_biased_reactant_selection import LengthBiasedReactantSelection
from uniform_reactant_selection import UniformReactantSelection
from spatial_reactant_selection import SpatialReactantSelection
from local_reactant_selection import LocalReactantSelection
from chem_molecule import ChemMolecule
from semi_realistic_chemistry import SemiRealisticChemistry
from state import State
import weighting_functions
import bond_energies


BASE_DIR = "C:\Users\Thom\Dropbox\Experiments"

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


def get_ar_timeseries(theta, sd, initial, generations):

    """
    Return an integer AR(1) timeseries.

    :param theta:
    :param sd:
    :param initial:
    :param generations:
    :return:
    """

    ts = []
    value = 0
    for i in range(0, generations+1):
        value = theta * value + random.gauss(0, sd)
        ts.append(max(0, int(value + initial)))  # lower bound of zero
    return ts

def get_ar_timeseries2(theta, sd, generations):

    """
    Return an integer AR(1) timeseries.

    :param theta:
    :param sd:
    :param initial:
    :param generations:
    :return:
    """

    ts = []
    value = 0
    for i in range(0, generations+1):
        value = theta * value + random.gauss(0, sd)
        ts.append(value)
    return ts

def run_experiment(filename, population, experiment, generations, environment=None):
    reactor = experiment[0](population=copy.deepcopy(population))
    initial_population = dict(Counter([str(x) for x in reactor.get_population()]))

    tw = ToyWorld2(reactor=reactor, product_selection=experiment[1], chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies))

    state = State(filename=filename, initial_population=initial_population)
    state = tw.run(generations=generations, state=state, environment=environment)

    final_population = dict(Counter([str(x) for x in reactor.get_population()]))
    state.close(final_population=final_population)


def runner(population, factors, generations, number_of_repeats, number_of_environments):

    filebase = int(time.time())
    experiment_number = 0
    total_experiments = number_of_repeats * number_of_environments * len(factors['REACTANT_SELECTION'])*len(factors['PRODUCT_SELECTION'])

    with open(os.path.join(BASE_DIR, '{}-metadata.csv'.format(filebase)), 'w', 0) as f:

        for experiment in itertools.product(factors['REACTANT_SELECTION'], factors['PRODUCT_SELECTION']):

            for environment_number in range(0, number_of_environments):

                # environment_specification = (theta, sd, bias): f(t) = theta * f(t-1) + random.gauss(0, sd)
                if environment_number == 0:
                    environment_specification = (0, 0)
                else:
                    environment_specification = random.uniform(0.8, 1.0), random.uniform(0, 10)
                # environment = get_ar_timeseries(*environment_specification, initial=len(population), generations=generations)
                environment = get_ar_timeseries2(*environment_specification, generations=generations)
                print(environment)

                metadata = [str(experiment_number), experiment[0].__name__, experiment[1].__name__]
                metadata.extend([str(x) for x in environment_specification])
                metadata.append(str(nolds.hurst_rs(environment)))
                metadata.append(str(nolds.dfa(environment)))
                metadata.append(str(nolds.sampen(environment)))

                f.write(','.join(metadata) + "\n")

                for repeat_number in range(0, number_of_repeats):
                    print("{0}/{1} with {2}".format((experiment_number*number_of_repeats) + repeat_number + 1, total_experiments, environment_specification))
                    filename = "{}-energy-{}-{}-{}.json".format(filebase, experiment_number, environment_number, repeat_number)
                    run_experiment(os.path.join(BASE_DIR, filename), population, experiment, generations, environment=environment[:])

            experiment_number += 1


if __name__ == "__main__":

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL")
    parent_parser.add_argument('-f', '--log_filename', help="Filename for logging (relative to location of evaluation design file) (optional)")
    parent_parser.add_argument('-g', '--generations', type=int, help="Number of generations", default=50000)
    parser = argparse.ArgumentParser(parents=[parent_parser])

    args = parser.parse_args()
    initialise_logging(args, os.getcwd())

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

    # adenine, guanine, cytosine, uracil
    adenine = "Nc1ncnc2[nH]cnc12"  # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16708
    cytosine = "Nc1cc[nH]c(=O)n1" # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=16040
    guanine = "Nc1nc2[nH]cnc2c(=O)[nH]1"  # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16235
    uracil = "O=c1cc[nH]c(=O)[nH]1"  # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17568

    rna_bases = [adenine, cytosine, guanine, uracil]
    #defn = {x: 100 for x in rna_bases}

    defn = {"[H][H]": 10, "FO": 10, "O": 20, "[O-][N+](=O)[N+]([O-])=O": 10, "N(=O)[O]": 10, "O=C=O": 20}
    population = []
    for symbol, quantity in defn.iteritems():
        for i in range(quantity):
            population.append(ChemMolecule(symbol))

    logging.info("Generations: {}".format(args.generations))

    factors = {
        'REACTANT_SELECTION': [LocalReactantSelection],
        'PRODUCT_SELECTION': [weighting_functions.least_energy_weighting],
    }

    #experiment = [LocalReactantSelection, weighting_functions.least_energy_weighting]
    #run_experiment('data/local_unchanging.json', experiment=experiment, population=population, generations=args.generations)

    runner(population, factors, generations=args.generations, number_of_repeats=1, number_of_environments=10)

