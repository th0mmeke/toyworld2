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
from kinetic_reactant_selection import KineticReactantSelection
from uniform_reactant_selection import UniformReactantSelection
from chem_molecule import ChemMolecule
from semi_realistic_chemistry import SemiRealisticChemistry
from state import State
import weighting_functions
import bond_energies


BASE_DIR = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(BASE_DIR):
    BASE_DIR = '/home/cosc/guest/tjy17/Dropbox/Experiments'
BASE_DIR = '.'

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


def get_ar_timeseries2(theta, sd, n):

    """
    Return an AR(1) timeseries.

    :param theta:
    :param sd:
    :param n:
    :return:
    """

    ts = []
    value = 0
    for i in range(0, n+1):
        value = theta * value + random.gauss(0, sd)
        ts.append(value)
    return ts


def run_experiment(filename, population, factors, generations, environment):
    reactor = factors['REACTANT_SELECTION'](population=copy.deepcopy(population))
    initial_population = dict(Counter([str(x) for x in reactor.get_population()]))
    state = State(filename=filename, initial_population=initial_population)

    tw = ToyWorld2(reactor=reactor, product_selection=factors['PRODUCT_SELECTION'], chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies))
    state = tw.run(generations=generations, state=state,  environment_target=factors['ENVIRONMENT_TARGET'], environment_shape=environment)

    final_population = dict(Counter([str(x) for x in reactor.get_population()]))
    state.close(final_population=final_population)


def runner(population, factors, generations, number_of_repeats, number_of_environments):

    filebase = int(time.time())
    experiment_number = 0

    experiments = (dict(itertools.izip(factors, x)) for x in itertools.product(*factors.itervalues()))
    total_experiments = number_of_repeats * number_of_environments * len(list(experiments))

    with open(os.path.join(BASE_DIR, '{}-metadata.csv'.format(filebase)), 'w', buffering=0) as f:

        for experiment in (dict(itertools.izip(factors, x)) for x in itertools.product(*factors.itervalues())):
            for environment_number in range(0, number_of_environments):

                if experiment['ENVIRONMENT_SHAPE'] == 'BISTATE' or environment_number == 0:
                    environment_specification = (0, 0)
                else:
                    environment_specification = random.uniform(0.0, 1.0), random.uniform(0, 20)

                if experiment['ENVIRONMENT_SHAPE'] == 'BISTATE' and environment_number > 0:
                    environment = list(itertools.chain.from_iterable([[-20]*10000 + [20]*10000] * (int(generations/20000)+1)))
                else:
                    environment = get_ar_timeseries2(*environment_specification, n=generations)
                assert len(environment) >= generations

                metadata = [str(experiment_number)]
                metadata.extend([experiment['REACTANT_SELECTION'].__name__, experiment['PRODUCT_SELECTION'].__name__])
                metadata.extend([experiment['ENVIRONMENT_TARGET'], experiment['ENVIRONMENT_SHAPE']])
                metadata.extend([str(x) for x in environment_specification])
                metadata.append(str(nolds.dfa(environment)))
                metadata.append(str(nolds.sampen(environment)))

                f.write(','.join(metadata) + "\n")

                for repeat_number in range(0, number_of_repeats):
                    print("{0}/{1}".format((experiment_number*number_of_repeats) + repeat_number + 1, total_experiments))
                    filename = "{}-{}-{}-{}.json".format(filebase, experiment_number, environment_number, repeat_number)
                    run_experiment(os.path.join(BASE_DIR, filename), population, experiment, generations, environment=environment)

            experiment_number += 1

if __name__ == "__main__":

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-l', '--log_level', default='INFO', help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL")
    parent_parser.add_argument('-f', '--log_filename', help="Filename for logging (relative to location of evaluation design file) (optional)")
    parent_parser.add_argument('-g', '--generations', type=int, help="Number of generations", default=25000)
    parser = argparse.ArgumentParser(parents=[parent_parser])

    args = parser.parse_args()
    initialise_logging(args, os.getcwd())

    defn = {"[H][H]": 10, "FO": 10, "O": 20, "[O-][N+](=O)[N+]([O-])=O": 10, "N(=O)[O]": 10, "O=C=O": 20}
    population = []
    for symbol, quantity in defn.iteritems():
        for i in range(quantity):
            population.append(ChemMolecule(symbol))

    logging.info("Generations: {}".format(args.generations))

    factors = {
        'REACTANT_SELECTION': [KineticReactantSelection],
        'PRODUCT_SELECTION': [weighting_functions.least_energy_weighting, weighting_functions.uniform_weighting],
        'ENVIRONMENT_SHAPE': ['AR', 'BISTATE'],
        'ENVIRONMENT_TARGET': ['POPULATION', 'KE']
    }

    runner(population, factors, generations=args.generations, number_of_repeats=1, number_of_environments=3)

