from collections import Counter
import os
import logging
import argparse
import random
import itertools
import nolds

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

MAX_SD = 0.4
N_REPEATS = 5
N_ENVIRONMENTS = 10


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


def get_ar_timeseries(theta, sd, bias, generations):
    ts = []
    value = 0
    t = random.randint(-100, -50)  # initial burn-in period
    for i in range(t, generations):
        value = theta * value + random.gauss(0, sd) + bias  # error term with mean = 0 and sd = sd
        if i >= 0:
            ts.append(value)
    return ts

def get_environment_specification():
    # Environment change is modelled as a change in the relationship between entity and environment.

    # Entity is represented by a fitness and a fidelity (to parent)
    # The difficulty is in how to connect fidelity and fitness (both composites) to environmental change.
    # Without knowing the specifics of the change, and modelling the reaction of each entity to that change,
    # we can only model abstracts.

    # Parental relationships provide a lineage, where related entities should have a similar response to change
    # In other words, environmental change is a change in fitness, but one where the change is conditioned
    # by lineage.

    for i in range(N_ENVIRONMENTS):
        yield random.uniform(-MAX_SD, MAX_SD), random.uniform(0, MAX_SD), random.uniform(-MAX_SD/10, MAX_SD/10)


def runner(population, generations, number_of_repeats, number_of_environments):

    factor_defns = {
        'REACTANT_SELECTION': [LocalReactantSelection, UniformReactantSelection, SpatialReactantSelection, LengthBiasedReactantSelection],
        'PRODUCT_SELECTION': [weighting_functions.least_energy_weighting, weighting_functions.uniform_weighting],
    }

    experiment_number = run_number = 0
    total_experiments = number_of_repeats * number_of_environments * len(factor_defns['REACTANT_SELECTION'])*len(factor_defns['PRODUCT_SELECTION'])

    with open('data/metadata.csv', 'w') as f:
        for experiment in itertools.product(factor_defns['REACTANT_SELECTION'], factor_defns['PRODUCT_SELECTION']):

            for environment_specification in get_environment_specification():

                environment = get_ar_timeseries(*environment_specification, generations=generations)
                metadata = [str(experiment_number), experiment[0].__name__, experiment[1].__name__]
                metadata.extend([str(x) for x in environment_specification])
                metadata.append(str(nolds.hurst_rs(environment)))
                metadata.append(str(nolds.dfa(environment)))
                metadata.append(str(nolds.sampen(environment)))
                metadata.append(str(nolds.corr_dim(environment, 2)))
                print(','.join(metadata))
                #f.write(','.join(metadata))

                for repeat in range(0, number_of_repeats):

                    print("{0}/{1}".format((experiment_number*number_of_repeats) + repeat + 1, total_experiments))

                    # reactor = experiment[0](population=copy.deepcopy(population))
                    # initial_population = dict(Counter([str(x) for x in reactor.get_population()]))

                    # tw = ToyWorld2(reactor=reactor, product_selection=experiment[1], chemistry=SemiRealisticChemistry(bond_energies=bond_energies.bond_energies))

                    # state = State(filename="data/env-{}-{}.json".format(experiment_number, repeat), initial_population=initial_population)
                    # state = tw.run(generations=args.generations, state=state, environments=environment)

                    # final_population = dict(Counter([str(x) for x in reactor.get_population()]))
                    # state.close(final_population=final_population)

                experiment_number += 1


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

    # adenine, guanine, cytosine, uracil
    adenine = "Nc1ncnc2[nH]cnc12"  # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16708
    cytosine = "Nc1cc[nH]c(=O)n1" # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=16040
    guanine = "Nc1nc2[nH]cnc2c(=O)[nH]1"  # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:16235
    uracil = "O=c1cc[nH]c(=O)[nH]1"  # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17568

    rna_bases = [adenine, cytosine, guanine, uracil]
    #defn = {x: 100 for x in rna_bases}

    population = []
    for symbol, quantity in defn.iteritems():
        for i in range(quantity):
            population.append(ChemMolecule(symbol))

    logging.info("Generations: {}".format(args.generations))
    runner(population, generations=args.generations, number_of_repeats=N_REPEATS, number_of_environments=N_ENVIRONMENTS)


