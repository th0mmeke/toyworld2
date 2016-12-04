import os
import json

from evaluator_actual_cycles import EvaluatorActualCycles

datadir = "..\data"

for filename in os.listdir(datadir):
    basename, ext = os.path.splitext(filename)
    evaluator_filename = "{}-actual.json".format(basename)
    if ext == '.json' and not os.path.exists(evaluator_filename):  # in working directory
        with open(os.path.join(datadir, filename)) as f:
            reactions = json.load(f)
        e = EvaluatorActualCycles(reactions=reactions['reactions'])

        print(filename, evaluator_filename, e.g.number_of_nodes(), e.g.number_of_edges())

        population_stoichiometry = []
        count = 0
        for reactant in e.reactants:
            count += 1
            if len(reactant) >= 10:
                print("{}/{}: {}".format(count, len(e.reactants), reactant))
                population_stoichiometry.extend(e.get_reactant_stoichiometry(reactant, minimum_stoichiometry=3, max_depth=10))
        all_cycles = [x['cycle'] for x in population_stoichiometry]

        with open(evaluator_filename, mode='w') as f:
            json.dump(all_cycles, f)
