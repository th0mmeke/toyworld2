import os
import json

from evaluator_actual_cycles import EvaluatorActualCycles



def evaluate(filename, datadir="C:\Users\Thom\Dropbox/Experiments"):
    basename, ext = os.path.splitext(filename)
    evaluator_filename = os.path.join(datadir, "{}-actual.json".format(basename))

    if ext == '.json' and not os.path.exists(evaluator_filename) and basename[-1:] != 'l' and basename[-1:] != 'y':
        data_filename = os.path.join(datadir, filename)
        print(data_filename)
        with open(data_filename) as f:
            state = json.load(f)
        e = EvaluatorActualCycles(reactions=state['reactions'])

        print(filename, evaluator_filename, e.g.number_of_nodes(), e.g.number_of_edges())

        population_stoichiometry = []
        count = 0
        for reactant in e.reactants:
            count += 1
            if len(reactant) >= 10:
                print("{}/{}: {}".format(count, len(e.reactants), reactant))
                s = e.get_reactant_stoichiometry(reactant, minimum_stoichiometry=2, max_depth=10)
                for item in s:
                    cycle = []
                    for step in item['cycle']:
                        if '+' not in step and '>' not in step and '<' not in step:
                            cycle.append(e.get_smiles(step))
                    print(item['cycle'], cycle)
                    population_stoichiometry.append({'stoichiometry': item['stoichiometry'], 'cycle': cycle})

        with open(evaluator_filename, mode='w') as f:
            json.dump(population_stoichiometry, f)

datadir = "C:\Users\Thom\Dropbox/Experiments"
# evaluate('1480963448-0-21-0.json')
for filename in sorted(os.listdir(datadir), reverse=True):
    evaluate(filename)
