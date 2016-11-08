import json
from evaluator_cycles import EvaluatorCycles

filename = '../data/toyworld2-500000.json'

with open(filename) as f:
    reactions = json.load(f)
e = EvaluatorCycles(reactions=reactions['reactions'])

e.get_population_stoichiometry(max_depth=5, minimum_length=10, minimum_stoichiometry=2)
