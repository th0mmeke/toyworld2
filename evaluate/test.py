import json
from evaluator_cycles import EvaluatorCycles

filename = '../data/toyworld2-500000.json'

with open(filename) as f:
    reactions = json.load(f)
e = EvaluatorCycles(reactions=reactions['reactions'])
