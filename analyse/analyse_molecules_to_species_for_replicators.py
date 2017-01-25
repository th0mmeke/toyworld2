import json
import os
import glob

def map_id_to_smiles(molecule_cycle, smiles):
    smiles_cycle = []
    for step in molecule_cycle:
        new_step = step
        for id, s in smiles.iteritems():
            new_step = new_step.replace(id, s)
        smiles_cycle.append(new_step)
    assert len(molecule_cycle) == len(smiles_cycle)
    return smiles_cycle


def load_smiles(reactions):
    mapping = {}
    for reaction in reactions:
        for id, smiles in reaction['reactants'].iteritems():
            mapping[id] = smiles
        for id, smiles in reaction['products'].iteritems():
            mapping[id] = smiles
    return mapping

datadir = 'C:\Users\Thom\Dropbox/Experiments'
if not os.path.isdir(datadir):
    datadir = '/home/cosc/guest/tjy17/Dropbox/Experiments'
filebase = '1484540618'

# Construct list of stable cycles per environment

for filename in glob.glob(os.path.join(datadir, filebase+'*replicators.json')):

    basename, ext = os.path.splitext(filename)
    print(filename)
    datetime, experiment, repeat, bistate, dummy2, dummy3 = basename.split('-')
    with open(os.path.join(datadir, filename)) as f:
        replicators = json.load(f)
    with open(os.path.join(datadir, '{}-{}-{}-{}.json'.format(datetime, experiment, repeat, bistate))) as f:
        state = json.load(f)
        smiles = load_smiles(state['reactions'])

    species_replicators = []
    for replicator in replicators:
        species_cycles = []
        for molecular_cycle in replicator:
            species_cycles.append(map_id_to_smiles(molecular_cycle, smiles))
        species_replicators.append(species_cycles)
    # print(species_replicator)
    # print([len(x) for x in replicators], [len(x) for x in species_replicators])

    evaluator_filename = os.path.join(datadir, '{}-species.json'.format(basename))
    print(evaluator_filename)
    with open(evaluator_filename, mode='w') as f:
        json.dump(species_replicators, f)