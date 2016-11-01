from rdkit.Chem import AllChem as Chem
import copy
import logging
import networkx as nx
import re
import os
import sys
from contextlib import contextmanager

from i_chemistry import IChemistry
from reaction import Reaction
from chem_molecule import ChemMolecule
from ulps import Ulps


@contextmanager
def stderr_redirected():

    """
    Suppress error messages from Chem.SanitizeMol(), produced even if rdkit logging is set to critical.
    Based on stdout_redirected() from http://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python/22434262#22434262
    """

    stderr = sys.stderr
    stderr_fd = sys.stderr.fileno()
    # copy stderr_fd before it is overwritten
    # NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stderr_fd), 'wb') as copied:
        stderr.flush()  # flush library buffers that dup2 knows nothing about
        with open(os.devnull, 'wb') as to_file:
            os.dup2(to_file.fileno(), stderr_fd)  # $ exec > to
        try:
            yield stderr  # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            # NOTE: dup2 makes stderr_fd inheritable unconditionally
            stderr.flush()
            os.dup2(copied.fileno(), stderr_fd)  # $ exec >&copied


class SemiRealisticChemistry(IChemistry):
    """
    A simple Chemistry based on RDKit. All molecules must be ChemMolecules to take advantage of bond manipulation
    methods in RDKit.
    """

    def __init__(self, **kwargs):

        try:
            self._bond_energies = kwargs['bond_energies']
        except KeyError:
            raise ValueError('Missing kwarg: bond_energies')

        count = {}
        self._default_bond_energies = {}
        for bond, energy in self._bond_energies.iteritems():
            key = int(re.match(r"(\w+)(\d)(\w+)", bond).groups()[1])
            try:
                count[key] += 1
                self._default_bond_energies[key] += energy
            except KeyError:
                count[key] = 1
                self._default_bond_energies[key] = energy

        for i in count.keys():
            self._default_bond_energies[i] = self._default_bond_energies[i] / count[i]

    def enumerate(self, partial_reaction):

        """
        Discover all possible options for reactions between the separate components of this ChemMolecule.
        The options are based on possible bond changes:

        1. Breaking of existing bonds
        2. Switch one bond type for another
        3. Addition of bonds

        :param partial_reaction: Reaction containing ChemMolecule reactants and reactant_value
        :rtype: [Reaction]
        """

        if partial_reaction is None:
            logging.debug("No reactants to enumerate")
            return None

        reaction_options = self._get_change_options(partial_reaction)
        addition_options = self._get_addition_options(partial_reaction)
        if len(addition_options) > 0:
            reaction_options.extend(addition_options)

        logging.debug("{} reaction options found".format(len(reaction_options)))

        return reaction_options

    def _get_change_options(self, partial_reaction):

        """
        :param partial_reaction: Reaction
        :return: [Reaction]
        """

        reactant = SemiRealisticChemistry._join(partial_reaction.get_reactants())

        options = []
        for bond in reactant.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            old_bond_order = int(bond.GetBondType())

            for new_bond_order in range(0, 3):  # all bond options between none (break bond) and one-less-than-current
                try:
                    product = SemiRealisticChemistry._change_bond(copy.deepcopy(reactant), begin_atom_idx, end_atom_idx, new_bond_order)
                except ValueError:
                    pass  # just ignore if this option isn't possible
                else:  # no exception, so managed to execute _change_bond()...
                    with stderr_redirected():
                        Chem.SanitizeMol(product)
                    bond_energy = self._get_bond_energy(product.GetAtomWithIdx(begin_atom_idx).GetSymbol(),
                                                        product.GetAtomWithIdx(end_atom_idx).GetSymbol(),
                                                        from_bond_type=old_bond_order,
                                                        to_bond_type=new_bond_order)
                    options.append(Reaction(reactants=partial_reaction.get_reactants(),
                                            reactant_value=partial_reaction.reactant_value,
                                            products=SemiRealisticChemistry._split(product),
                                            product_value=bond_energy))
        return options

    def _get_addition_options(self, partial_reaction):

        """
        :param partial_reaction: Reaction
        :return: [Reaction]
        """

        reactant = SemiRealisticChemistry._join(partial_reaction.get_reactants())
        options = []

        g = nx.Graph()
        for bond in reactant.GetBonds():
            g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        components = list(nx.connected_components(g))
        # Add single atoms as independent components
        for idx in range(reactant.GetNumAtoms(onlyExplicit=False)):
            if len(reactant.GetAtomWithIdx(idx).GetBonds()) == 0:
                components.append([idx])

        # Check for possible bonds between all atoms with the potential for one or more additional bonds...
        for begin_atom in reactant.GetAtoms():

            component = [i for i in components if begin_atom.GetIdx() in i][0]  # all idx in the same component as begin_atom_idx

            for end_atom_idx in range(begin_atom.GetIdx() + 1, reactant.GetNumAtoms(onlyExplicit=False)):
                if end_atom_idx not in component:  # begin_atom and end_atom must join different components
                    for bond_order in range(1, 4):  # all bond options up to and including max_bond_order
                        try:
                            product = self._change_bond(copy.deepcopy(reactant), begin_atom.GetIdx(), end_atom_idx, bond_order)

                        except ValueError:
                            pass  # just ignore invalid options
                        else:
                            bond_energy = self._get_bond_energy(begin_atom.GetSymbol(),
                                                                product.GetAtomWithIdx(end_atom_idx).GetSymbol(),
                                                                to_bond_type=bond_order)  # bond creation of order bond_order
                            options.append(Reaction(reactants=partial_reaction.get_reactants(),
                                                    reactant_value=partial_reaction.reactant_value,
                                                    products=SemiRealisticChemistry._split(product),
                                                    product_value=bond_energy))
        return options

    @staticmethod
    def _split(molecule):
        """
        Split a molecule at the '.' symbols in the SMILES representation.

        Raise ValueError if no valid split of molecule.

        :param molecule: RDKit.Mol with explicit Hs
        :return: ChemMolecule
        """

        initial_mass = sum([atom.GetMass() for atom in molecule.GetAtoms()])
        mols = Chem.GetMolFrags(molecule, asMols=True)

        final_mass = sum([sum([atom.GetMass() for atom in molecule.GetAtoms()]) for molecule in mols])
        assert Ulps.almost_equal(initial_mass, final_mass)

        mols = [ChemMolecule(Chem.MolToSmiles(x)) for x in mols]  # can raise ValueError here even if molecule has passed validity check in _change_bond()
        assert Ulps.almost_equal(initial_mass, sum([m.mass for m in mols]))

        return mols

    @staticmethod
    def _join(molecules):
        """
        Combine reactants to form a single RDKit Mol. Assumes that get_symbol() for each
        reactant returns a string containing the SMILES representation for that Molecule.

        :param molecules: Molecule
        :return: RDKit.Mol
        """

        if len(molecules) > 1:
            mols = map(lambda z: Chem.AddHs(Chem.MolFromSmiles(z.get_symbol())), molecules)
            mol = reduce(lambda x, y: Chem.CombineMols(x, y), mols)
        elif len(molecules) == 1:
            mol = Chem.AddHs(Chem.MolFromSmiles(molecules[0].get_symbol()))
            mols = [mol]

        assert sum([m.GetNumAtoms(onlyExplicit=False) for m in mols]) == mol.GetNumAtoms(onlyExplicit=False)
        assert Ulps.almost_equal(sum([atom.GetMass() for atom in mol.GetAtoms()]), sum([r.mass for r in molecules]))

        return mol

    def _get_bond_energy(self, atom_1, atom_2, to_bond_type=0, from_bond_type=0):
        """
        Returns the energy REQUIRED to make the bond change from start_bond_type (or existing type if not provided) to end_bond_type.
        Creation of a bond requires -e; breaking the bond +e
        Energies taken from http://www.cem.msu.edu/~reusch/OrgPage/bndenrgy.htm - Average Bond Dissociation Enthalpies in kcal per mole

        :param atom_1: Symbol for one end of the bond
        :param atom_2: Other end of the bond
        :param from_bond_type: Type of any initial bond, corresponding to index into Chem.BondType.values
        :param to_bond_type: Type of any resulting bond, corresponding to index into Chem.BondType.values
        """

        # Energy to release current bond state
        if from_bond_type <= 0:
            start_energy = 0
        else:
            from_bond_type = min(3, from_bond_type)
            try:
                start_energy = self._bond_energies[atom_1 + str(from_bond_type) + atom_2]
            except KeyError:
                start_energy = self._default_bond_energies[from_bond_type]

        # Energy to create desired bond state
        if to_bond_type <= 0:
            end_energy = 0
        else:
            to_bond_type = min(3, to_bond_type)
            try:
                end_energy = self._bond_energies[atom_1 + str(to_bond_type) + atom_2]
            except KeyError:
                end_energy = self._default_bond_energies[to_bond_type]

        return start_energy - end_energy

    @staticmethod
    def _change_bond(mol, begin_atom_idx, end_atom_idx, new_bond_order):

        """
        Change type of bond in place.
        Raise ValueError if cannot be changed.

        :return: RDKit.mol
        """

        if new_bond_order > 3:
            raise ValueError  # to meet RDKit restriction from organic reactions that maximum likely bond is triple-bond

        bond = mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
        new_mol = Chem.RWMol(mol)

        begin_atom = new_mol.GetAtomWithIdx(begin_atom_idx)
        end_atom = new_mol.GetAtomWithIdx(end_atom_idx)
        begin_atom_fc = begin_atom.GetFormalCharge()
        end_atom_fc = end_atom.GetFormalCharge()

        if bond is not None:  # remove any existing bond
            if begin_atom.GetAtomicNum() == 1 and end_atom.GetAtomicNum() != 1:
                begin_atom.SetFormalCharge(1)
                end_atom.SetFormalCharge(-1)
            elif end_atom.GetAtomicNum() == 1 and begin_atom.GetAtomicNum() != 1:
                begin_atom.SetFormalCharge(-1)
                end_atom.SetFormalCharge(1)
            new_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        if new_bond_order != 0:  # add in any new bond
            new_mol.AddBond(begin_atom_idx, end_atom_idx, Chem.BondType.values[new_bond_order])
            # Adjust formal charges when bonding ions...goal is to make fc of both zero

            if cmp(begin_atom_fc, 0) != cmp(end_atom_fc, 0):  # must be opposite sign for this to work
                adjustment = min(abs(begin_atom_fc), abs(end_atom_fc), new_bond_order)
                begin_atom.SetFormalCharge(begin_atom_fc - adjustment * cmp(begin_atom_fc, 0))
                end_atom.SetFormalCharge(end_atom_fc - adjustment * cmp(end_atom_fc, 0))

        new_mol = new_mol.GetMol()

        with stderr_redirected():
            Chem.SanitizeMol(new_mol)  # Raises ValueError if an invalid molecule can't be corrected

        # Removing a bond between H and another atom can add in extra Hs...for conservation of mass reject these:
        if new_mol.GetNumAtoms(onlyExplicit=False) != mol.GetNumAtoms(onlyExplicit=False):
            raise ValueError(Chem.MolToSmiles(new_mol))
        if Chem.MolFromSmiles(Chem.MolToSmiles(new_mol)) is None:
            raise ValueError

        return new_mol
