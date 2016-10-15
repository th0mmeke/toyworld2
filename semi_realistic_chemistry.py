from rdkit.Chem import AllChem as Chem
import copy
import logging
import networkx as nx

from i_chemistry import IChemistry
from reaction import Reaction
from molecule import Molecule


class SemiRealisticChemistry(IChemistry):
    """
    A simple Chemistry based on RDKit.
    """

    def __init__(self, **kwargs):

        try:
            self.bond_energies = kwargs['bond_energies']
        except KeyError:
            raise ValueError('Missing kwarg: bond_energies')

        count = {}
        self.default_bond_energies = {}
        for bond, energy in self.bond_energies.iteritems():
            key = int(bond[1])
            try:
                count[key] += 1
                self.default_bond_energies[key] += energy
            except KeyError:
                count[key] = 1
                self.default_bond_energies[key] = energy

        for i in count.keys():
            self.default_bond_energies[i] = self.default_bond_energies[i] / count[i]

    def enumerate(self, reactants):

        """
        Discover all possible options for reactions between the separate components of this Molecule, based on the
        Chemistry provided. The options are based on possible bond changes:

        1. Breaking of existing bonds
        2. Switch one bond type for another
        3. Addition of bonds

        In each of these cases we update the formal charges to reflect the changes, and test the resulting product for sanity

        If any Hs have been added earlier this will preserve them.

        Note: this is surprisingly difficult in RDKit. RDKit's transformation method - EditableMol - doesn't play well in subclasses as it GetMol() method
        returns an object of class Mol rather than the original type. So using EditableMol results in a different type from the original.
        This means that we either 1) avoid using EditableMol, 2) restore the type, or 3) base everything on Mol rather than subclasses
        3) is impractical - KineticMolecule and Molecule are natural extensions of Mol
        1) is impossible - we need to break bonds, which requires EditableMol (SetBondType() can't remove a bond as Chem.BondType doesn't have a 'nothing' option)
        2) is difficult. We need transformed copies of our original reactants (which are subclasses of Mol), but Mol cannot be copied (no copy method, doesn't
        play with copy.copy or copy.deepcopy as doesn't support Pickle)

        :param reactants: [Molecule]
        :rtype: [Reaction]
        """

        reaction_options = self.get_change_options(reactants)
        addition_options = self.get_addition_options(reactants)
        if len(addition_options) > 0:
            reaction_options.append(addition_options)

        logging.debug("{} reaction options found".format(len(reaction_options)))
        return reaction_options

    def get_change_options(self, reactants):

        """
        :param reactants: [Molecule]
        :return:
        """

        reactant = SemiRealisticChemistry.join(reactants)

        options = []
        for bond in reactant.GetBonds():
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            old_bond_order = int(bond.GetBondType())

            for new_bond_order in range(0, old_bond_order):  # all bond options between none (break bond) and one-less-than-current
                reactant_copy = copy.deepcopy(reactant)  # operate on copy of self, so options are not cumulative
                try:
                    SemiRealisticChemistry._change_bond(reactant_copy, begin_atom_idx, end_atom_idx, new_bond_order)
                except ValueError:
                    pass  # just ignore if this option isn't possible
                else:  # no exception, so managed to execute _change_bond()...
                    bond_energy = self._get_bond_energy(reactant_copy.GetAtomWithIdx(begin_atom_idx).GetSymbol(),
                                                        reactant_copy.GetAtomWithIdx(end_atom_idx).GetSymbol(),
                                                        from_bond_type=old_bond_order,
                                                        to_bond_type=new_bond_order)
                    options.append(Reaction(reactants=reactants,
                                            products=SemiRealisticChemistry.split(reactant_copy),
                                            weight=bond_energy))
        return options

    def get_addition_options(self, reactants):

        """
        :param reactants: Molecule
        :return:
        """

        reactant = SemiRealisticChemistry.join(reactants)
        options = []

        bond_potential = map(lambda x: SemiRealisticChemistry._get_bond_potential(x), reactant.GetAtoms())

        g = nx.Graph()
        for bond in reactant.GetBonds():
            g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            components = list(nx.connected_components(g))
        # Add single atoms as independent components
        for idx in range(reactant.GetNumAtoms()):
            if len(reactant.GetAtomWithIdx(idx).GetBonds()) == 0:
                components.append([idx])

        # Check for possible bonds between all atoms with the potential for one or more additional bonds...
        for begin_atom_idx in range(len(bond_potential)):

            component = [i for i in components if begin_atom_idx in i][0]  # all idx in the same component as begin_atom_idx

            for end_atom_idx in range(begin_atom_idx + 1, len(bond_potential)):

                if end_atom_idx not in component:  # begin_atom and end_atom must join different components

                    # Add options for bonds of various order between begin_atom and end_atom...
                    max_bond_order = min(bond_potential[begin_atom_idx], bond_potential[end_atom_idx])
                    # to meet RDKit restriction from organic self that maximum likely bond is triple-bond
                    max_bond_order = min(max_bond_order, 3)

                    if max_bond_order > 0:

                        for bond_order in range(1, max_bond_order + 1):  # all bond options up to and including max_bond_order
                            reactant_copy = copy.deepcopy(reactant)  # operate on fresh copy so options don't accumulate
                            try:
                                self._change_bond(reactant_copy, begin_atom_idx, end_atom_idx, bond_order)
                            except ValueError:
                                pass  # just ignore invalid options
                            else:
                                bond_energy = self._get_bond_energy(reactant_copy.GetAtomWithIdx(begin_atom_idx).GetSymbol(),
                                                                    reactant_copy.GetAtomWithIdx(end_atom_idx).GetSymbol(),
                                                                    end_bond_type=bond_order)  # bond creation of order bond_order
                                options.append(Reaction(reactants=reactants,
                                                        products=SemiRealisticChemistry.split(reactant_copy),
                                                        weight=bond_energy))
        return options

    @staticmethod
    def split(molecule):
        return [Molecule(smiles) for smiles in Chem.MolToSmiles(molecule).split(".")]

    @staticmethod
    def join(reactants):

        if len(reactants) > 1:
            mols = map(lambda z: Chem.MolFromSmiles(z.get_symbol()), reactants)
            reactant = reduce(lambda x, y: Chem.CombineMols(x, y), mols)
        elif len(reactants) == 1:
            reactant = Chem.MolFromSmiles(reactants[0].get_symbol())

        reactant = Chem.AddHs(reactant)

        return reactant

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
                start_energy = self.bond_energies[atom_1 + str(from_bond_type) + atom_2]
            except KeyError:
                start_energy = self.default_bond_energies[from_bond_type]

        # Energy to create desired bond state
        if to_bond_type <= 0:
            end_energy = 0
        else:
            to_bond_type = min(3, to_bond_type)
            try:
                end_energy = self.bond_energies[atom_1 + str(to_bond_type) + atom_2]
            except KeyError:
                end_energy = self.default_bond_energies[to_bond_type]

        return start_energy - end_energy

    @staticmethod
    def _get_bond_potential(atom):
        """Requires Explicit Hs!

        Simple method based on standard Lewis dot-structures e.g., http://library.thinkquest.org/C006669/data/Chem/bonding/lewis.html

        Bond calculation:
        FC = V - N - B (where single bond = 1, double = 2, etc) to make N = 8
        """

        if atom.GetAtomicNum() == 1:
            if len(atom.GetBonds()) == 0:  # if not already bound...
                return 1
            else:
                return 0
        else:
            bonded_electrons = 0
            for bond in atom.GetBonds():
                bonded_electrons += bond.GetBondType()  # relies on Chem.BondType mapping to int...

            valence_electrons = Chem.GetPeriodicTable().GetNOuterElecs(atom.GetAtomicNum())

            return 8 - (valence_electrons + bonded_electrons + atom.GetFormalCharge())

    @staticmethod
    def _change_bond(mol, begin_atom_idx, end_atom_idx, new_bond_order):

        """
        Change type of bond in place.
        Raise ValueError if cannot be changed.
        """

        if new_bond_order > 3:
            raise ValueError  # to meet RDKit restriction from organic reactions that maximum likely bond is triple-bond

        bond = mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
        e_mol = Chem.EditableMol(mol)
        if bond is not None:  # remove any existing bond
            e_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if new_bond_order != 0:  # add in any new bond
            e_mol.AddBond(begin_atom_idx, end_atom_idx, Chem.BondType.values[new_bond_order])

        # Adjust formal charges if using bond based on charge_order...goal is to make fc of both zero
        begin_atom = mol.GetAtomWithIdx(begin_atom_idx)
        end_atom = mol.GetAtomWithIdx(end_atom_idx)
        begin_atom_fc = begin_atom.GetFormalCharge()
        end_atom_fc = end_atom.GetFormalCharge()

        if cmp(begin_atom_fc, 0) != cmp(end_atom_fc, 0):  # must be opposite sign for this to work
            adjustment = min(abs(begin_atom_fc), abs(end_atom_fc), new_bond_order)
            begin_atom.SetFormalCharge(begin_atom_fc - adjustment * cmp(begin_atom_fc, 0))
            end_atom.SetFormalCharge(end_atom_fc - adjustment * cmp(end_atom_fc, 0))
