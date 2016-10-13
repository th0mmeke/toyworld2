from rdkit.Chem import AllChem as Chem

from i_chemistry import IChemistry


class SemiRealisticChemistry(IChemistry):

    """A simple Chemistry based on real-world reactions."""

    def __init__(self, bond_energies):

        self.bond_energies = bond_energies

        count = {}
        self.default_bond_energies = {}
        for bond, energy in bond_energies.iteritems():
            key = int(bond[1])
            try:
                count[key] += 1
                self.default_bond_energies[key] += energy
            except KeyError:
                count[key] = 1
                self.default_bond_energies[key] = energy

        for i in count.keys():
            self.default_bond_energies[i] = self.default_bond_energies[i] / count[i]

    @classmethod
    def generate_reaction_options(cls, elements):
        """
        Discover all possible options for reactions between the separate components of this Molecule, based on the
        Chemistry provided. The options are based on possible bond changes; each option is returned in a Reaction object as a list of product Molecules along
        with the associated change in potential energy.

        That is, self.get_potential_energy() + rxn.get_energy_delta() = sum(product molecule potential energy)

        If any Hs have been added earlier this will preserve them.

        Note: this is surprisingly difficult in RDKit. RDKit's transformation method - EditableMol - doesn't play well in subclasses as it GetMol() method
        returns an object of class Mol rather than the original type. So using EditableMol results in a different type from the original.
        This means that we either 1) avoid using EditableMol, 2) restore the type, or 3) base everything on Mol rather than subclasses
        3) is impractical - KineticMolecule and Molecule are natural extensions of Mol
        1) is impossible - we need to break bonds, which requires EditableMol (SetBondType() can't remove a bond as Chem.BondType doesn't have a 'nothing' option)
        2) is difficult. We need transformed copies of our original reactants (which are subclasses of Mol), but Mol cannot be copied (no copy method, doesn't
        play with copy.copy or copy.deepcopy as doesn't support Pickle)

        :param element: Single reactant
        :rtype: list of Reaction
        """

        def _adjust_formal_charge(begin_atom, end_atom, bond_order):
            """Bring the functional charges on the atoms in a bond to as close to zero as possible..."""

            begin_atom_fc = begin_atom.GetFormalCharge()
            end_atom_fc = end_atom.GetFormalCharge()
            # logging.debug("Adjusting test on {} ({}) and {} ({})".format(begin_atom.GetSymbol(), begin_atom_fc, end_atom.GetSymbol(), end_atom_fc))

            if (cmp(begin_atom_fc, 0) != cmp(end_atom_fc, 0)):  # must be opposite sign for this to work
                adjustment = min(abs(begin_atom_fc), abs(end_atom_fc), bond_order)
                begin_atom.SetFormalCharge(begin_atom_fc - adjustment * cmp(begin_atom_fc, 0))
                end_atom.SetFormalCharge(end_atom_fc - adjustment * cmp(end_atom_fc, 0))

        def _change_bond(mol, begin_atom_idx, end_atom_idx, new_bond_order):
            """Change type of bond - throw exception if cannot legitimately be changed.

            :rtype: same type as caller"""

            # logging.debug("Attempting to change bond between {} and {} to {}".format(begin_atom_idx, end_atom_idx, new_bond_order))

            if new_bond_order > 3:
                raise ValueError  # to meet RDKit restriction from organic reactions that maximum likely bond is triple-bond

            bond = mol.GetBondBetweenAtoms(begin_atom_idx, end_atom_idx)
            e_mol = Chem.EditableMol(mol)
            if bond is not None:  # remove any existing bond
                e_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            if new_bond_order != 0:  # add in any new bond
                e_mol.AddBond(begin_atom_idx, end_atom_idx, Chem.BondType.values[new_bond_order])
            new_mol = e_mol.GetMol()

            # Adjust formal charges if using bond based on charge_order...goal is to make fc of both zero
            _adjust_formal_charge(new_mol.GetAtomWithIdx(begin_atom_idx), new_mol.GetAtomWithIdx(end_atom_idx), new_bond_order)

            try:
                new_mol = type(mol)(new_mol, kinetic_energy=mol.get_kinetic_energy())
            except:
                new_mol = type(mol)(new_mol)

            Chem.SanitizeMol(new_mol, catchErrors=False)  # possibly not required
            assert mol.GetNumAtoms(onlyExplicit=False) == new_mol.GetNumAtoms(onlyExplicit=False)  # reaction must preserve matter
            return new_mol  # return same type as caller, self defined in parent method

        element = elements[0].join(elements)
        reaction_options = []

        # Reaction options:
        # 1. breaking of existing bonds
        # 2. addition of bonds
        # 3. switch one bond type for another
        # In each of these cases we update the formal charges to reflect the changes, and test the resulting product for sanity

        # Bond breaking and substitution options...
        for bond in element.GetBonds():
            old_bond_order = int(bond.GetBondType())

            for new_bond_order in range(0,
                                        old_bond_order):  # all bond options between none (break bond) and one-less-than-current
                bond_energy = self._chemistry._get_bond_energy(bond.GetBeginAtomIdx(),
                                                               bond.GetEndAtomIdx(),
                                                               start_bond_type=old_bond_order,
                                                               end_bond_type=new_bond_order)  # change of energy from bond of GetBondType() to bond_order
                if self.weight_option(reaction_energy, bond_energy) > 0:
                    reactants_copy = copy.deepcopy(
                        element)  # operate on copy of self, so options are not cumulative
                    try:
                        new_mol = _change_bond(reactants_copy, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                                               new_bond_order)
                    except ValueError:
                        pass  # just ignore if this option isn't possible
                    else:  # no exception, so managed to execute _change_bond()...
                        reaction_options.append(Reaction(new_mol.split_molecule(), bond_energy))

        # Bond addition options...
        bond_potential = []
        for atom in element.GetAtoms():
            bond_potential.append(self._chemistry._get_bond_potential(atom))

        # Check for possible bonds between all atoms with the potential for one or more additional bonds...
        for begin_atom_idx in range(len(bond_potential)):
            for end_atom_idx in range(begin_atom_idx + 1, len(bond_potential)):
                # logging.debug('{} and {} are in the same component = {}'.format(begin_atom_idx, end_atom_idx, reactant_mol.same_component(begin_atom_idx, end_atom_idx)))
                if not element.same_component(begin_atom_idx, end_atom_idx):

                    # Add options for bonds of various order between begin_atom and end_atom...
                    begin_atom_valence = bond_potential[begin_atom_idx]
                    end_atom_valence = bond_potential[end_atom_idx]
                    max_bond_order = min(begin_atom_valence, end_atom_valence)
                    # to meet RDKit restriction from organic self._chemistry that maximum likely bond is triple-bond
                    max_bond_order = min(max_bond_order, 3)

                    if max_bond_order > 0:

                        for bond_order in range(1, max_bond_order + 1):  # all bond options up to and including max_bond_order
                            reactants_copy = copy.deepcopy(element)  # operate on fresh copy so options don't accumulate
                            try:
                                new_mol = _change_bond(reactants_copy, begin_atom_idx, end_atom_idx, bond_order)
                            except ValueError:
                                pass  # just ignore invalid options
                            else:
                                bond_energy = self._chemistry._get_bond_energy(
                                    new_mol.GetAtomWithIdx(begin_atom_idx), new_mol.GetAtomWithIdx(
                                        end_atom_idx),
                                    end_bond_type=bond_order)  # bond creation of order bond_order
                                logging.debug(
                                    "Reaction option: {} to {} - addition of {} bond from {} to {} (e={})".format(
                                        Chem.MolToSmiles(element), Chem.MolToSmiles(new_mol), bond_order,
                                        begin_atom_idx, end_atom_idx, bond_energy))
                                assert Ulps.almost_equal(element.get_potential_energy(self._chemistry) + bond_energy, Molecule(new_mol).get_potential_energy(self._chemistry))
                                reaction_options.append(Reaction(new_mol.split_molecule(), bond_energy))

        logging.debug("{} reaction options found".format(len(reaction_options)))

        return reaction_options


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

    def _get_bond_energy(self, atom_1, atom_2, to_bond_type=0, from_bond_type=0):
        """
        Returns the energy REQUIRED to make the bond change from start_bond_type (or existing type if not provided) to end_bond_type.
        Creation of a bond requires -e; breaking the bond +e
        Energies taken from http://www.cem.msu.edu/~reusch/OrgPage/bndenrgy.htm - Average Bond Dissociation Enthalpies in kcal per mole

        :param atom_1: One end of the bond
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
                start_energy = self.bond_energies[atom_1.GetSymbol() + str(from_bond_type) + atom_2.GetSymbol()]
            except KeyError:
                start_energy = self.default_bond_energies[from_bond_type]

        # Energy to create desired bond state
        if to_bond_type <= 0:
            end_energy = 0
        else:
            to_bond_type = min(3, to_bond_type)
            try:
                end_energy = self.bond_energies[atom_1.GetSymbol() + str(to_bond_type) + atom_2.GetSymbol()]
            except KeyError:
                end_energy = self.default_bond_energies[to_bond_type]

        return start_energy - end_energy
