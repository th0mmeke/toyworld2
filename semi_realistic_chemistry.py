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

    @staticmethod
    def get_bond_potential(atom):

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

    def get_bond_energy(self, atom_1, atom_2, to_bond_type=0, from_bond_type=0):

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
