import logging
import math

from rdkit.Chem import AllChem as Chem
import networkx as nx
import pymunk as pm


class KineticMolecule(Chem.Mol):

    """A base representation of a Molecule with potential and kinetic energy.

    * Potential energy is the energy required to form all bonds in the Molecule (therefore always negative as bond formation releases energy).
    * Kinetic energy is, as usual, equal to 1/2 * mass * velocity ^ 2."""

    def __init__(self, source, internal_energy=0, kinetic_energy=0, canonize=True, components=None):
        """
        :param internal_energy: initial internal energy for this Molecule
        :param kinetic_energy: initial kinetic energy for this Molecule
        :param canonize: make implicit Hs explicit? Default True, but when copying Molecules we don't want to chance that these changes might be introduced
        :param components: for molecules that consist of multiple disjoint components, a mapping of atoms to component
        """

        # Make all H explicit for later processing

        if not isinstance(source, Chem.Mol):
            source = Chem.MolFromSmiles(source)

        if canonize:
            source = Chem.AddHs(source)

        Chem.Mol.__init__(self, source.ToBinary())

        self._smiles = Chem.MolToSmiles(self)

        if components is None and self._smiles.find(".") == -1:  # if simple molecule with only one component, set components manually
            components = [set(range(source.GetNumAtoms()))]

        self.radius = 1
        self._components = components
        self._mass = sum([atom.GetMass() for atom in self.GetAtoms()])

        inertia = pm.moment_for_circle(self._mass, 0, self.radius, (0, 0))
        self.body = pm.Body(self._mass, inertia)
        shape = pm.Circle(self.body, self.radius, (0, 0))
        shape.elasticity = 0.999  # required for standard PyMunk collision handler to do a 'perfect' bounce
        shape.friction = 0.0
        self._shapes = [shape]

        self.set_kinetic_energy(kinetic_energy)
        self.set_internal_energy(internal_energy)

    def set_position(self, x, y):
        self.body.position = x, y

    def get_position(self):
        return self.body.position

    def set_kinetic_energy(self, ke):
        new_speed = math.sqrt(2.0 * ke / self.get_mass())  # v = sqrt(2ke/m)
        try:
            self.body.velocity = self.body.velocity * new_speed / self.get_speed()
        except:
            self.body.velocity = Kinetics2D.radial_to_xy(r=new_speed)  # random direction

    def get_kinetic_energy(self):
        # Three different ways that the kinetic energy of this molecule can be changed:
        # 1. Explicitly through set_kinetic_energy()
        # 2. By modifying the velocity of the molecule through set_velocity()
        # 3. Implicitly by pymunk integrating force impulses on the molecule
        # For these reasons, must be careful velocity has not changed if optimizing this calculation by caching a value
        speed = self.get_speed()
        return 0.5 * self.get_mass() * speed * speed

    def get_speed(self):
        return self.body.velocity.get_length()
        # return Kinetics2D.get_speed(*self.body.velocity)

    def get_velocity(self):
        return self.body.velocity

    def set_velocity(self, vel):
        self.body.velocity = (vel.x,vel.y)

    def get_mass(self):
        return self._mass

    def get_potential_energy(self, chemistry):

        """
        Return the energy required to form all of the molecule's bonds (therefore a negative quantity as bond formation releases energy)
        WARNING: assumes formation energy = energy of breaking (symmetrical)

        :rtype: float
        """

        return sum([chemistry.get_bond_energy(bond.GetBeginAtom(), bond.GetEndAtom(), end_bond_type=int(bond.GetBondType())) for bond in self.GetBonds()])

    def get_internal_energy(self):
        return self._internal_energy

    def set_internal_energy(self, value):
        if value < 0:
            raise ValueError
        self._internal_energy = value

    def set_kinetic_energy(self, value):
        if value < 0:
            raise ValueError
        self._kinetic_energy = value

    def get_kinetic_energy(self):
        return self._kinetic_energy

    def same_component(self, idx1, idx2):
        # Components is a list of sets of atom indexes in molecule - same set = same component
        if self._components is None:  # Catch-all - components should have been set in init() and multiple components initialized in previous call to combine_molecules
            self._components = self._get_strongly_connected_components()

        for component in self._components:
            if idx1 in component and idx2 in component:
                return True
        return False

    def get_total_formal_charge(self):
        return sum([i.GetFormalCharge() for i in self.GetAtoms()])

    def _assign_formal_charge(self):
        """
        OpenEye Charge Model - http://www.eyesopen.com/docs/toolkits/current/html/OEChem_TK-python/valence.html#subsection-valence-openeye-charge
        The OpenEye formal charge model assigns formal charges to elements based upon their total valence.
        In OEChem, this functionality is invoked by the OEAssignFormalCharges function.
        If the formal charge on an atom is non-zero, it is left unchanged.

        Hydrogen: If the valence isn't one, the formal charge is +1.
        Carbon: If the valence is three, the formal charge is +1 if the atom has a polar neighbor, i.e. N, O or S, and formal charge -1 otherwise.
        Nitrogen: If the valence is two, the formal charge is -1, and if the valence is four the formal charge is +1.
        Oxygen: If the valence is one, the formal charge is -1, and if the valence is three the formal charge is +1.
        """

        for i in self.GetAtoms():

            valence = i.GetDegree()
            formal_charge = 0

            if i.GetAtomicNum() == 1:  # H
                if valence != 1:
                    formal_charge = 1

            if i.GetAtomicNum() == 6:  # C
                if valence == 3:
                    formal_charge = -1
                    for neighbor in i.GetNeighbors():
                        if neighbor.GetAtomicNum() == 7 or neighbor.GetAtomicNum() == 8:
                            formal_charge = 1

            if i.GetAtomicNum() == 7:  # N
                if valence == 2:
                    formal_charge = -1
                if valence == 3:
                    formal_charge = 1

            if i.GetAtomicNum() == 8:  # O
                if valence == 1:
                    formal_charge = -1
                elif valence == 3:
                    formal_charge = 1

            i.SetFormalCharge(formal_charge)

    def _get_strongly_connected_components(self):
        logging.info("Call to _get_strongly_connected_components")
        g = nx.Graph()
        for bond in self.GetBonds():
            g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        connected_components = list(nx.connected_components(g))
        # Add single atoms as independent components
        for idx in range(self.GetNumAtoms()):
            if len(self.GetAtomWithIdx(idx).GetBonds()) == 0:
                connected_components.append([idx])

        return connected_components

    def __deepcopy__(self, memo):
        # don't mess with current structure - just leave it exactly as it is
        return Molecule(self, internal_energy=self.get_internal_energy(), kinetic_energy=self.get_kinetic_energy(), canonize=False)

    def __str__(self):
        return Chem.MolToSmiles(self)
