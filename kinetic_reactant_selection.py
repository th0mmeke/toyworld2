import random
import math
import logging

import pymunk as pm
from i_reactant_selection import IReactantSelection
from kinetics_2D import Kinetics2D
from reaction import Reaction
from chem_molecule import ChemMolecule


class KineticReactantSelection(IReactantSelection):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    """

    REACTION_VESSEL_SIZE = 500  # -500->500
    REACTANT = 1
    PRODUCT = 2
    WALL = 3
    BASE_MOLECULE_RADIUS = REACTION_VESSEL_SIZE / 500

    def __init__(self, population, **kwargs):

        """
        Population consists of ChemMolecules with mass property.

        :param population: [ChemMolecule]
        :param kwargs: Standard python kwargs dict; expect a 'ke' property with the value of the initial KE for the reactor
        """

        if not isinstance(population, list):
            raise TypeError

        self.shape2mol = {}  # dictionary shape:mol
        self.mol2body = {}  # dictionary mol:body
        self.reactant_list = []

        try:
            ke = kwargs['ke']
        except KeyError:
            ke = 100

        self.ke = ke
        self.foodset = [mol.get_symbol() for mol in population]

        self.space = pm.Space()
        self.space.gravity = pm.Vec2d(0, 0)

        wall_thickness = 10000  # nice and thick so that fast moving molecules don't tunnel straight through
        wall_end_point = KineticReactantSelection.REACTION_VESSEL_SIZE + wall_thickness
        self._walls = [pm.Segment(self.space.static_body, (-wall_end_point, -wall_end_point), (-wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (-wall_end_point, wall_end_point), (wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (wall_end_point, wall_end_point), (wall_end_point, -wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (-wall_end_point, -wall_end_point), (wall_end_point, -wall_end_point), wall_thickness)]
        for wall in self._walls:  # can't set these in the Segment constructor
            wall.collision_type = KineticReactantSelection.WALL
            wall.elasticity = 0.999  # using 1.0 is not recommended by pymunk
            wall.friction = 0
            self.space.add(wall)

        h = self.space.add_collision_handler(KineticReactantSelection.REACTANT, KineticReactantSelection.REACTANT)
        h.begin = self._begin_handler
        h = self.space.add_collision_handler(KineticReactantSelection.PRODUCT, KineticReactantSelection.PRODUCT)
        h.separate = self._end_handler

        locations = [pm.Vec2d([random.uniform(-KineticReactantSelection.REACTION_VESSEL_SIZE, KineticReactantSelection.REACTION_VESSEL_SIZE) for i in range(2)]) for mol in population]

        # Add them into the main space using the new locations

        for molecule, location in zip(population, locations):
            velocity = pm.Vec2d(math.sqrt(2.0 * ke / molecule.mass), 0)
            velocity.angle = random.uniform(-math.pi, math.pi)

            self._add_molecule(molecule, location=location, velocity=velocity, collision_type=KineticReactantSelection.REACTANT)

        self.step_size = self._calculate_step_size()
        self._previous_value = None

    def update_environment(self, target, value):
        if self._previous_value is not None:
            delta = self._previous_value - int(value)
            if delta != 0:
                if target == "POPULATION":
                    self.adjust_population(delta)
                elif target == "KE":
                    self.adjust_ke(delta)

        self._previous_value = int(value)

    def adjust_ke(self, delta):
        for body in self.mol2body.values():
            # Can't set ke directly, so adjust velocity
            new_ke = self.ke + delta
            old_ke = body.kinetic_energy
            body.velocity = body.velocity * new_ke / old_ke
        # With changes in velocity, recalculate simulation rate
        self.step_size = self._calculate_step_size()

    def adjust_population(self, delta):
        """
        Adjust the foodset to the new size. If bigger, extend by a sample from the current foodset; if smaller,
        take a sample of the current set.
        By bounding the sample size to the size of the foodset we potentially vary from any input timeseries. But the effects
        are likely to be very minor.

        :param new_size: New food set size
        :return: [ChemMolecule]
        """

        if delta < 0:
            # Remove random molecules
            delta = min(len(self.shape2mol), -delta)

            kill_list = random.sample(self.shape2mol.values(), delta)
            self._remove_molecules(kill_list)
        else:
            # Add new molecules from foodset with ke=self.ke, with random direction from random wall
            delta = min(len(self.foodset), delta)
            add_list = [ChemMolecule(smiles) for smiles in random.sample(self.foodset, delta)]

            for molecule in add_list:

                velocity = pm.Vec2d(math.sqrt(2.0 * self.ke / molecule.mass), 0)
                wall_location = random.uniform(-KineticReactantSelection.REACTION_VESSEL_SIZE, KineticReactantSelection.REACTION_VESSEL_SIZE)
                velocity.angle = random.uniform(0, 2*math.pi)

                if 1.0/4.0*math.pi < velocity.angle <= 3.0/4.0*math.pi:  # bottom wall
                    location = pm.Vec2d(wall_location, -KineticReactantSelection.REACTION_VESSEL_SIZE)
                elif 3.0/4.0*math.pi < velocity.angle <= 5.0/4.0*math.pi:  # left wall
                    location = pm.Vec2d(-KineticReactantSelection.REACTION_VESSEL_SIZE, wall_location)
                elif 5.0/4.0*math.pi < velocity.angle <= 7.0/4.0*math.pi:  # top wall
                    location = pm.Vec2d(wall_location, KineticReactantSelection.REACTION_VESSEL_SIZE, )
                else:  # right wall
                    location = pm.Vec2d(KineticReactantSelection.REACTION_VESSEL_SIZE, wall_location)

                self._add_molecule(molecule, location=location, velocity=velocity, collision_type=KineticReactantSelection.REACTANT)

        logging.info("Population size = {}".format(len(self.get_population())))

    def get_reactants(self):

        i = 0
        while len(self.reactant_list) == 0:  # reactant_list maintained by _begin_handler()
            i += 1

            # Now step simulation
            self.space.step(self.step_size)  # trigger _begin_handler on collision
            if i > 30 and (len(self.reactant_list) == 0 or len(self.reactant_list) > 10):
                self.step_size = self._calculate_step_size()  # Throws ValueError if populaton size == 0
                if self.step_size > 1E4:
                    raise ValueError
                i = 0

        # Now weed out any reactions that involve a molecule that no longer exists because of a prior reaction
        molecules = self.shape2mol.values()  # {shape: Molecule}
        while len(self.reactant_list) > 0:
            r = self.reactant_list.pop(0)  # {Molecule:pm.Body}
            try:
                [molecules.index(reactant) for reactant in r.keys()]
            except ValueError:
                pass
            else:
                # Energy available for the reaction = KE of molecules - KE of centre of mass
                bodies = r.values()
                initial_ke = sum([body.kinetic_energy for body in bodies])
                reaction_energy = initial_ke - Kinetics2D.get_cm_energy(bodies)
                return Reaction(reactants=r.keys(), reactant_value=reaction_energy)

        return None

    def react(self, reaction):

        """
        :param reaction: Reaction
        """

        try:
            reactant_bodies = [self.mol2body[reactant] for reactant in reaction.get_reactants()]
        except KeyError:
            raise ValueError

        initial_bodies = len(self.space.bodies)
        initial_shapes = len(self.space.shapes)
        assert initial_bodies == initial_shapes - 4

        # Add in product bodies to middle point of reaction
        product_masses = [molecule.mass for molecule in reaction.get_products()]

        try:
            delta_ke, out_v = Kinetics2D.inelastic_collision(reactant_bodies, product_masses)  # might throw ValueError
        except ValueError:
            raise ValueError

        # Find middle point of reactant bodies
        midpoint = sum([b.position for b in reactant_bodies]) / len(reactant_bodies)  # Vec2d

        self._remove_molecules(reaction.get_reactants())

        assert initial_bodies - len(reactant_bodies) == len(self.space.bodies)
        assert initial_shapes - len(reactant_bodies) == len(self.space.shapes)

        for molecule, velocity in zip(reaction.get_products(), out_v):
            self._add_molecule(molecule, location=midpoint, velocity=velocity, collision_type=KineticReactantSelection.PRODUCT)

        assert initial_bodies - len(reactant_bodies) + len(reaction.get_products()) == len(self.space.bodies)
        assert initial_shapes - len(reactant_bodies) + len(reaction.get_products()) == len(self.space.shapes)

        return reaction.as_dict()

    def get_population(self):
        return self.shape2mol.values()

    @classmethod
    def _clamp_to_vessel(cls, location):
        """
        Limit to within the reaction vessel
        :param location: Vec2d
        :return: Vec2d
        """

        location.x = min(cls.REACTION_VESSEL_SIZE, max(-cls.REACTION_VESSEL_SIZE, location.x))
        location.y = min(cls.REACTION_VESSEL_SIZE, max(-cls.REACTION_VESSEL_SIZE, location.y))
        return location

    def _add_molecule(self, molecule, location, velocity, collision_type):

        """
        Add a single Molecule into the reactor.

        :param molecule: ChemMolecule
        :param location: pm.Vec2d
        :param velocity: pm.Vec2d
        :param collision_type: Int
        :return:
        """

        location = self._clamp_to_vessel(location)
        assert abs(location.x) <= self.REACTION_VESSEL_SIZE and abs(location.y) <= self.REACTION_VESSEL_SIZE
        # If fail assertion, then molecules are either skipping over walls, or the walls aren't reflecting them, or the walls are in the wrong place

        inertia = pm.moment_for_circle(molecule.mass, 0, KineticReactantSelection.BASE_MOLECULE_RADIUS, (0, 0))
        body = pm.Body(molecule.mass, inertia)
        body.position = location
        body.velocity = velocity

        shape = pm.Circle(body, KineticReactantSelection.BASE_MOLECULE_RADIUS, (0, 0))
        shape.elasticity = 0.999  # required for standard collision handler to do a 'perfect' bounce
        shape.friction = 0.0
        shape.collision_type = collision_type

        self.space.add(shape)
        self.space.add(body)

        self.shape2mol[shape] = molecule
        self.mol2body[molecule] = body

    def _remove_molecules(self, molecules):

        bodies = [self.mol2body[mol] for mol in molecules]

        # Remove reactant bodies+shapes
        for body in bodies:
            for shape in body.shapes:
                del self.shape2mol[shape]  # Remove shape:Molecule from the lookup table
            self.space.remove(body.shapes)

        self.space.remove(bodies)
        for mol in molecules:
            del self.mol2body[mol]

    def _calculate_step_size(self):

        """
        Calculate an appropriate time step for the pymunk space updater based on the mean velocity of the bodies.
        Too small a step size will be slow, while too big a step size might miss collisions.
        Some missed collisions however are acceptable.

        :return: float
        """

        if len(self.shape2mol) == 0:
            raise ValueError

        velocity_distribution = [shape.body.velocity.length for shape in self.shape2mol]
        return 1.0 * len(velocity_distribution) / sum(velocity_distribution)

    def _end_handler(self, arbiter, space, data):

        """
        Called when two molecules separate. Mark them as potential reactants.
        :param arbiter:
        :param space:
        :param data:
        :return:
        """

        for shape in arbiter.shapes:
            shape.collision_type = KineticReactantSelection.REACTANT
        return False

    def _begin_handler(self, arbiter, space, data):

        """
        Add {Molecule: pm.Body} for all reactants to reactant_list
        :param arbiter:
        :param space:
        :param data:
        :return:
        """

        try:
            reactants = {self.shape2mol[shape]: shape.body for shape in arbiter.shapes}
        except KeyError:
            return False  # one or more reactants were in collision with another molecule previously in this time-step

        self.reactant_list.append(reactants)

        return True


