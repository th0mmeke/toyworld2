import random
import math
import logging
from collections import Counter

import pymunk as pm
from rdkit.Chem import AllChem as Chem
from i_reactant_selection import IReactantSelection
from kinetics_2D import Kinetics2D
from reaction import Reaction
from molecule import Molecule


class SpatialReactantSelection(IReactantSelection):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    """

    REACTION_VESSEL_SIZE = 500  # -500->500
    REACTANT = 1
    PRODUCT = 2
    WALL = 3
    BASE_MOLECULE_RADIUS = REACTION_VESSEL_SIZE / 50

    def __init__(self, population, **kwargs):

        """
        Population consists of Molecules with mass property.

        :param population: [Molecule with mass]
        :param kwargs: Standard python kwargs dict; expect a 'ke' property with the value of the initial KE for the reactor
        """

        if not isinstance(population, list):
            raise TypeError

        self.bodies = {}  # dictionary body:mol
        self.current_reactions = []
        self.reactant_list = []

        self.space = pm.Space()
        self.space.gravity = pm.Vec2d(0, 0)

        wall_thickness = 10000  # nice and thick so that fast moving molecules don't tunnel straight through
        wall_end_point = SpatialReactantSelection.REACTION_VESSEL_SIZE + wall_thickness
        self._walls = [pm.Segment(self.space.static_body, (-wall_end_point, -wall_end_point), (-wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (-wall_end_point, wall_end_point), (wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (wall_end_point, wall_end_point), (wall_end_point, -wall_end_point), wall_thickness),
                       pm.Segment(self.space.static_body, (-wall_end_point, -wall_end_point), (wall_end_point, -wall_end_point), wall_thickness)]
        for wall in self._walls:  # can't set these in the Segment constructor
            wall.collision_type = SpatialReactantSelection.WALL
            wall.elasticity = 0.999  # using 1.0 is not recommended by pymunk
            wall.friction = 0
            self.space.add(wall)

        h = self.space.add_collision_handler(SpatialReactantSelection.REACTANT, SpatialReactantSelection.REACTANT)
        h.begin = self._begin_handler
        h.separate = self._end_handler

        locations = [pm.Vec2d([random.uniform(-SpatialReactantSelection.REACTION_VESSEL_SIZE, SpatialReactantSelection.REACTION_VESSEL_SIZE) for i in range(2)]) for mol in population]

        # Add them into the main space using the new locations

        for molecule, location in zip(population, locations):
            velocity = pm.Vec2d(math.sqrt(2.0 * kwargs['ke'] / molecule.mass), 0)
            velocity.angle = random.uniform(-math.pi, math.pi)

            self._add_molecule(molecule, location=location, velocity=velocity, collision_type=SpatialReactantSelection.REACTANT)

        logging.info("Initial population: {}".format(Counter([str(x) for x in self.get_population()])))
        self.step_size = self._calculate_step_size()

    def get_reactants(self):

        i = 0
        while len(self.reactant_list) == 0:  # reactant_list maintained by _begin_handler()
            i += 1
            self.space.step(self.step_size)  # trigger _begin_handler on collision
            if i > 10 and (len(self.reactant_list) == 0 or len(self.reactant_list) > 10):
                self.step_size = self._calculate_step_size()
                i = 0

        # Now weed out any reactions that involve a molecule that no longer exists because of a prior reaction
        molecules = self.bodies.values()  # {shape: Molecule}

        while len(self.reactant_list) > 0:
            r = self.reactant_list.pop(0)  # {Molecule:pm.Body}
            try:
                [molecules.index(reactant) for reactant in r.keys()]
            except ValueError:
                pass
            else:
                self.current_reactions.append(r)

                # Energy available for the reaction = KE of molecules - KE of centre of mass

                bodies = r.values()
                initial_ke = sum([body.kinetic_energy for body in bodies])
                reaction_energy = initial_ke - Kinetics2D.get_cm_energy(bodies)
                return Reaction(reactants=r.keys(), reactant_value=reaction_energy)

        return None

    def react(self, reaction):

        # Match reaction to one in current_reactions - shouldn't be a very long list at all, so looping is acceptable
        for i in range(len(self.current_reactions)):  # [{Molecule: pm.Body}], loop by index so can later delete easily
            if set(reaction.reactants) == set(self.current_reactions[i].keys()):
                r = self.current_reactions[i]
                break
        else:
            raise ValueError

        # Find middle point of reactant bodies
        reactant_bodies = r.values()
        midpoint = sum([b.position for b in reactant_bodies]) / len(reactant_bodies)  # Vec2d

        # Remove reactant bodies+shapes
        for body in reactant_bodies:
            for shape in body.shapes:
                del self.bodies[shape]  # Remove shape:Molecule from the lookup table
            self.space.remove(body.shapes)

        self.space.remove(reactant_bodies)

        # Add in product bodies to middle point of reaction
        product_masses = [sum([atom.GetMass() for atom in Chem.AddHs(Chem.MolFromSmiles(molecule.get_symbol())).GetAtoms()]) for molecule in reaction.products]
        out_v = Kinetics2D.inelastic_collision(reactant_bodies, product_masses)
        for molecule, velocity, mass in zip(reaction.products, out_v, product_masses):
            molecule.mass = mass
            self._add_molecule(molecule, location=midpoint, velocity=velocity, collision_type=SpatialReactantSelection.PRODUCT)

        del self.current_reactions[i]

    def get_population(self):
        return self.bodies.values()

    def _add_molecule(self, molecule, location, velocity, collision_type):

        """
        Add a single Molecule into the reactor.

        :param molecule: Molecule with mass property
        :param location: pm.Vec2d
        :param velocity: pm.Vec2d
        :param collision_type: Int
        :return:
        """

        assert abs(location.x) <= self.REACTION_VESSEL_SIZE and abs(location.y) <= self.REACTION_VESSEL_SIZE
        # If fail assertion, then molecules are either skipping over walls, or the walls aren't reflecting them, or the walls are in the wrong place

        inertia = pm.moment_for_circle(molecule.mass, 0, SpatialReactantSelection.BASE_MOLECULE_RADIUS, (0, 0))
        body = pm.Body(molecule.mass, inertia)
        body.position = location
        body.velocity = velocity

        shape = pm.Circle(body, SpatialReactantSelection.BASE_MOLECULE_RADIUS, (0, 0))
        shape.elasticity = 0.999  # required for standard collision handler to do a 'perfect' bounce
        shape.friction = 0.0
        shape.collision_type = collision_type

        self.space.add(shape)
        self.space.add(body)

        self.bodies[shape] = molecule

    def _calculate_step_size(self):

        """
        Calculate an appropriate time step for the pymunk space updater based on the mean velocity of the bodies.
        Too small a step size will be slow, while too big a step size might miss collisions.
        Some missed collisions however are acceptable.

        :return: float
        """

        velocity_distribution = [shape.body.velocity.length for shape in self.bodies]
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
            shape.collision_type = SpatialReactantSelection.REACTANT
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
            reactants = {self.bodies[shape]: shape.body for shape in arbiter.shapes}
        except KeyError:
            return False  # one or more reactants were in collision with another molecule previously in this time-step

        self.reactant_list.append(reactants)

        return True


