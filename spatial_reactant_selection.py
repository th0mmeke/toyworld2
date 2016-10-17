import random
import math

import pymunk as pm
from rdkit.Chem import AllChem as Chem
from reactant_selection import ReactantSelection
from kinetics_2D import Kinetics2D

from molecule import Molecule


class SpatialReactantSelection(ReactantSelection):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    With KineticMolecules of size 1, ratio of vessel size:molecule size = 10^6:1
    """

    REACTION_VESSEL_SIZE = 500  # -500->500
    REACTANT = 1
    PRODUCT = 2
    WALL = 3
    BASE_MOLECULE_RADIUS = REACTION_VESSEL_SIZE / 500
    bodies = {}  # dictionary body:mol
    space = pm.Space()

    reactant_list = []
    lookup = {}

    def __init__(self, population, **kwargs):

        super(SpatialReactantSelection, self).__init__(population)

        self.current_reactions = []

        SpatialReactantSelection.space.gravity = pm.Vec2d(0, 0)

        wall_thickness = 100  # nice and thick so that fast moving molecules don't tunnel straight through
        wall_end_point = SpatialReactantSelection.REACTION_VESSEL_SIZE + wall_thickness
        self._walls = [pm.Segment(SpatialReactantSelection.space.static_body, (-wall_end_point, -wall_end_point), (-wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(SpatialReactantSelection.space.static_body, (-wall_end_point, wall_end_point), (wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(SpatialReactantSelection.space.static_body, (wall_end_point, wall_end_point), (wall_end_point, -wall_end_point), wall_thickness),
                       pm.Segment(SpatialReactantSelection.space.static_body, (-wall_end_point, -wall_end_point), (wall_end_point, -wall_end_point), wall_thickness)
                       ]
        for wall in self._walls:  # can't set these in the Segment constructor
            wall.collision_type = SpatialReactantSelection.WALL
            wall.elasticity = 0.9999
            wall.friction = 0

        h = SpatialReactantSelection.space.add_collision_handler(SpatialReactantSelection.REACTANT, SpatialReactantSelection.REACTANT)
        h.begin = SpatialReactantSelection._begin_handler
        h.separate = SpatialReactantSelection._end_handler

        locations = [[random.uniform(-SpatialReactantSelection.REACTION_VESSEL_SIZE, SpatialReactantSelection.REACTION_VESSEL_SIZE) for i in range(2)] for mol in population]

        # Add them into the main space using the new locations

        for molecule, location in zip(population, locations):
            mol = Chem.MolFromSmiles(molecule.get_symbol())  # purely to calculate molecule mass
            mass = sum([atom.GetMass() for atom in mol.GetAtoms()])
            velocity = pm.Vec2d(math.sqrt(2.0 * kwargs['ke'] / mass), 0)
            velocity.angle = random.uniform(-math.pi, math.pi)

            self.add_molecule(molecule, mass, location, velocity, SpatialReactantSelection.REACTANT)

    def get_reactants(self):

        while len(SpatialReactantSelection.reactant_list) == 0:  # reactant_list maintained by _begin_handler()
            SpatialReactantSelection.space.step(0.1)  # trigger _begin_handler on collision
        if len(SpatialReactantSelection.reactant_list) > 10:
            print(len(SpatialReactantSelection.reactant_list))

        # Now weed out any reactions that involve a molecule that no longer exists because of a prior reaction
        molecules = SpatialReactantSelection.lookup.values()  # {shape: Molecule}

        while len(SpatialReactantSelection.reactant_list) > 0:
            r = SpatialReactantSelection.reactant_list.pop(0)  # {Molecule:pm.Body}
            try:
                idx = [molecules.index(reactant) for reactant in r.keys()]
            except ValueError:
                pass
            else:
                self.current_reactions.append(r)
                return r.keys()

        return []

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
        midpoint = sum([b.position for b in reactant_bodies])/len(reactant_bodies)  # Vec2d

        # Remove reactant bodies+shapes
        for body in reactant_bodies:
            for shape in body.shapes:
                del SpatialReactantSelection.lookup[shape]  # Remove shape:Molecule from the lookup table
            self.space.remove(body.shapes)

        self.space.remove(reactant_bodies)

        # Add in product bodies to middle point of reaction
        product_masses = [sum([atom.GetMass() for atom in Chem.MolFromSmiles(molecule.get_symbol()).GetAtoms()]) for molecule in reaction.products]
        out_v = Kinetics2D.inelastic_collision(reactant_bodies, product_masses)
        for molecule, velocity, mass in zip(reaction.products, out_v, product_masses):
            self.add_molecule(molecule, mass, midpoint, velocity, SpatialReactantSelection.PRODUCT)

        del self.current_reactions[i]

    @classmethod
    def get_population(cls):
        return SpatialReactantSelection.lookup.values()

    @classmethod
    def add_molecule(cls, molecule, mass, location, velocity, collision_type):

        """

        :param molecule: Molecule
        :param location: pm.Vec2d
        :param velocity: pm.Vec2d
        :param collision_type: Int
        :return:
        """

        inertia = pm.moment_for_circle(mass, 0, cls.BASE_MOLECULE_RADIUS, (0, 0))
        body = pm.Body(mass, inertia)
        body.position = location
        body.velocity = velocity

        shape = pm.Circle(body, cls.BASE_MOLECULE_RADIUS, (0, 0))
        shape.elasticity = 0.999  # required for standard collision handler to do a 'perfect' bounce
        shape.friction = 0.0
        shape.collision_type = collision_type

        SpatialReactantSelection.space.add(shape)
        SpatialReactantSelection.space.add(body)

        SpatialReactantSelection.lookup[shape] = molecule

    @classmethod
    def _end_handler(cls, arbiter, space, data):
        '''Called when two molecules separate. Mark them as potential reactants.
        :param arbiter:
        :param space:
        :param data:
        :return:
        '''

        for shape in arbiter.shapes:
            shape.collision_type = cls.REACTANT
        return False

    @classmethod
    def _begin_handler(cls, arbiter, space, data):

        """
        Add {Molecule: pm.Body} for all reactants to cls.reactant_list
        :param arbiter:
        :param space:
        :param data:
        :return:
        """

        try:
            reactants = {SpatialReactantSelection.lookup[shape]: shape.body for shape in arbiter.shapes}
        except KeyError:
            return False  # one or more reactants were in collision with another molecule previously in this timestep

        cls.reactant_list.append(reactants)

        return True


