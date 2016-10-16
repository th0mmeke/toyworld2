import importlib
import logging
import random

import pymunk as pm
from rdkit.Chem import AllChem as Chem
from reactant_selection import ReactantSelection


class SpatialReactor(ReactantSelection):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    With KineticMolecules of size 1, ratio of vessel size:molecule size = 10^6:1
    """

    reaction_vessel_size = 500  # -500->500
    REACTANT = 1
    PRODUCT = 2
    WALL = 3
    bodies = {}  # dictionary body:mol
    space = pm.Space()

    reactant_list = []
    lookup = {}

    def __init__(self, population):

        super(SpatialReactor, self).__init__(population)

        SpatialReactor.space.gravity = pm.Vec2d(0, 0)

        wall_thickness = 100  # nice and thick so that fast moving molecules don't tunnel straight through
        wall_end_point = SpatialReactor.reaction_vessel_size + wall_thickness
        self._walls = [pm.Segment(SpatialReactor.space.static_body, (-wall_end_point, -wall_end_point), (-wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(SpatialReactor.space.static_body, (-wall_end_point, wall_end_point), (wall_end_point, wall_end_point), wall_thickness),
                       pm.Segment(SpatialReactor.space.static_body, (wall_end_point, wall_end_point), (wall_end_point, -wall_end_point), wall_thickness),
                       pm.Segment(SpatialReactor.space.static_body, (-wall_end_point, -wall_end_point), (wall_end_point, -wall_end_point), wall_thickness)
                       ]
        for wall in self._walls:  # can't set these in the Segment constructor
            wall.collision_type = SpatialReactor.WALL
            wall.elasticity = 0.9999
            wall.friction = 0

        h = SpatialReactor.space.add_collision_handler(SpatialReactor.REACTANT, SpatialReactor.REACTANT)
        h.begin = SpatialReactor._begin_handler
        h.separate = SpatialReactor._end_handler

        locations = [[random.uniform(-SpatialReactor.reaction_vessel_size, SpatialReactor.reaction_vessel_size) for i in range(2)] for mol in population]

        # Add them into the main space using the new locations
        base_molecule_radius = SpatialReactor.reaction_vessel_size/100

        for molecule, (x, y) in zip(population, locations):
            mol = Chem.MolFromSmiles(molecule.get_symbol())
            mass = sum([atom.GetMass() for atom in mol.GetAtoms()])

            inertia = pm.moment_for_circle(mass, 0, base_molecule_radius, (0, 0))
            body = pm.Body(mass, inertia)
            body.position = x, y

            shape = pm.Circle(body, base_molecule_radius, (0, 0))
            shape.elasticity = 0.999  # required for standard collision handler to do a 'perfect' bounce
            shape.friction = 0.0
            shape.collision_type = SpatialReactor.REACTANT  # Mark molecule as potential reactant

            body.velocity=pm.Vec2d(math.sqrt(2.0 * ke / mass), 0)
            body.velocity.angle = random.uniform(-math.pi, math.pi)

            SpatialReactor.space.add(shape)
            SpatialReactor.space.add(body)

            SpatialReactor.lookup[shape] = molecule


    def get_reactants(self):

        while len(SpatialReactor.reactant_list) == 0:  # reactant_list maintained by _begin_handler()
            SpatialReactor.space.step(1)  # trigger _begin_handler on collision

        return SpatialReactor.reactant_list.pop(0)

    def react(self, reaction):
        # Find middle point of reactant bodies
        # Remove reactant bodies+shapes
        # Add in product bodies to middle point of reaction
        # Set product velocities to preserve momentum

    @classmethod
    def _end_handler(cls, arbiter, space, data):
        '''Called when two molecules separate. Mark them as potential reactants.
        :param arbiter:
        :param space:
        :param data:
        :return:
        '''

        for shape in arbiter.shapes:
            shape.collision_type = SpatialReactor.REACTANT
        return False

    @classmethod
    def _begin_handler(cls, arbiter, space, data):

        try:
            molecules = [SpatialReactor.lookup[shape] for shape in arbiter.shapes]
        except KeyError:
            return False  # one or more reactants were in collision with another molecule previously in this timestep

        SpatialReactor.reactant_list.append(molecules)


