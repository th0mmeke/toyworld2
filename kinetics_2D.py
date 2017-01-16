import logging
import math
import random

import pymunk as pm
from ulps import Ulps


class Kinetics2D(object):

    @classmethod
    def get_ke(cls, m, x, y):
        return 0.5 * m * (x * x + y * y)

    @classmethod
    def get_cm_energy(cls, mols):
        """
        Return KE of Centre of Mass: _ke = 1/2mv^2, where mv for the centre of mass = sum (mi * vi) for all particles i

        :param mols: [pm.Body]
        :rtype: float
        """

        total_mass = sum([mol.mass for mol in mols])
        return cls.get_ke(total_mass, *cls.get_cm_velocity(mols))

    @classmethod
    def get_cm_velocity(cls, mols):

        """
        Return the momentum (mdx,mdy) of the centre of mass for these particles
        :param mols: [pm.Body]
        :rtype: pm.Vec2d
        """

        cm_momentum = pm.Vec2d(0, 0)
        total_mass = sum([mol.mass for mol in mols])
        for mol in mols:
            cm_momentum += mol.velocity * mol.mass
        cm_velocity = cm_momentum / total_mass
        return cm_velocity

    @classmethod
    def inelastic_collision(cls, reactants, out_mass):
        """
        Determine velocities of product molecules following a collision of reactant molecules, for between one and three product molecules.

        Model as a collision, followed by an explosion, meaning that the total momentum of the system is conserved - if two particles, each has equal and opposite momentum in CoM frame
        Assume an impulse, or force, splitting the particles apart, acting equally on each particle
        Then impulse J = mv2-mv1 and so momentum change will be the same for all particles
        Implies that for two particles, equal and opposite mv in CoM frame, and for three particles, mv arranged in equilateral triangle

        Post-conditions:
        1. Sum in_mass = Sum out_mass - although #in_molecules ne #out_molecules
        2. Vector speed and direction of CoM remains constant
        3. in_KE + in_PE + in_IE = Sum out_KE + out_PE + out_IE or in_KE - delta_KE = out_KE

        :param reactants: [pm.Body]. Reactants - must have total KE > 0
        :param out_mass:[float]. Masses of the products of the reaction - must be 1, 2 or 3 products only
        :rtype: [pm.Vec2d]
        """

        def total_mv(mv):
            totals = [0, 0]
            for mv_ in mv:
                for dim in range(len(totals)):
                    totals[dim] += mv_[dim]
            return totals

        if len(out_mass) < 1 or len(out_mass) > 3:
            raise ValueError()

        in_v = [mol.velocity for mol in reactants]
        in_mass = [mol.mass for mol in reactants]
        in_mv = [[m * v_ for v_ in v] for m, v in zip(in_mass, in_v)]
        in_ke = sum([mol.kinetic_energy for mol in reactants])

        assert Ulps.almost_equal(sum(in_mass), sum(out_mass))

        # Velocity of centre of mass after collision
        # Momentums add to zero in the CoM frame
        cm_in_v = cls.get_cm_velocity(reactants)

        # Bound energy_of_collision to above zero (rounding errors for small values)
        # consistent sense with that in discover_reaction - final_PE = initial_PE + energy_delta => final_KE = initial_KE - energy_delta
        energy_of_collision = max(0, in_ke - cls.get_cm_energy(reactants))
        if energy_of_collision <= 0:
            raise ValueError

        if len(out_mass) == 1:
            # One out particle is stationary in out_CoM frame
            out_v_in_com_frame = [pm.Vec2d(0, 0)]

        elif len(out_mass) == 2:
            ke_in_cm_frame = random.uniform(0, energy_of_collision)
            mv = math.sqrt((2.0 * ke_in_cm_frame * out_mass[0] * out_mass[1]) / (out_mass[0] + out_mass[1]))
            v1 = pm.Vec2d([mv,0])
            v2 = pm.Vec2d([mv,0])
            v1.angle = cm_in_v.angle + math.pi * 0.5
            v2.angle = cm_in_v.angle + math.pi * 1.5
            out_v_in_com_frame = [v1, v2]

        elif len(out_mass) == 3:
            # Sum of vector momentums = 0, and in centre of momentum frame arranged as equilateral triangle, side mv
            # Must then convert to velocities by dividing by particle mass, which means no longer equilateral...but unimportant, as only needed equilateral to initially arrange
            ke_in_cm_frame = random.uniform(0, energy_of_collision)  # The energy of the collision - over and above the energy of the centre of mass, which is invariant
            mv = math.sqrt((2.0 * ke_in_cm_frame * out_mass[0] * out_mass[1] * out_mass[2]) / (out_mass[0] * out_mass[1] + out_mass[1] * out_mass[2] + out_mass[0] * out_mass[2]))
            v1 = pm.Vec2d([mv,0])
            v2 = pm.Vec2d([mv,0])
            v3 = pm.Vec2d([mv,0])
            v1.angle = cm_in_v.angle + math.pi / 3.0
            v2.angle = cm_in_v.angle - math.pi / 3.0
            v3.angle = cm_in_v.angle + math.pi
            out_v_in_com_frame = [v1, v2, v3]

        # Now convert from momentums to velocities by scaling by 1/mass
        out_v_in_com_frame = [particle_mv/mass for particle_mv, mass in zip(out_v_in_com_frame, out_mass)]
        # Finally convert back from CoM frame to lab frame
        out_v = [v+cm_in_v for v in out_v_in_com_frame]

        #########################
        # Confirm post-conditions
        # 1. Mass
        assert Ulps.almost_equal(sum(in_mass), sum(out_mass))  # should be safe

        # 2. Momentum

        out_mv = [[m * v_ for v_ in v] for m, v in zip(out_mass, out_v)]
        in_mv_total = total_mv(in_mv)
        out_mv_total = total_mv(out_mv)
        logging.debug("IN MV = {}, OUT MV = {}".format(in_mv_total, out_mv_total))
        for in_, out_ in zip(in_mv_total, out_mv_total):
            if not Ulps.almost_equal(in_, out_):
                raise ValueError

        # 3. Energy
        out_ke = sum([cls.get_ke(m, *v) for m, v in zip(out_mass, out_v)])
        if not Ulps.almost_equal(in_ke, out_ke, max_diff=3000):
            print(in_v, in_mass, out_v, out_mass)
            raise ValueError  # if a calculation problem with the velocities, throw it back to the caller

        return out_v
