from i_product_selection import IProductSelection
from kinetic_molecule import KineticMolecule


class LeastEnergy(IProductSelection):

    """
    Select the reaction which has the least energy required to transform reactants into products.


    """

    @classmethod
    def get_product(cls, reactions):
        for reaction in reactions:
            pass
    return reaction

    @staticmethod
    def weight_option(reaction_energy, energy_delta):
        if energy_delta < 0:
            return -energy_delta
        else:
            return max(0, reaction_energy - energy_delta)
