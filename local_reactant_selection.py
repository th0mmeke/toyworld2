import random
from i_reactant_selection import IReactantSelection
from reaction import Reaction
import pymunk as pm


class LocalReactantSelection(IReactantSelection):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    """

    REACTION_VESSEL_SIZE = 500  # -500->500
    FOOD_SET = pm.Vec2d(REACTION_VESSEL_SIZE+1, REACTION_VESSEL_SIZE+1)
    BASE_MOLECULE_RADIUS = REACTION_VESSEL_SIZE / 500

    def __init__(self, population):

        """
        Population represents the Food Set of ChemMolecules. Only important thing is the relative proportions of the elements as this determines the
        the probability of a particular food set molecule taking part in a reaction.

        :param population: [ChemMolecule]
        """

        if not isinstance(population, list):
            raise TypeError

        self.population = {x: LocalReactantSelection.FOOD_SET for x in population}  # {ChemMolecule: (x,y)}

    def get_reactants(self):
        if len(self.population) > 1:
            reactants = random.sample(self.population.keys(), 2)
        else:
            reactants = self.population.keys()
        assert type(reactants) == list
        return Reaction(reactants=reactants)

    def react(self, reaction):

        """
        :param reaction: Reaction
        """

        locations = []
        for x in reaction.get_reactants():
            if x not in self.population:
                raise ValueError
            if self.population[x] != LocalReactantSelection.FOOD_SET:
                locations.append(self.population[x])
                del self.population[x]

        if len(locations) == 0:
            location = pm.Vec2d([random.uniform(-LocalReactantSelection.REACTION_VESSEL_SIZE, LocalReactantSelection.REACTION_VESSEL_SIZE) for i in range(2)])
        else:
            location = sum(locations)/len(locations)
        for x in reaction.get_products():
            self.population[x] = location

        return reaction.as_dict()

    def get_population(self):

        """
        Return the list of all elements in the population, including the Food Set elements.
        :return: [ChemMolecule]
        """

        return self.population.keys()
