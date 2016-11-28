import random
import logging
from i_reactant_selection import IReactantSelection
from reaction import Reaction
import pymunk as pm
from collections import Counter


class LocalReactantSelection(IReactantSelection):

    """
    n-dimensional structure of fixed size (ranging [-reaction_vessel_size,reaction_vessel_size] in each dimension).
    Molecules bounce off virtual walls.
    """

    REACTION_VESSEL_SIZE = 500  # -500->500
    FOOD_SET = pm.Vec2d(REACTION_VESSEL_SIZE+1, REACTION_VESSEL_SIZE+1)
    BASE_MOLECULE_RADIUS = REACTION_VESSEL_SIZE / 500
    VESICLE_SIZE_SQRD = BASE_MOLECULE_RADIUS * BASE_MOLECULE_RADIUS * 16

    def __init__(self, population, **kwargs):

        """
        Population represents the Food Set of ChemMolecules. Only important thing is the relative proportions of the elements as this determines the
        the probability of a particular food set molecule taking part in a reaction.

        :param population: [ChemMolecule]
        """

        if not isinstance(population, list):
            raise TypeError

        try:
            self.ke = kwargs['ke']
        except KeyError:
            self.ke = 100

        self.food_set = population
        self.population = {}  # {ChemMolecule: (x,y)}

    def _clamp_range(self, n):
        lower_bound = 0
        upper_bound = len(self.food_set)
        return min(upper_bound, max(lower_bound, int(n * len(self.food_set))))

    def update_environment(self, delta):
        if delta > 0:
            # add molecules to foodset
            n = min(self._clamp_range(delta), max(0, len(self.population)-len(self.food_set)))  # limit foodset size to twice population to prevent unrestricted growth scenarios
            self.food_set.extend(random.sample(self.food_set, n ))
        else:
            # remove molecules from foodset
            self.food_set = random.sample(self.food_set, self._clamp_range(1 - delta))
        logging.info("Food set size = {}".format(len(self.food_set)))

    def get_reactants(self):
        """
        Return the reactants for a reaction, essentially through uniform selection.
        However to mimic the effects of a membrane, all non-foodset reactants must be within a given distance of each other.
        We do this by first selecting a reactant with uniform probability from all molecules. If it is a foodset molecule,
        we can select a second reactant without constraint in the same manner. If it is not however, our second reactant can
        either be a foodset molecule or a molecule from the same location.

        :return: Reaction
        """

        reactants = random.sample(self.get_population(), 2)

        if reactants[0] not in self.food_set and reactants[1] not in self.food_set:

            # Reactant is within vesicle - second reactant must then be either 1) foodset or 2) from same location
            if random.random() <= len(self.food_set) * 1.0 / (len(self.population) + len(self.food_set)):
                sample_population = self.food_set
                logging.info("Food...")
            else:
                sample_population = self.population.copy()
                del sample_population[reactants[0]]
                sample_population = [x for x in sample_population if sample_population[x].get_dist_sqrd(self.population[reactants[0]]) <= self.VESICLE_SIZE_SQRD]
                logging.info("Vesicle of size {}".format(len(sample_population)))
                if len(sample_population) < 1:
                    # vesicle is tiny, and only holds reactant[0], so revert to food_set
                    sample_population = self.food_set

            reactants[1] = random.sample(sample_population, 1)[0]

        assert type(reactants) == list
        return Reaction(reactants=reactants, reactant_value=self.ke)

    def react(self, reaction):

        """
        Foodset molecules are treated as inexhaustible. The actual number of foodset molecules of a particular species
        is only important for get_reactants(), which selects reactants in proportion to their abundance.
        All products of a reaction are placed into the reaction vessel at the reaction location, with zero KE.
        This is a simple way to approximate the effects of a membrane - products stay in close proximity to each other.

        :param reaction: Reaction
        :return: Reaction
        """

        locations = []
        for x in reaction.get_reactants():
            if x in self.population:  # not in foodset
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

        return self.food_set + self.population.keys()
