import itertools
import Levenshtein


def uniform_weighting(reaction):

    """
    Uniform selection.

    :param reaction: Reaction
    :return: float
    """

    return 1.0


def replicant_weighting(reaction):

    """
    Selection pressure towards A+X -> 2A + Y

    Heuristic to weight reaction options biased towards:

    - 3 or more products
    - One product contains a section of a reactant, the more the better
    - Longer A's are better

    :param reaction: Reaction
    :return: float
    """

    best_weight = 1 # default is a low but selectable weight

    # Only consider reactions with 3 products
    # This also greatly reduces unproductive null "reactions" where the products == reactants
    if len(reaction.get_products()) > 2:

        # Similarity between reactants and products
        for reactant in reaction.get_reactants():

            # Find best combination of products - 2 products with the smallest combined distance from this reactant
            for combination in itertools.combinations(reaction.get_products(), 2):
                weight = 2*len(reactant.get_symbol())  # start with length of A
                weight -= sum([Levenshtein.distance(reactant.get_symbol(), product.get_symbol()) for product in combination])
                if weight > best_weight:
                    best_weight = weight

    return best_weight


def least_energy_weighting(reaction):

    """
    Weight by the least energy strategy.

    :param reaction: Reaction
    :return: float
    """

    if reaction.product_value < 0:
        return -reaction.product_value
    else:
        return max(0, reaction.reactant_value - reaction.product_value)
