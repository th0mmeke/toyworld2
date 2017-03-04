import itertools
import Levenshtein


def uniform_weighting(dummy1, dummy2):

    """
    Uniform selection.

    :param reaction: Reaction
    :return: float
    """

    return 1.0


def least_energy_weighting(reaction, dummy):

    """
    Weight by the least energy strategy.

    :param reaction: Reaction
    :return: float
    """

    if reaction.product_value < 0:
        return -reaction.product_value
    else:
        return max(0, reaction.reactant_value - reaction.product_value)


def biased_least_energy_weighting(reaction, bias):

    """
    Weight by the least energy strategy.

    :param reaction: Reaction
    :return: float
    """

    if reaction.product_value < 0:
        return -reaction.product_value
    else:
        return max(0, reaction.reactant_value - reaction.product_value + bias)
