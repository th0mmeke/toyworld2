import random
import itertools


def choice(elements, weights):

    """
    Choose one item from a list of weighted alternatives.

    Return None if weights == 0

    :param elements: []
    :param weights: [Float] of same length as elements; total does not have to add to 1.0
    :return:
    """

    if len(elements) != len(weights):
        raise ValueError

    total = sum(weights)

    if total != 0:
        bound = random.random() * total

        cumulative_sum = 0
        for value, p in zip(elements, weights):
            cumulative_sum += p
            if cumulative_sum > bound:
                return value

    return None


def weighted_selection(reactions, weight_function, bias):
    if reactions is None or len(reactions) == 0:
        return None

    return choice(elements=reactions, weights=map(weight_function, reactions, itertools.repeat(bias, len(reactions))))
