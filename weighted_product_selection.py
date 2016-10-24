import random


def choice(elements, weights):

    """
    Return None if weights == 0

    :param elements:
    :param weights:
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


def product_selection(reactions):
    if reactions is None or len(reactions) == 0:
        return None

    return choice(elements=reactions, weights=[reaction.product_value for reaction in reactions])
