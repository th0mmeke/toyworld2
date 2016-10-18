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


def weight(reaction):
    if reaction.product_value < 0:
        return -reaction.product_value
    else:
        return max(0, reaction.reactant_value - reaction.product_value)


def least_energy_product_selection(reactions):
    if len(reactions) == 0:
        return []

    return choice(elements=reactions, weights=[weight(reaction) for reaction in reactions])
