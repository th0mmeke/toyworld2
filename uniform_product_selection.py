import random


def product_selection(reactions):
    if len(reactions) == 0:
        return []

    return random.sample(reactions, 1)[0]
