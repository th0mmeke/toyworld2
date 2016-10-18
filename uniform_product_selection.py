import random


def product_selection(reactions):
    if reactions is None or len(reactions) == 0:
        return None

    return random.sample(reactions, 1)[0]
