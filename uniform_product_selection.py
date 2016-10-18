import random


def product_selection(reactions):
    if len(reactions) == 0 or reactions is None:
        return None

    return random.sample(reactions, 1)[0]
