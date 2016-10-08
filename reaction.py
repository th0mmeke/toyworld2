class Reaction(object):

    def __init__(self, reactants, products):

        if not isinstance(reactants, list) or not isinstance(products, list):
            raise ValueError

        self.reactants = reactants
        self.products = products
