class Reaction(object):

    def __init__(self, reactants, products):

        if not isinstance(reactants, list) or not isinstance(products, list):
            raise TypeError

        self.reactants = reactants
        self.products = products

    def as_dict(self):
        return {'reactants': self.reactants, 'products': self.products}