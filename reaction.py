class Reaction(object):

    """
    Reactions record a transformation from initial reactant Molecules to the resulting product Molecules.
    Associated with the transformation is a 'weight', whose meaning is specific to the IChemistry which generated the Reaction.
    """

    def __init__(self, reactants, products, weight=1):

        if not isinstance(reactants, list) or not isinstance(products, list):
            raise TypeError

        self.reactants = reactants
        self.products = products
        self.weight = weight

    def get_reactants(self):
        return self.reactants

    def get_products(self):
        return self.products

    def get_weight(self):
        return self.weight

    def as_dict(self):
        return {'reactants': self.reactants, 'products': self.products, 'weight': self.weight}