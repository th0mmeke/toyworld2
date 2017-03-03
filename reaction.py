class Reaction(object):

    """
    Reactions record a transformation from initial reactant Molecules to the resulting product Molecules.
    Associated with the transformation is a 'weight', whose meaning is specific to the IChemistry which generated the Reaction.
    """

    def __init__(self, reactants, products=None, reactant_value=0, product_value=0):

        if not isinstance(reactants, list) or (products is not None and not isinstance(products, list)):
            raise TypeError

        if len(reactants) == 0:
            raise ValueError

        self.reactants = reactants
        self.products = products
        self.reactant_value = reactant_value
        self.product_value = product_value

    def get_reactants(self):
        return self.reactants

    def get_products(self):
        return self.products

    def get_reactant_value(self):
        return self.reactant_value

    def get_product_value(self):
        return self.product_value

    def as_dict(self):
        return {'reactants': {x.get_id(): x.get_symbol() for x in self.reactants},
                'products': {x.get_id(): x.get_symbol() for x in self.products},
                'reactant_value': self.reactant_value,
                'product_value': self.product_value}

    def __str__(self):
        return str(self.as_dict())
