class Reactor(object):

    def __init__(self, population, reactant_selection):

        if not isinstance(population, list):
            raise TypeError
        self.population = population
        self.reactant_selection = reactant_selection

    def get_reactants(self):

        """
        Idempotent method to identify the next reactants from the reactor population.
        Just calls the class reactant_selection method.
        :return: [IElement]
        """

        return self.reactant_selection(self.population)

    def react(self, reaction):
        '''
        Action the provided Reaction by replacing the reactants in the population with the products.
        Raises ValueError if any of the IElements provided as reactants doesn't match an item in the population.

        :param reaction: Reaction
        :return: StateRecord

        '''

        for x in reaction.reactants:
            self.population.remove(x)  # raises ValueError if this reactant is not in the population
        self.population.extend(reaction.products)

        return reaction.as_dict()

    def __str__(self):
        return ",".join([str(x) for x in self.population])
