class IReactantSelection(object):

    def get_reactants(self):

        """
        Idempotent method to identify the next reactants from the reactor population.

        :return: Reaction, for the reactants only, or None if no reactants can be found.
        """

        pass

    def react(self, reaction):

        """
        Action the provided Reaction by replacing the reactants in the population with the products.
        Raises ValueError if any of the IElements provided as reactants doesn't match an item in the population.

        :param reaction: Reaction
        :return: StateRecord

        """

        pass