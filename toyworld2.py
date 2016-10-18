from collections import Counter
import logging


class ToyWorld:

    def __init__(self, reactor, chemistry, product_selection):
        self.product_selection = product_selection
        self.chemistry = chemistry
        self.reactor = reactor

    def run(self, generations, state):

        logging.info("Initial population: {}".format(Counter([str(x) for x in self.reactor.get_population()])))

        for i in range(generations):
            partial_reaction = self.reactor.get_reactants()
            reactions = self.chemistry.enumerate(partial_reaction)
            if reactions is not None:
                reaction = self.product_selection(reactions)
                self.reactor.react(reaction)

                state.add(reaction.as_dict())

        logging.info("Final population: {}".format(Counter([str(x) for x in self.reactor.get_population()])))

        return state
